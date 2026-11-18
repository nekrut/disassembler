#!/usr/bin/env python3
"""Build the site-by-accession frequency matrix and heatmap from a cohort's VCFs.

The walker writes one VCF per accession. This reads them directly rather than merging
first: the matrix is what a merge is wanted for here, and building it from the per-accession
files avoids depending on bcftools and on a merge's INFO, which carries one sample's values
rather than the cohort's.

    python3 make_af_matrix.py <dir of *.vcf> [out-prefix]

Writes `<prefix>.af_matrix.tsv` and `<prefix>.af_matrix.png`. Only matplotlib is needed.
"""
from __future__ import annotations

import csv
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_cohort(vcf_dir: pathlib.Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    accessions = sorted(p.stem for p in vcf_dir.glob("*.vcf"))
    sites: dict[str, dict[str, float]] = {}
    for acc in accessions:
        for line in (vcf_dir / f"{acc}.vcf").read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            sites.setdefault(f"{f[1]}:{f[3]}>{f[4]}", {})[acc] = float(info["AF"])
    return accessions, sites


def write_tsv(path: pathlib.Path, accessions, sites, order) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["site"] + accessions + ["n_samples"])
        for key in order:
            row = sites[key]
            w.writerow(
                [key]
                + [f"{row[a]:.3f}" if a in row else "" for a in accessions]
                + [len(row)]
            )


def draw(path: pathlib.Path, accessions, sites, order) -> None:
    m = np.full((len(order), len(accessions)), np.nan)
    for i, key in enumerate(order):
        for j, acc in enumerate(accessions):
            if acc in sites[key]:
                m[i, j] = sites[key][acc]

    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("#ebebeb")
    fig, ax = plt.subplots(figsize=(1.35 * len(accessions) + 4, 0.42 * len(order) + 3))
    im = ax.imshow(m, cmap=cmap, vmin=0.10, vmax=0.40, aspect="auto")

    import os

    # Label each column by the part that differs, keeping at least four characters so the
    # accession stays recognisable when a cohort shares a long prefix.
    common = len(os.path.commonprefix(accessions))
    cut = min(common, max(0, min(len(a) for a in accessions) - 4))
    prefix = accessions[0][:cut]
    ax.set_xticks(range(len(accessions)))
    ax.set_xticklabels([a[cut:] for a in accessions], rotation=45, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel(f"accession ({prefix}..)" if prefix else "accession")
    ax.set_ylabel("Zika NS5 site (query coord)")
    ax.set_title(
        f"Zika NS5 intra-host SNP frequencies across a {len(accessions)}-sample cohort\n"
        "(cell = minor-allele ka:f ratio; grey = not called)"
    )
    for i in range(len(order)):
        for j in range(len(accessions)):
            if not np.isnan(m[i, j]):
                ax.text(
                    j, i, f"{m[i, j]:.2f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if m[i, j] > 0.26 else "black",
                )
    fig.colorbar(im, ax=ax, shrink=0.55, label="minor-allele frequency (ka:f ratio)")
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def main() -> int:
    if not 2 <= len(sys.argv) <= 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    vcf_dir = pathlib.Path(sys.argv[1])
    prefix = pathlib.Path(sys.argv[2]) if len(sys.argv) == 3 else vcf_dir / "zika_ns5"

    accessions, sites = read_cohort(vcf_dir)
    if not accessions:
        print(f"no VCFs in {vcf_dir}", file=sys.stderr)
        return 1
    order = sorted(sites, key=lambda k: int(k.split(":")[0]))

    write_tsv(prefix.with_suffix(".af_matrix.tsv"), accessions, sites, order)
    draw(prefix.with_suffix(".af_matrix.png"), accessions, sites, order)
    shared = sum(1 for k in order if len(sites[k]) >= 2)
    print(
        f"{len(order)} sites across {len(accessions)} accessions, "
        f"{shared} seen in two or more"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
