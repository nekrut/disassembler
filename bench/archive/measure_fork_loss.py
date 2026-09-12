#!/usr/bin/env python3
"""Measure the haplotypes discarded at fork choice during path walking.

Reproduces every number in `work/issue_intra_accession_variance_lost.md`.

Usage:
    measure_fork_loss.py --accession SRR8209082 --query queries/query_hiv1_pol.fasta
    measure_fork_loss.py --accession ERR2028820 --query queries/pfcrt_only.fasta
    measure_fork_loss.py --all                       # both, as in the issue doc
    measure_fork_loss.py --all --json                # machine-readable

`logan_path_walker.py` deletes its GFA and GAF in a `finally` block (:327-330), so the
evidence this needs does not survive a normal run. Stages 1 to 5 here are copied from the
walker and keep those files; the decompressed unitigs are cached in --workdir, so a second
run of the same accession takes seconds.

Deliberately self-contained rather than importing logan_path_walker, for two reasons: the
walker pulls in Biopython for one `SeqIO.parse` call, and phasing is reimplemented here with
the orientation fix (see rust-walker/src/stages/phase.rs), which the walker does not have.
Because it is a copy, --self-check asserts it still reproduces the golden fixtures.
"""
import sys, os, re, json, gzip, argparse, subprocess, collections
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
K, MIN_ABUND, MIN_VAF, MIN_COV, HOPS = 31, 5.0, 0.10, 0.85, 5

CASES = [
    ("ERR2028820", ROOT / "queries/pfcrt_only.fasta"),
    ("SRR8209082", ROOT / "queries/query_hiv1_pol.fasta"),
]
# From bench/golden/. Guards the copied stages against drift.
GOLDEN_GRAPH = {
    "ERR2028820": {"n_unitigs": 1277362, "n_anchors": 76, "subgraph_nodes": 356},
    "SRR8209082": {"n_unitigs": 81534, "n_anchors": 4102, "subgraph_nodes": 26323},
}

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]

def find_bin(name: str) -> str:
    """Same resolution order as logan_path_walker.find_bin, plus the repo's conda env."""
    mamba = Path(os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")))
    for c in (mamba / "envs/codon/bin" / name, ROOT / ".condaenv-walker/bin" / name):
        if c.exists():
            return str(c)
    w = subprocess.run(["which", name], capture_output=True, text=True)
    return w.stdout.strip() if w.returncode == 0 else name

# --------------------------------------------------------------------------
# stages 1-5, copied from logan_path_walker.py but retaining the artifacts
# --------------------------------------------------------------------------

def stream_unitigs(acc, out_fa, timeout=900):
    if out_fa.exists() and out_fa.stat().st_size > 0:
        pass
    else:
        url = f"https://logan-pub.s3.amazonaws.com/u/{acc}/{acc}.unitigs.fa.zst"
        with open(out_fa, "wb") as fh:
            p1 = subprocess.Popen(["curl", "-s", "--fail", "--max-time", str(timeout), url],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p2 = subprocess.Popen(["zstd", "-dc"], stdin=p1.stdout, stdout=fh,
                                  stderr=subprocess.PIPE)
            p1.stdout.close()
            _, err = p2.communicate()
            p1.wait()
        if p2.returncode != 0 or not out_fa.exists() or out_fa.stat().st_size == 0:
            if out_fa.exists():
                out_fa.unlink()
            raise RuntimeError(f"stream failed for {acc}: {err.decode(errors='ignore')[:200]}")
    return sum(1 for l in open(out_fa, errors="ignore") if l.startswith(">"))

def find_anchors_blast(unitigs_fa, query_fa, wd, acc, min_pident=80.0, evalue="1e-10"):
    db = wd / f"{acc}_bdb"
    subprocess.run([find_bin("makeblastdb"), "-in", str(unitigs_fa), "-dbtype", "nucl",
                    "-out", str(db)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=True)
    r = subprocess.run([find_bin("blastn"), "-task", "megablast", "-query", str(query_fa),
                        "-db", str(db), "-perc_identity", str(min_pident), "-evalue", evalue,
                        "-max_target_seqs", "100000", "-outfmt", "6 sseqid pident length bitscore"],
                       capture_output=True, text=True, check=True)
    m = {l.split("\t")[0].split("_")[-1] for l in r.stdout.strip().split("\n") if l}
    for p in wd.glob(f"{acc}_bdb.*"):
        p.unlink()
    return sorted(m)

def extract_local_subgraph(unitigs_fa, anchor_ids, gfa_out, hops=HOPS, kmer_len=K):
    adj = collections.defaultdict(set)
    links = collections.defaultdict(list)
    abundances = {}
    with open(unitigs_fa, errors="ignore") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            tokens = line[1:].split()
            uid = tokens[0].split("_")[-1]
            for t in tokens[1:]:
                if t.startswith("ka:f:"):
                    try:
                        abundances[uid] = float(t.split(":")[-1])
                    except ValueError:
                        pass
                elif t.startswith("L:"):
                    parts = t.split(":")
                    if len(parts) == 4:
                        _, s1, tgt, s2 = parts
                        adj[uid].add(tgt)
                        adj[tgt].add(uid)
                        links[uid].append((s1, tgt, s2))
    local = set(anchor_ids)
    for _ in range(hops):
        frontier = set()
        for n in local:
            frontier.update(adj.get(n, set()))
        local.update(frontier)
    seqs, cur, keep, buf = {}, None, False, []
    with open(unitigs_fa, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].split("_")[-1]
                keep, buf = cur in local, []
            elif keep:
                buf.append(line.strip())
        if keep and cur:
            seqs[cur] = "".join(buf)
    overlap = f"{kmer_len - 1}M"
    with open(gfa_out, "w") as o:
        o.write("H\tVN:Z:1.0\n")
        for uid, s in seqs.items():
            o.write(f"S\t{uid}\t{s}\n")
        for uid in seqs:
            for s1, tgt, s2 in links.get(uid, []):
                if tgt in seqs:
                    o.write(f"L\t{uid}\t{s1}\t{tgt}\t{s2}\t{overlap}\n")
    return {"num_nodes": len(seqs), "abundances": abundances, "seqs": seqs, "links": links}

def run_graphaligner(gfa, query_fa, gaf, threads=12):
    subprocess.run([find_bin("GraphAligner"), "-g", str(gfa), "-f", str(query_fa),
                    "-a", str(gaf), "-x", "dbg", "-t", str(threads)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return [l.rstrip("\n").split("\t") for l in open(gaf) if l.strip()]

# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------

def oriented_adjacency(seqs, links):
    """Read `L:s1:t:s2` in u's header as the traversal (u,s1) -> (t,s2).

    This is the orientation fix. `phase_superbubbles` discards s1 and s2, which turns a
    chain into a fork and a bidirected edge into a reconvergence. Validated by
    validate_convention(): every step of the walked path must be an edge here.
    """
    out = collections.defaultdict(set)
    for u, edge_list in links.items():
        if u not in seqs:
            continue
        for s1, t, s2 in edge_list:
            if t in seqs:
                out[(u, s1)].add((t, s2))
    return out

def validate_convention(path, out):
    """Consecutive path steps that are real oriented edges, under this reading and the
    flipped-target alternative. The conclusions are only meaningful if this passes."""
    ok = sum(1 for i in range(len(path) - 1) if path[i + 1] in out[path[i]])
    flip = collections.defaultdict(set)
    for k, vs in out.items():
        for t, s in vs:
            flip[k].add((t, "-" if s == "+" else "+"))
    bad = sum(1 for i in range(len(path) - 1) if path[i + 1] in flip[path[i]])
    return ok, bad, len(path) - 1

def stitch(steps, seqs):
    parts = []
    for i, (orient, nid) in enumerate(steps):
        s = seqs[nid]
        if orient == "<":
            s = revcomp(s)
        parts.append(s if i == 0 else s[K - 1:])
    return "".join(parts)

def phase_oriented(seqs, abundances, out):
    """Bubbles, with orientation honoured. Mirrors rust-walker/src/stages/phase.rs."""
    flip = lambda x: (x[0], "-" if x[1] == "+" else "+")
    bubbles = []
    for fork, branches in out.items():
        if len(branches) < 2:
            continue
        hits = collections.Counter()
        for b in sorted(branches):
            for c in sorted(out.get(b, ())):
                hits[c] += 1
        reconverged = sorted(c for c, n in hits.items() if n >= 2 and c != fork)
        if not reconverged:
            continue
        merge = reconverged[0]
        if fork > flip(merge):          # report each bidirected bubble once
            continue
        total = sum(abundances.get(b[0], 0.0) for b in branches)
        if total <= 0:
            continue
        kept = [(b, abundances.get(b[0], 0.0), abundances.get(b[0], 0.0) / total)
                for b in sorted(branches)
                if abundances.get(b[0], 0.0) >= MIN_ABUND
                and abundances.get(b[0], 0.0) / total >= MIN_VAF]
        if len(kept) >= 2:
            bubbles.append({"fork": fork, "merge": merge, "branches": sorted(kept)})
    return bubbles

def analyse(acc, query_fa, wd, max_detour):
    wd.mkdir(parents=True, exist_ok=True)
    uni = wd / f"{acc}.unitigs.fa"
    gfa, gaf = wd / f"{acc}.sub.gfa", wd / f"{acc}.sub.gaf"

    n_unitigs = stream_unitigs(acc, uni)
    anchors = find_anchors_blast(uni, query_fa, wd, acc)
    sub = extract_local_subgraph(uni, anchors, gfa)
    rows = run_graphaligner(gfa, query_fa, gaf)
    seqs, ab, links = sub["seqs"], sub["abundances"], sub["links"]

    qseq = "".join(l.strip() for l in open(query_fa) if not l.startswith(">"))
    qlen = len(qseq)

    R = {"accession": acc, "query": query_fa.name, "query_len": qlen,
         "graph": {"n_unitigs": n_unitigs, "n_anchors": len(anchors),
                   "subgraph_nodes": sub["num_nodes"]},
         "gaf_rows": len(rows), "gaf_rows_discarded": max(0, len(rows) - 1)}
    if not rows:
        R["error"] = "no alignment"
        return R

    best = sorted(rows, key=lambda r: int(r[9]), reverse=True)[0]
    steps = re.findall(r"([<>])(\d+)", best[5])
    ps, pe = int(best[7]), int(best[8])
    full = stitch(steps, seqs)
    primary = full[ps:pe] if pe <= len(full) else full[ps:]
    path = [(n, "+" if o == ">" else "-") for o, n in steps]
    pos = {p: i for i, p in enumerate(path)}
    out = oriented_adjacency(seqs, links)

    ok, bad, total_steps = validate_convention(path, out)
    R["orientation_check"] = {"steps_that_are_edges": ok, "flipped_reading": bad,
                              "total": total_steps, "passes": ok == total_steps and bad == 0}
    R["primary"] = {"length": len(primary), "cov": round(len(primary) / qlen, 6)}
    R["path_nodes"] = len(path)
    R["forks_on_path"] = sum(1 for i in range(len(path) - 1) if len(out[path[i]]) >= 2)

    bubbles = phase_oriented(seqs, ab, out)
    R["bubbles_corrected"] = len(bubbles)
    R["bubbles_with_fork_on_path"] = sum(1 for b in bubbles if b["fork"] in pos)

    pa = [ab.get(n, 0.0) for _, n in steps]
    sa = [v for k, v in ab.items() if k in seqs]
    q = lambda x, p: sorted(x)[int(p * (len(x) - 1))]
    R["coverage"] = {
        "path_min": min(pa), "path_median": q(pa, .5), "path_max": max(pa),
        "subgraph_median": q(sa, .5),
        "path_below_min_abund": sum(1 for v in pa if v < MIN_ABUND),
        "path_below_10x": sum(1 for v in pa if v < 10),
    }

    # Every fork the path traverses: compare the branch taken with those not taken.
    shape = collections.Counter()
    cands, records = [], [(f"{acc}_primary", primary)]
    for i in range(len(path) - 1):
        fork = path[i]
        branches = out[fork]
        taken = path[i + 1]
        if len(branches) < 2 or taken not in branches:
            continue
        total = sum(ab.get(b[0], 0.0) for b in branches)
        if total <= 0:
            continue
        tk = ab.get(taken[0], 0.0)
        hits = collections.Counter()
        for b in sorted(branches):
            for c in sorted(out.get(b, ())):
                hits[c] += 1
        reconverges = any(n >= 2 and c != fork for c, n in hits.items())
        for alt in sorted(branches):
            if alt == taken:
                continue
            a_ab = ab.get(alt[0], 0.0)
            a_vaf = a_ab / total
            if not (a_ab >= MIN_ABUND and a_vaf >= MIN_VAF):
                continue
            tk_passes = tk >= MIN_ABUND and tk / total >= MIN_VAF
            shape["reported_as_bubble" if reconverges and tk_passes else
                  "bubble_shape_suppressed_by_vaf_gate" if reconverges else
                  "no_reconvergence_at_depth_1"] += 1
            # BFS back onto the walked path: if it rejoins, the detour is an allele
            dq = collections.deque([(alt, [alt])])
            seen, hit = {alt}, None
            while dq:
                node, trail = dq.popleft()
                if len(trail) > max_detour:
                    continue
                j = pos.get(node)
                if j is not None and j > i + 1:
                    hit = (j, trail)
                    break
                for nxt in sorted(out.get(node, ())):
                    if nxt not in seen:
                        seen.add(nxt)
                        dq.append((nxt, trail + [nxt]))
            if not hit:
                cands.append({"rejoins": False, "ka_taken": tk, "ka_alt": a_ab, "vaf_alt": a_vaf})
                continue
            j, detour = hit
            spliced = steps[:i + 1] + [(">" if x[1] == "+" else "<", x[0])
                                       for x in detour[:-1]] + steps[j:]
            try:
                af = stitch(spliced, seqs)
            except KeyError:
                continue
            aseq = af[ps:pe] if pe <= len(af) else af[ps:]
            nd = (sum(1 for x, y in zip(aseq, primary) if x != y)
                  if len(aseq) == len(primary) else -1)
            name = f"{acc}_alt_fork{fork[0]}{fork[1]}"
            if any(r[0] == name for r in records):
                continue
            records.append((name, aseq))
            cands.append({"rejoins": True, "name": name, "fork": f"{fork[0]}{fork[1]}",
                          "path_index": i, "rejoin_index": j, "detour_nodes": len(detour) - 1,
                          "ka_taken": tk, "ka_alt": a_ab, "vaf_alt": a_vaf,
                          "length": len(aseq), "nt_from_primary": nd})

    R["fork_shapes"] = dict(shape)
    R["credible_untaken_branches"] = len(cands)
    R["rejoin_the_path"] = sum(1 for c in cands if c["rejoins"])

    # Score every constructed haplotype against the reference, as the pipeline would.
    fa = wd / f"{acc}.alts.fasta"
    with open(fa, "w") as fh:
        for n, s in records:
            fh.write(f">{n}\n{s}\n")
    bl = subprocess.run([find_bin("blastn"), "-task", "megablast", "-query", str(fa),
                         "-subject", str(query_fa), "-outfmt",
                         "6 qseqid pident length"], capture_output=True, text=True)
    hit = {}
    for line in bl.stdout.strip().split("\n"):
        if line:
            p = line.split("\t")
            hit.setdefault(p[0], (float(p[1]), int(p[2])))
    ppid, paln = hit.get(f"{acc}_primary", (0.0, 0))
    R["primary"]["pident"] = ppid
    clean = []
    for c in cands:
        if not c["rejoins"]:
            continue
        pid, aln = hit.get(c["name"], (0.0, 0))
        c["pident"], c["blast_cov"] = pid, round(aln / qlen, 4)
        if c["blast_cov"] >= MIN_COV and 0 < c["nt_from_primary"] <= 50:
            clean.append(c)
    R["clean_alternatives"] = len(clean)
    # Rejoining candidates that did not survive the clean filter, with the reason, so the
    # cases that fail can be cited as readily as the ones that pass.
    R["rejected"] = [
        {**{k: c[k] for k in ("name", "fork", "nt_from_primary", "pident", "blast_cov")},
         "reason": ("below min_cov" if c["blast_cov"] < MIN_COV else
                    "identical to primary" if c["nt_from_primary"] == 0 else
                    "not a local swap (>50 nt from primary)")}
        for c in cands if c["rejoins"] and c not in clean
    ]
    if clean:
        ratios = sorted(c["ka_alt"] / c["ka_taken"] for c in clean if c["ka_taken"] > 0)
        R["findings"] = {
            "pident_primary": ppid,
            "pident_alt_min": min(c["pident"] for c in clean),
            "pident_alt_max": max(c["pident"] for c in clean),
            "alts_less_identical_to_reference": sum(1 for c in clean if c["pident"] < ppid),
            "alts_more_identical_to_reference": sum(1 for c in clean if c["pident"] > ppid),
            "discarded_better_supported": sum(1 for c in clean if c["ka_alt"] > c["ka_taken"]),
            "discarded_is_majority_allele": sum(1 for c in clean if c["vaf_alt"] >= 0.5),
            "discarded_over_90pct_coverage": sum(1 for c in clean if c["vaf_alt"] >= 0.9),
            "coverage_ratio_median": round(ratios[len(ratios) // 2], 1) if ratios else None,
            "coverage_ratio_max": round(max(ratios), 1) if ratios else None,
            "walked_under_10x_vs_discarded_over_10000x":
                sum(1 for c in clean if c["ka_taken"] < 10 and c["ka_alt"] >= 10000),
        }
    R["alternatives"] = clean
    R["alts_fasta"] = str(fa)
    return R

# --------------------------------------------------------------------------

def report(R):
    print(f"=== {R['accession']}  ({R['query']}, {R['query_len']} nt reference) ===")
    g = R["graph"]
    exp = GOLDEN_GRAPH.get(R["accession"])
    tag = ""
    if exp:
        tag = "  [matches golden]" if all(g[k] == exp[k] for k in exp) else "  [DIFFERS FROM GOLDEN]"
    print(f"  unitigs {g['n_unitigs']:,}  anchors {g['n_anchors']:,}  "
          f"subgraph {g['subgraph_nodes']:,} nodes{tag}")
    if "error" in R:
        print(f"  {R['error']}")
        return
    oc = R["orientation_check"]
    print(f"  orientation check: {oc['steps_that_are_edges']}/{oc['total']} path steps are "
          f"oriented edges (flipped reading: {oc['flipped_reading']}/{oc['total']}) "
          f"-> {'PASS' if oc['passes'] else 'FAIL'}")
    print(f"  GAF rows {R['gaf_rows']}, discarded by 'best = rows[0]': {R['gaf_rows_discarded']}")
    print(f"  primary {R['primary']['length']} nt, cov {R['primary']['cov']:.4f}, "
          f"pident {R['primary']['pident']:.3f}%")
    c = R["coverage"]
    print(f"  ka:f along the {R['path_nodes']} walked nodes: min {c['path_min']:.1f}, "
          f"median {c['path_median']:.1f}, max {c['path_max']:.1f} "
          f"(subgraph median {c['subgraph_median']:.1f})")
    print(f"  walked nodes below {MIN_ABUND}x: {c['path_below_min_abund']}, below 10x: "
          f"{c['path_below_10x']}")
    print(f"  bubbles after the orientation fix: {R['bubbles_corrected']}, "
          f"with fork on the walked path: {R['bubbles_with_fork_on_path']}")
    print(f"  forks on the walked path: {R['forks_on_path']}")
    print(f"  credible untaken branches: {R['credible_untaken_branches']}  "
          f"(shapes: {R['fork_shapes']})")
    print(f"    ... rejoin the path within the detour limit: {R['rejoin_the_path']}")
    print(f"    ... clean alternative haplotypes: {R['clean_alternatives']}")
    f = R.get("findings")
    if not f:
        print("  no clean alternative haplotypes; nothing is being discarded here")
        for c in R["rejected"]:
            print(f"    rejected: fork {c['fork']:>9}  cov {c['blast_cov']:.3f}  "
                  f"{c['nt_from_primary']} nt from primary  ({c['reason']})")
        return
    n = R["clean_alternatives"]
    print(f"\n  identity to the reference: primary {f['pident_primary']:.3f}%, "
          f"alternatives {f['pident_alt_min']:.3f}% to {f['pident_alt_max']:.3f}%")
    print(f"    alts LESS identical than the primary: {f['alts_less_identical_to_reference']}/{n}")
    print(f"    alts MORE identical:                  {f['alts_more_identical_to_reference']}/{n}")
    print(f"  discarded branch better supported than the one walked: "
          f"{f['discarded_better_supported']}/{n}")
    print(f"  discarded branch is the majority allele (VAF >= 0.5): "
          f"{f['discarded_is_majority_allele']}/{n}")
    print(f"  discarded branch carries >= 90% of local coverage:    "
          f"{f['discarded_over_90pct_coverage']}/{n}")
    print(f"  coverage ratio discarded:walked  median {f['coverage_ratio_median']}x, "
          f"max {f['coverage_ratio_max']}x")
    print(f"  walked < 10x while discarded >= 10,000x: "
          f"{f['walked_under_10x_vs_discarded_over_10000x']} cases")
    worst = sorted([c for c in R["alternatives"] if c["ka_taken"] > 0],
                   key=lambda c: c["ka_taken"])[:5]
    print(f"\n  worst cases by coverage of the branch actually walked:")
    for c in worst:
        print(f"    fork {c['fork']:>9} (path idx {c['path_index']:>5}): walked {c['ka_taken']:>7.1f}x, "
              f"discarded {c['ka_alt']:>9.1f}x (VAF {c['vaf_alt']:.3f}), "
              f"{c['nt_from_primary']} nt from primary, pident {c['pident']:.3f}%")
    print(f"\n  alternative haplotypes written to {R['alts_fasta']}")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", "-a")
    ap.add_argument("--query", "-q")
    ap.add_argument("--all", action="store_true", help="both accessions from the issue doc")
    ap.add_argument("--workdir", default="/tmp/measure_fork_loss",
                    help="caches decompressed unitigs; reruns are fast")
    ap.add_argument("--max-detour", type=int, default=8,
                    help="hops allowed for an alternative to rejoin the path (default 8)")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--self-check", action="store_true",
                    help="exit 1 unless the copied stages still reproduce the golden graph "
                         "shape and the orientation check passes")
    args = ap.parse_args()

    if args.all:
        cases = CASES
    elif args.accession and args.query:
        cases = [(args.accession, Path(args.query))]
    else:
        ap.error("need --accession with --query, or --all")

    results, bad = [], 0
    for acc, q in cases:
        R = analyse(acc, Path(q), Path(args.workdir) / acc, args.max_detour)
        results.append(R)
        exp = GOLDEN_GRAPH.get(acc)
        if exp and any(R["graph"][k] != exp[k] for k in exp):
            bad += 1
        if not R.get("orientation_check", {}).get("passes", False):
            bad += 1
        if not args.json:
            report(R)
            print()
    if args.json:
        print(json.dumps(results, indent=2))
    if args.self_check:
        print(f"self-check: {'FAIL' if bad else 'PASS'}", file=sys.stderr)
        return 1 if bad else 0
    return 0

sys.exit(main())
