//! The run's private scratch directory.
//!
//! GraphAligner reads and writes files, so each accession in flight has a GFA, a GAF and a
//! one-record query FASTA in scratch. They used to go straight into a shared
//! `/dev/shm/logan-walker`. Created by one user with mode 0755, that directory failed every
//! other user's run with a bare EACCES; a run stopped with SIGTERM left its files there (12
//! were measured), holding RAM until reboot; and `/dev/shm` does not exist on macOS.
//!
//! Each run now makes its own directory, `logan-walker.<pid>.<suffix>` with mode 0700, under
//! `--scratch-dir` (default `/dev/shm` where it exists, else the system temp directory). It is
//! removed when the run ends, on SIGINT, SIGTERM or SIGHUP, and, if a SIGKILL or the OOM
//! killer left it behind, by the next run over the same base once that PID is gone.
//! `LOGAN_KEEP_SCRATCH` keeps it, for debugging.

use crate::error::{file_err, Error, Result};
use std::path::{Path, PathBuf};

const PREFIX: &str = "logan-walker.";

pub struct Scratch {
    path: PathBuf,
    keep: bool,
}

impl Scratch {
    /// Make this run's directory under `base`, or under the default base when `None`.
    pub fn create(base: Option<&Path>) -> Result<Scratch> {
        let base = base.map(Path::to_path_buf).unwrap_or_else(default_base);
        std::fs::create_dir_all(&base).map_err(file_err("create the scratch directory", &base))?;
        let swept = sweep_stale(&base);
        if swept > 0 {
            eprintln!(
                "[logan-walker] removed {swept} scratch director{} left in {} by runs that \
                 were killed",
                if swept == 1 { "y" } else { "ies" },
                base.display()
            );
        }
        let keep = std::env::var_os("LOGAN_KEEP_SCRATCH").is_some();
        let pid = std::process::id();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos())
            .unwrap_or(0);
        for attempt in 0..16u32 {
            let path = base.join(format!("{PREFIX}{pid}.{nanos:08x}{attempt:x}"));
            match make_private_dir(&path) {
                Ok(()) => return Ok(Scratch { path, keep }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => {
                    return Err(Error::Usage(format!(
                        "cannot create a scratch directory in {}: {e}. Pass --scratch-dir \
                         <DIR> to use another",
                        base.display()
                    )))
                }
            }
        }
        Err(Error::Usage(format!(
            "could not find a free scratch directory name in {}",
            base.display()
        )))
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// On SIGINT, SIGTERM or SIGHUP, remove the directory, then end the process as the
    /// signal would have, so the caller still sees "terminated by SIGTERM".
    pub fn remove_on_signal(&self) -> Result<()> {
        use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};
        let mut signals = signal_hook::iterator::Signals::new([SIGINT, SIGTERM, SIGHUP])?;
        let (path, keep) = (self.path.clone(), self.keep);
        std::thread::spawn(move || {
            if let Some(sig) = signals.forever().next() {
                if !keep {
                    remove(&path);
                }
                let _ = signal_hook::low_level::emulate_default_handler(sig);
                std::process::exit(128 + sig);
            }
        });
        Ok(())
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        if self.keep {
            eprintln!("[logan-walker] LOGAN_KEEP_SCRATCH: kept {}", self.path.display());
        } else {
            remove(&self.path);
        }
    }
}

/// `remove_dir_all`, tried twice: a worker can create a file in the directory while it is
/// being emptied, which makes the first attempt fail with "directory not empty".
fn remove(path: &Path) {
    if std::fs::remove_dir_all(path).is_err() {
        let _ = std::fs::remove_dir_all(path);
    }
}

/// tmpfs `/dev/shm` where it exists (Linux), otherwise the system temp directory (macOS).
fn default_base() -> PathBuf {
    let shm = Path::new("/dev/shm");
    if shm.is_dir() {
        shm.to_path_buf()
    } else {
        std::env::temp_dir()
    }
}

#[cfg(unix)]
fn make_private_dir(path: &Path) -> std::io::Result<()> {
    use std::os::unix::fs::DirBuilderExt;
    std::fs::DirBuilder::new().mode(0o700).create(path)
}

#[cfg(not(unix))]
fn make_private_dir(path: &Path) -> std::io::Result<()> {
    std::fs::create_dir(path)
}

/// Remove run directories in `base` that are ours and whose process no longer exists.
/// Another user's directory, or one whose PID is alive (even if reused by another program),
/// is left alone. Returns how many were removed.
#[cfg(unix)]
fn sweep_stale(base: &Path) -> usize {
    use std::os::unix::fs::MetadataExt;
    let Ok(entries) = std::fs::read_dir(base) else { return 0 };
    // SAFETY: getuid cannot fail and has no preconditions.
    let uid = unsafe { libc::getuid() };
    let mut removed = 0;
    for entry in entries.flatten() {
        let name = entry.file_name();
        let Some(pid) = name
            .to_str()
            .and_then(|n| n.strip_prefix(PREFIX))
            .and_then(|rest| rest.split('.').next())
            .and_then(|p| p.parse::<i32>().ok())
        else {
            continue;
        };
        let Ok(meta) = entry.metadata() else { continue };
        if !meta.is_dir() || meta.uid() != uid || pid as u32 == std::process::id() {
            continue;
        }
        if !process_exists(pid) {
            remove(&entry.path());
            removed += 1;
        }
    }
    removed
}

#[cfg(not(unix))]
fn sweep_stale(_base: &Path) -> usize {
    0
}

/// `kill(pid, 0)`: it fails with ESRCH only when no such process exists. EPERM means it
/// exists under another user.
#[cfg(unix)]
fn process_exists(pid: i32) -> bool {
    // SAFETY: signal 0 delivers nothing; it only checks that the process exists.
    let r = unsafe { libc::kill(pid, 0) };
    r == 0 || std::io::Error::last_os_error().raw_os_error() != Some(libc::ESRCH)
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    fn base(name: &str) -> PathBuf {
        let b = std::env::temp_dir().join(format!("logan_scratch_{}_{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&b);
        b
    }

    #[test]
    fn a_private_directory_removed_on_drop() {
        let b = base("drop");
        let s = Scratch::create(Some(&b)).unwrap();
        let p = s.path().to_path_buf();
        assert!(p.starts_with(&b) && p.is_dir());
        assert_eq!(std::fs::metadata(&p).unwrap().permissions().mode() & 0o777, 0o700);
        std::fs::write(p.join("x.sub.gfa"), b"H\n").unwrap();
        drop(s);
        assert!(!p.exists(), "the run's directory must be removed with its files");
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn runs_get_separate_directories() {
        let b = base("two");
        let (s1, s2) = (Scratch::create(Some(&b)).unwrap(), Scratch::create(Some(&b)).unwrap());
        assert_ne!(s1.path(), s2.path());
        drop((s1, s2));
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn a_dead_runs_directory_is_swept_and_a_live_ones_kept() {
        let b = base("sweep");
        std::fs::create_dir_all(&b).unwrap();
        // A PID that has just exited stands in for a run killed with SIGKILL.
        let mut child = std::process::Command::new("true").spawn().unwrap();
        let dead = child.id();
        child.wait().unwrap();
        let stale = b.join(format!("{PREFIX}{dead}.0"));
        std::fs::create_dir(&stale).unwrap();
        std::fs::write(stale.join("left.sub.gfa"), b"H\n").unwrap();
        let live = b.join(format!("{PREFIX}1.0")); // PID 1 always exists
        std::fs::create_dir(&live).unwrap();
        let other = b.join("not-ours");
        std::fs::create_dir(&other).unwrap();

        let s = Scratch::create(Some(&b)).unwrap();
        assert!(!stale.exists(), "the dead run's directory is removed");
        assert!(live.exists() && other.exists(), "only dead runs' directories are touched");
        drop(s);
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn an_unwritable_base_names_itself() {
        let b = base("ro");
        std::fs::create_dir_all(&b).unwrap();
        std::fs::set_permissions(&b, std::fs::Permissions::from_mode(0o500)).unwrap();
        let err = Scratch::create(Some(&b)).err().map(|e| e.to_string());
        std::fs::set_permissions(&b, std::fs::Permissions::from_mode(0o700)).unwrap();
        let _ = std::fs::remove_dir_all(&b);
        if unsafe { libc::geteuid() } == 0 {
            return; // root writes anywhere
        }
        let err = err.expect("an unwritable base must fail");
        assert!(err.contains(&b.display().to_string()) && err.contains("--scratch-dir"), "{err}");
    }
}
