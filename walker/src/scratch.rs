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
//!
//! The same signal also stops the run's GraphAligner children (WF-15). They are the only
//! processes the walker starts, they are what the scratch files are for, and without this
//! they outlived the run: see [`watch_child`].

use crate::error::{file_err, Error, Result};
use std::path::{Path, PathBuf};

/// Whether the run should keep its scratch files, `LOGAN_KEEP_SCRATCH`.
///
/// TW-08. Four sites read the variable directly and any value, even `0`, turned cleanup off.
/// A developer who follows `docs/walker.md` and exports it then saw four cleanup tests fail
/// on a correct build. `cfg!(test)` makes the unit tests ignore the environment; the CLI
/// tests clear it from the child's environment instead, since they run a real binary.
pub fn keep_scratch() -> bool {
    !cfg!(test) && std::env::var_os("LOGAN_KEEP_SCRATCH").is_some()
}

/// A scratch directory owned by one unit test, created now and removed when the test's guard
/// drops (TW-10).
///
/// Several test helpers built a `Config` around a path under the system temp directory and
/// never removed it. `write_gfa` creates the directory, the test removed only the GFA, and a
/// `cargo test` run left about seven empty `logan_*_test_*` directories behind. On a machine
/// where the temp directory is `/dev/shm` those are not even on disk.
///
/// The name carries a per-call counter as well as the PID, so two tests running in parallel
/// never share a directory and one guard dropping cannot delete the other's files.
#[cfg(test)]
pub struct TestDir(PathBuf);

#[cfg(test)]
impl TestDir {
    /// A fresh directory tagged with `tag`.
    ///
    /// It is created here. `extract::write_gfa` used to create a missing scratch directory
    /// itself, which is how these tests got away with a path that did not exist; it now
    /// refuses instead (WF-17), because the one thing that removes the directory mid-run is
    /// this run's own signal handler, and a worker that puts it back undoes the cleanup.
    /// Tests that want a scratch directory make their own, as a run does.
    pub fn new(tag: &str) -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static N: AtomicU64 = AtomicU64::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "logan_{tag}_test_{}_{n}",
            std::process::id()
        ));
        let _ = std::fs::create_dir_all(&path);
        Self(path)
    }
    pub fn path(&self) -> PathBuf {
        self.0.clone()
    }
}

#[cfg(test)]
impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// PIDs of the GraphAligner children this run has started and not yet reaped.
///
/// WF-15. The signal thread used to remove the scratch directory and exit, and nothing knew
/// about the aligners: after `kill <walker pid>` or `timeout N logan-walker`, up to `-w`
/// GraphAligner processes kept burning CPU and RAM, reparented to init, working on files that
/// had just been deleted. Ctrl-C in a terminal hid it, because the tty signals the whole
/// foreground process group, and so do runners that kill a process group or a cgroup; a bare
/// `kill` to the walker alone does not.
///
/// Each entry is also a **process group id**, because `align::run` spawns with
/// `process_group(0)`, which makes the child its own group leader. Signalling the group rather
/// than the process is what reaches a `--graphaligner` wrapper script's own children.
///
/// The one race left is a PID reused between the child being reaped and [`ChildGuard`]
/// dropping, a window of microseconds. Nothing cheap closes it, and the cost of losing is one
/// stray SIGTERM to an unrelated group of this user's, at a moment when this process is
/// already exiting.
static ALIGNERS: std::sync::Mutex<Vec<u32>> = std::sync::Mutex::new(Vec::new());

/// Register a spawned aligner, and unregister it when the returned guard drops.
///
/// The caller holds the guard until it has reaped the child, so the list only ever names
/// processes that are, as far as this process knows, still running.
pub fn watch_child(pid: u32) -> ChildGuard {
    lock_aligners().push(pid);
    ChildGuard(pid)
}

pub struct ChildGuard(u32);

impl Drop for ChildGuard {
    fn drop(&mut self) {
        lock_aligners().retain(|&p| p != self.0);
    }
}

/// The aligner list, taking a poisoned lock rather than panicking.
///
/// A panicking compute worker must not disarm the cleanup for every other worker's child.
fn lock_aligners() -> std::sync::MutexGuard<'static, Vec<u32>> {
    ALIGNERS.lock().unwrap_or_else(|e| e.into_inner())
}

/// SIGTERM every running aligner's process group.
///
/// SIGTERM, not SIGKILL: this is the signal a scheduler would have sent them itself had they
/// been in the job's process group, and GraphAligner has nothing to flush. The walker exits
/// immediately afterwards, so it does not wait for them.
#[cfg(unix)]
fn stop_children() {
    for pid in lock_aligners().drain(..) {
        // `kill(-pgid)` with pgid 0 or 1 would mean "my own group" and "every process I may
        // signal". Neither can be a child of ours, so neither is worth the risk.
        if pid > 1 {
            // SAFETY: `kill` has no preconditions beyond a valid signal number; an exited
            // group is ESRCH, which is ignored.
            unsafe { libc::kill(-(pid as i32), libc::SIGTERM) };
        }
    }
}

#[cfg(not(unix))]
fn stop_children() {}

const PREFIX: &str = "logan-walker.";
/// Marks a scratch directory as claimed by a live run (WF-01).
const LOCK: &str = ".lock";
/// Written beside it under `LOGAN_KEEP_SCRATCH`, so a kept directory is never swept.
const KEEP: &str = ".keep";

pub struct Scratch {
    path: PathBuf,
    keep: bool,
    /// The run's `.lock`, held open and `flock`ed for the life of the run (WF-01).
    ///
    /// Liveness used to be `kill(pid, 0)`, which sees only the caller's PID namespace on the
    /// caller's host, so a live run in another container or on another node sharing the base
    /// looked dead and had its in-flight directory **deleted**. The victim did not notice:
    /// the directory was silently recreated and the vanished GAF became the biological
    /// answer `no_path`, with exit 0. A held lock is a fact about the filesystem rather than
    /// about one process table.
    _lock: Option<std::fs::File>,
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
        let keep = keep_scratch();
        let pid = std::process::id();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos())
            .unwrap_or(0);
        for attempt in 0..16u32 {
            let path = base.join(format!("{PREFIX}{pid}.{nanos:08x}{attempt:x}"));
            match make_private_dir(&path) {
                Ok(()) => {
                    // Claim it before anything is written into it.
                    let lock = claim(&path, keep);
                    return Ok(Scratch { path, keep, _lock: lock });
                }
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

    /// On a signal that would end the run, stop the run's GraphAligner children and remove
    /// the directory, then end the process as the signal would have, so the caller still sees
    /// "terminated by SIGTERM". The list is in the body, with what sends each one.
    ///
    /// **Signals the process inherited as ignored are left ignored** (WF-02). Registering a
    /// handler for one is a 0.1.1 regression that breaks `nohup`: bash sends SIGHUP to every
    /// job when an SSH session drops, `nohup` exists to set that signal to ignored, and
    /// catching it anyway killed the run with exit 129. Since outputs are written at the end,
    /// the whole run was lost. 0.1.0 had no signal code at all, so this never affected it.
    ///
    /// It matters a second time for the children: a *caught* signal is reset to the default
    /// across `exec`, so GraphAligner stopped inheriting the ignore and died of the same
    /// hangup. Leaving the disposition alone passes the ignore through.
    pub fn remove_on_signal(&self) -> Result<()> {
        use signal_hook::consts::{
            SIGHUP, SIGINT, SIGPIPE, SIGQUIT, SIGTERM, SIGUSR1, SIGUSR2, SIGXCPU, SIGXFSZ,
        };
        // WF-17. Only INT, TERM and HUP were watched, so every other way a run ends killed
        // the process with the directory still in `/dev/shm`, holding RAM until some later
        // run over the same base swept it. Each of these was measured leaving one behind:
        //
        // - SIGPIPE (141), `logan-walker ... | head -1`. `main::default_sigpipe` undoes the
        //   `SIG_IGN` Rust installs at startup, which is what lets this list claim it;
        //   `crate::outln!` then keeps the writing thread from panicking on the EPIPE that
        //   now comes back before the signal is acted on.
        // - SIGQUIT (131), Ctrl-\.
        // - SIGUSR1 (138) and SIGUSR2 (140), which is how SGE and LSF announce a soft
        //   wall-time limit, and SIGXCPU (152), which is a CPU ulimit or a cgroup.
        // - SIGXFSZ, a file-size ulimit, included for the same reason though it was not
        //   measured: the walker writes files.
        //
        // SIGALRM, SIGVTALRM and SIGPROF are deliberately **not** here, though the review
        // listed them. Nothing sends them to announce a limit; what does send them is a timer
        // or a sampling profiler belonging to the process itself, and catching those would
        // end a profiled run rather than clean up after one.
        let wanted: Vec<i32> =
            [SIGINT, SIGTERM, SIGHUP, SIGPIPE, SIGQUIT, SIGUSR1, SIGUSR2, SIGXCPU, SIGXFSZ]
                .into_iter()
                .filter(|&s| !is_ignored(s))
                .collect();
        if wanted.is_empty() {
            // Every signal we would clean up on is ignored, so there is nothing to hook.
            // `Drop` still removes the directory on an ordinary exit.
            return Ok(());
        }
        let mut signals = signal_hook::iterator::Signals::new(&wanted)?;
        let (path, keep) = (self.path.clone(), self.keep);
        std::thread::spawn(move || {
            if let Some(sig) = signals.forever().next() {
                // Before the files, because the aligners are reading and writing them: a
                // child left running after the directory goes is working on nothing.
                stop_children();
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

/// Whether this process inherited `SIG_IGN` for `sig`.
///
/// WF-02. `sigaction` with a null action reads the current disposition without changing it.
fn is_ignored(sig: libc::c_int) -> bool {
    // SAFETY: `sigaction` with a null `act` only reads. `old` is a valid, owned struct, and
    // the call is made before any handler is installed.
    unsafe {
        let mut old: libc::sigaction = std::mem::zeroed();
        if libc::sigaction(sig, std::ptr::null(), &mut old) != 0 {
            return false; // cannot tell; register as before
        }
        old.sa_sigaction == libc::SIG_IGN
    }
}

/// This kernel's boot id, or the hostname if the kernel does not publish one.
///
/// WF-01. A `flock` is only meaningful on the machine that took it: NFS mounted `nolock`, or
/// Lustre mounted `localflock`, does not carry it between nodes. Recording the machine makes
/// the sweep fail SAFE there, leaking directories rather than deleting a live run's.
fn machine_id() -> String {
    if let Ok(s) = std::fs::read_to_string("/proc/sys/kernel/random/boot_id") {
        let s = s.trim().to_string();
        if !s.is_empty() {
            return s;
        }
    }
    std::process::Command::new("hostname")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

/// Take this run's `.lock`, and write `.keep` when the directory is meant to survive.
///
/// A failure here is not fatal: the run proceeds unlocked and the only consequence is that
/// its directory is never swept later. Leaking a directory is the safe side of this trade.
#[cfg(unix)]
fn claim(dir: &Path, keep: bool) -> Option<std::fs::File> {
    use std::io::Write;
    use std::os::unix::io::AsRawFd;
    if keep {
        // `LOGAN_KEEP_SCRATCH` must outlive this process, and a released lock is not enough.
        let _ = std::fs::write(dir.join(KEEP), b"");
    }
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(dir.join(LOCK))
        .ok()?;
    // SAFETY: a valid open fd owned by `file`; LOCK_NB cannot block.
    if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) } != 0 {
        return None;
    }
    file.write_all(machine_id().as_bytes()).ok()?;
    file.flush().ok()?;
    Some(file)
}

#[cfg(not(unix))]
fn claim(_dir: &Path, _keep: bool) -> Option<std::fs::File> {
    None
}

/// Whether `dir`'s lock says it belongs to this machine and nobody is holding it.
///
/// A missing lock means the directory predates WF-01, and is left alone. A lock naming
/// another machine is left alone: we cannot reason about it from here.
#[cfg(unix)]
fn lock_is_free(dir: &Path, me: &str) -> bool {
    use std::os::unix::io::AsRawFd;
    let path = dir.join(LOCK);
    let Ok(owner) = std::fs::read_to_string(&path) else { return false };
    if owner.trim() != me {
        return false;
    }
    let Ok(file) = std::fs::OpenOptions::new().read(true).write(true).open(&path) else {
        return false;
    };
    // SAFETY: a valid open fd owned by `file`; LOCK_NB cannot block.
    if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) } != 0 {
        return false; // someone is running in there
    }
    // SAFETY: the same fd, releasing the probe lock.
    unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_UN) };
    true
}

/// `remove_dir_all`, retried: a worker can create a file in the directory while it is being
/// emptied, which makes an attempt fail with "directory not empty".
///
/// Two attempts back to back were not always enough on the signal path, where several
/// compute workers may be mid-write; 20 tries over 100 ms covers the gap between a worker's
/// last write and its death without making an ordinary exit wait for anything (the first
/// attempt succeeds, and the loop stops as soon as the directory is gone).
fn remove(path: &Path) {
    for attempt in 0..20 {
        if std::fs::remove_dir_all(path).is_ok() || !path.exists() {
            return;
        }
        if attempt < 19 {
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
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

/// Remove run directories in `base` that are ours, on this machine, and not in use.
///
/// WF-01. This used to decide with `kill(pid, 0)`, which answers only about the caller's PID
/// namespace on the caller's host: a live run in another container, or on another node
/// sharing the base over Lustre or NFS, looked dead and had its directory deleted mid-flight.
/// The victim recreated the directory silently and reported the affected accessions as
/// `no_path`, the biological answer, with exit 0. Reproduced twice; one reproduction lost 4
/// of 8 accessions.
///
/// Four conditions now, all required:
///
/// 1. the name is exactly what `create` builds, `logan-walker.<digits>.<hex>`;
/// 2. it is a directory we own;
/// 3. it holds a `.lock` naming **this** machine, so another node's directory is never
///    touched even where its lock is invisible to us;
/// 4. that lock is free, which `flock(LOCK_NB)` establishes, and there is no `.keep`.
///
/// This run protects itself for free: it holds its own lock, so condition 4 refuses it. A
/// directory left by 0.1.x has no lock and is never swept, so the release notes should ask
/// for a one-time manual cleanup. `process_exists` is gone: a PID is not evidence.
#[cfg(unix)]
fn sweep_stale(base: &Path) -> usize {
    use std::os::unix::fs::MetadataExt;
    let Ok(entries) = std::fs::read_dir(base) else { return 0 };
    // SAFETY: getuid cannot fail and has no preconditions.
    let uid = unsafe { libc::getuid() };
    let me = machine_id();
    let mut removed = 0;
    for entry in entries.flatten() {
        let name = entry.file_name();
        let Some(rest) = name.to_str().and_then(|n| n.strip_prefix(PREFIX)) else { continue };
        // `<digits>.<hex>`, which is what `create` builds and nothing else.
        let Some((pid, suffix)) = rest.split_once('.') else { continue };
        if pid.is_empty()
            || !pid.bytes().all(|b| b.is_ascii_digit())
            || suffix.is_empty()
            || !suffix.bytes().all(|b| b.is_ascii_hexdigit())
        {
            continue;
        }
        let Ok(meta) = entry.metadata() else { continue };
        if !meta.is_dir() || meta.uid() != uid {
            continue;
        }
        let dir = entry.path();
        if dir.join(KEEP).exists() || !lock_is_free(&dir, &me) {
            continue;
        }
        remove(&dir);
        removed += 1;
    }
    removed
}

#[cfg(not(unix))]
fn sweep_stale(_base: &Path) -> usize {
    0
}


#[cfg(all(test, unix))]
mod signal_tests {
    use super::is_ignored;
    use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};

    /// WF-02. `nohup` works by setting SIGHUP to ignored before exec. The walker registered a
    /// handler for it regardless, so bash's hangup on an SSH drop killed the run with exit
    /// 129 and, because outputs are written at the end, lost all of it. The fix is to read
    /// the inherited disposition first, which is what this asserts.
    #[test]
    fn an_inherited_ignore_is_detected() {
        assert!(!is_ignored(SIGTERM), "SIGTERM is not ignored by default");

        // Ignore SIGHUP the way `nohup` does, then check we can see it.
        // SAFETY: a test-local disposition change on a signal this test then restores.
        unsafe {
            let mut act: libc::sigaction = std::mem::zeroed();
            act.sa_sigaction = libc::SIG_IGN;
            let mut old: libc::sigaction = std::mem::zeroed();
            assert_eq!(libc::sigaction(SIGHUP, &act, &mut old), 0);
            assert!(is_ignored(SIGHUP), "an inherited SIG_IGN must be visible");
            assert_eq!(libc::sigaction(SIGHUP, &old, std::ptr::null_mut()), 0);
        }
        assert!(!is_ignored(SIGHUP), "restored to the default");
        assert!(!is_ignored(SIGINT));
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    /// A base directory for one test, removed when the returned guard drops (TW-10).
    ///
    /// The manual `remove_dir_all` these tests used to end with ran only when the test
    /// reached it, so an assertion failure leaked the directory, and the sweep test leaked
    /// one on every run: it deliberately leaves behind the directories the sweep must not
    /// touch, and so could never have ended by removing its base.
    fn base(name: &str) -> TestDir {
        TestDir::new(name)
    }

    #[test]
    fn a_private_directory_removed_on_drop() {
        let d = base("drop");
        let b = d.path();
        let s = Scratch::create(Some(&b)).unwrap();
        let p = s.path().to_path_buf();
        assert!(p.starts_with(&b) && p.is_dir());
        assert_eq!(std::fs::metadata(&p).unwrap().permissions().mode() & 0o777, 0o700);
        std::fs::write(p.join("x.sub.gfa"), b"H\n").unwrap();
        drop(s);
        assert!(!p.exists(), "the run's directory must be removed with its files");
    }

    #[test]
    fn runs_get_separate_directories() {
        let d = base("two");
        let b = d.path();
        let (s1, s2) = (Scratch::create(Some(&b)).unwrap(), Scratch::create(Some(&b)).unwrap());
        assert_ne!(s1.path(), s2.path());
        drop((s1, s2));
    }

    /// WF-01. Liveness is a held lock, not a PID.
    ///
    /// The old rule was `kill(pid, 0)`, which answers only about the caller's PID namespace
    /// on the caller's host, so a live run in another container or on another node sharing
    /// the base looked dead and had its in-flight directory deleted. Each case below turns on
    /// something that survives a namespace boundary; none turns on a PID.
    #[test]
    fn sweeping_is_decided_by_the_lock_not_the_pid() {
        use std::io::Write;
        use std::os::unix::io::AsRawFd;
        let d = base("sweep");
        let b = d.path();
        std::fs::create_dir_all(&b).unwrap();
        let me = super::machine_id();

        // A run that ended: lock free. Its PID is 1, which is always alive, to prove the PID
        // is no longer consulted.
        let ended = b.join(format!("{PREFIX}1.0a"));
        std::fs::create_dir(&ended).unwrap();
        std::fs::write(ended.join(LOCK), me.as_bytes()).unwrap();
        std::fs::write(ended.join("left.sub.gfa"), b"H\n").unwrap();

        // A run in flight: lock HELD. Its PID cannot exist.
        let live = b.join(format!("{PREFIX}99999999.0b"));
        std::fs::create_dir(&live).unwrap();
        let mut held = std::fs::OpenOptions::new()
            .write(true).create(true).truncate(true).open(live.join(LOCK)).unwrap();
        held.write_all(me.as_bytes()).unwrap();
        held.flush().unwrap();
        // SAFETY: a valid open fd owned by `held`, held for the rest of the test.
        assert_eq!(unsafe { libc::flock(held.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) }, 0);

        // A run on another node sharing this base over a cluster filesystem.
        let elsewhere = b.join(format!("{PREFIX}4242.0c"));
        std::fs::create_dir(&elsewhere).unwrap();
        std::fs::write(elsewhere.join(LOCK), b"another-kernel-boot-id").unwrap();

        // Left by 0.1.x: no lock at all.
        let legacy = b.join(format!("{PREFIX}4243.0d"));
        std::fs::create_dir(&legacy).unwrap();

        // LOGAN_KEEP_SCRATCH asked for this one to survive.
        let kept = b.join(format!("{PREFIX}4244.0e"));
        std::fs::create_dir(&kept).unwrap();
        std::fs::write(kept.join(LOCK), me.as_bytes()).unwrap();
        std::fs::write(kept.join(KEEP), b"").unwrap();

        // Names that are not ours, including one that only looks close.
        let other = b.join("not-ours");
        std::fs::create_dir(&other).unwrap();
        let nearly = b.join(format!("{PREFIX}notanumber.0f"));
        std::fs::create_dir(&nearly).unwrap();

        let s = Scratch::create(Some(&b)).unwrap();

        assert!(!ended.exists(), "a free lock on this machine is swept, whatever the PID");
        assert!(live.exists(), "a HELD lock is a live run, even at an impossible PID");
        assert!(elsewhere.exists(), "another machine's directory is never touched");
        assert!(legacy.exists(), "a 0.1.x directory has no lock and is left for the user");
        assert!(kept.exists(), "LOGAN_KEEP_SCRATCH survives the next run's sweep");
        assert!(other.exists() && nearly.exists(), "only our own name pattern is considered");
        assert!(s.path().join(LOCK).exists(), "and this run claims its own directory");

        drop(held);
        drop(s);
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn an_unwritable_base_names_itself() {
        let d = base("ro");
        let b = d.path();
        std::fs::create_dir_all(&b).unwrap();
        std::fs::set_permissions(&b, std::fs::Permissions::from_mode(0o500)).unwrap();
        let err = Scratch::create(Some(&b)).err().map(|e| e.to_string());
        std::fs::set_permissions(&b, std::fs::Permissions::from_mode(0o700)).unwrap();
        if unsafe { libc::geteuid() } == 0 {
            return; // root writes anywhere
        }
        let err = err.expect("an unwritable base must fail");
        assert!(err.contains(&b.display().to_string()) && err.contains("--scratch-dir"), "{err}");
    }
}
