//! Writing to stdout without dying of it.

/// `println!` that does not panic when stdout is a closed pipe.
///
/// WF-17, WP-26. Rust sets SIGPIPE to `SIG_IGN` at startup, so a write to a closed pipe
/// returns `EPIPE` rather than killing the process, and `println!` turns that into a panic:
/// `failed printing to stdout: Broken pipe`. `main` therefore restored the default action, so
/// `logan-walker ... | head -1` ended the way any Unix tool does. That worked, but it ended the
/// process **immediately**, in the middle of a write, with nothing run on the way out: the
/// run's scratch directory stayed in `/dev/shm` holding RAM until a later run swept it, and
/// the GraphAligner children kept going and kept writing into it.
///
/// SIGPIPE is now handled like SIGINT and SIGTERM (see [`crate::scratch`]), which means the
/// write returns `EPIPE` first and the cleanup thread ends the process a moment later. In that
/// moment the writing thread must not panic, so every progress line goes through here and
/// simply drops the error. The next thing to happen is the process exiting 141.
///
/// Only progress output uses it. Warnings and errors go to stderr, which is a different
/// stream with a different failure mode, and the summary and FASTA are files.
#[macro_export]
macro_rules! outln {
    ($($arg:tt)*) => {{
        use std::io::Write;
        let _ = writeln!(std::io::stdout().lock(), $($arg)*);
    }};
}
