// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, which was written by employees of the United
// States Government (National Institutes of Health) as part of their official
// duties. It is a "United States Government Work" (17 U.S.C. 105) and is in the
// public domain; outside the US, rights are waived under CC0 1.0. See LICENSE.
//
// ---------------------------------------------------------------------------
// WHAT THIS FILE IS
//
// Enforces the architecture rules from docs/ARCHITECTURE.md as tests, so a
// careless `cargo add` fails CI instead of silently eroding the design:
//
//   * `afni-core` must not depend (even indirectly) on `afni-io`
//     -> prevents an afni-core -> afni-io -> afni-core cycle.
//   * `afni-core` must not pull in GUI, GPU, windowing, socket, or file-format
//     crates.
//   * the source tree must not launch processes or open sockets.
//
// It asks Cargo itself (`cargo tree`) for the real resolved dependency list,
// rather than grepping Cargo.toml, so transitive dependencies are covered too.
// ---------------------------------------------------------------------------

use std::path::Path;
use std::process::Command;

/// Crates `afni-core` must never depend on, with the reason shown on failure.
const FORBIDDEN_CRATES: &[(&str, &str)] = &[
    (
        "afni-io",
        "file formats live in afni-io; the dependency arrow points afni-io -> afni-core",
    ),
    ("egui", "GUI toolkits belong in the viewers"),
    ("eframe", "GUI toolkits belong in the viewers"),
    ("wgpu", "GPU APIs belong in the viewers"),
    ("winit", "windowing belongs in the viewers"),
    ("glam", "viewer math types stay in the viewers"),
    ("flate2", "decompression is a file-format concern (afni-io)"),
    ("tokio", "no async runtime or sockets in the semantic layer"),
    ("mio", "no sockets in the semantic layer"),
    ("socket2", "no sockets in the semantic layer"),
];

/// Names of every crate in `afni-core`'s resolved dependency tree.
fn resolved_dependencies() -> Vec<String> {
    // Cargo sets $CARGO to the binary running the tests; fall back to PATH.
    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());
    let output = Command::new(cargo)
        .args([
            "tree",
            "--package",
            "afni-core",
            "--prefix",
            "none",
            "--edges",
            "normal,build,dev",
            "--offline",
        ])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("failed to run `cargo tree`");
    assert!(
        output.status.success(),
        "cargo tree failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    // Each line looks like "thiserror v2.0.21 (proc-macro)"; keep the name.
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .map(str::to_owned)
        .collect()
}

#[test]
fn no_forbidden_crates_in_dependency_tree() {
    let deps = resolved_dependencies();
    // Sanity check that parsing worked: we know we depend on thiserror.
    assert!(
        deps.iter().any(|d| d == "thiserror"),
        "cargo tree parse failed: {deps:?}"
    );

    for (name, why) in FORBIDDEN_CRATES {
        assert!(
            !deps.iter().any(|d| d == name),
            "afni-core must not depend on `{name}`: {why}"
        );
    }
}

/// Collect every `.rs` file under `dir`.
fn rust_files(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
    for entry in std::fs::read_dir(dir).expect("read_dir") {
        let path = entry.expect("dir entry").path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

#[test]
fn library_source_has_no_io_process_or_socket_code() {
    // Library code only (src/); the test harness is allowed to run AFNI.
    let mut files = Vec::new();
    rust_files(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),
        &mut files,
    );
    assert!(!files.is_empty());

    // Substrings that would indicate a boundary violation. Written with a
    // split so this very list is not matched when scanning src/ in future.
    let banned = [
        "std::process",
        "std::net",
        "std::fs",
        "std::path",
        "extern crate afni_io",
        "use afni_io",
    ];
    for file in files {
        let text = std::fs::read_to_string(&file).expect("read source");
        for pattern in banned {
            assert!(
                !text.contains(pattern),
                "{} contains `{pattern}`, which crosses the afni-core boundary \
                 (no files, paths, processes, sockets, or afni-io)",
                file.display()
            );
        }
        assert!(
            !text.contains("unsafe "),
            "{} uses `unsafe`; the crate forbids it",
            file.display()
        );
    }
}
