# Validation

MD Workstation keeps validation close to the codebase: unit tests cover core
math and file formats, curated TOML examples cover user-facing workflows, and
release automation builds the CLI and optional desktop UI on supported
platforms.

## Local Checks

Run these before changing shared simulation behavior:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --locked
cargo run --locked -p md-cli -- validate-suite
```

For desktop UI work, also run:

```bash
cargo check --locked -p md-ui --features desktop
```

For release-oriented work, build both shipped binaries:

```bash
cargo build --locked --release -p md-cli --bin md
cargo build --locked --release -p md-ui --features desktop --bin md-ui
```

## Curated Validation Suite

The standard suite is intentionally compact:

```bash
cargo run -p md-cli -- validate-suite
```

It validates deterministic examples for:

- short NVE dynamics
- cold generated lattice setup
- periodic boundaries
- neighbor-list force setup
- topology-aware force setup
- optional parallel execution
- XYZ and PDB coordinate inputs
- bonded terms, Coulomb, exclusions, and force-field defaults
- Berendsen thermostat configuration
- minimization configuration
- binary checkpoint configuration

The long suite adds slower stress-oriented configs:

```bash
cargo run -p md-cli -- validate-suite --long
```

`--long` includes `examples/validation/long-nve-drift.toml` and
`examples/validation/larger-neighbor-compare.toml`.

## Manual Stress Checks

`examples/validation/long-nve-drift-100k.toml` is intentionally outside the
automated suite. Use it when you want a slower manual energy-drift check:

```bash
cargo run -p md-cli -- run examples/validation/long-nve-drift-100k.toml
```

For timestep sensitivity work, keep the baseline `dt = 0.0001` and vary one
parameter at a time. Useful manual axes are:

- `simulation.steps`: 10000, 100000, or larger
- `[neighbor].enabled`: `true` or `false`
- `[execution].parallel`: `true` or `false`
- `[box].periodic`: `true` or `false`

## Force-Path Benchmark Validation

Use `bench-neighbor` to compare serial, parallel, naive, and neighbor-list force
paths for one config:

```bash
cargo run -p md-cli -- bench-neighbor examples/validation/larger-neighbor-compare.toml --repeats 10
```

The benchmark reports timing plus energy and max-force deltas against the serial
naive baseline. Treat the timings as local-machine guidance, not as portable
performance claims.

## CI And Release Validation

GitHub Actions runs the CI workflow on pushes, pull requests, and manual
dispatches. It checks formatting, Clippy, workspace tests, the validation suite,
and the desktop UI feature build.

The release workflow is triggered by `v*` tags or manual dispatch for an
existing tag. It validates the workspace, builds `md` and `md-ui`, packages
Linux, Windows, and macOS artifacts with project docs and examples when present,
and uploads SHA-256 checksums.

## Interpreting Results

- The standard suite is a compact smoke test for curated configs and force
  setup; it is not exhaustive scientific validation.
- The long suite and manual stress configs are intended for slower regression
  checks that do not need to run on every edit.
- Force-path benchmarks should be read as consistency checks first and timing
  measurements second.
- Passing validation means the documented examples still load and exercise the
  expected paths; it does not turn the project into a production MD engine.
