# Limitations

MD Workstation is a transparent, small-scale Rust MD workstation. It is useful
for education, prototyping, force-field experiments, validation exercises, and
reproducible reports. It is not a production molecular simulation package.

## Scientific Scope

- Reduced Lennard-Jones units only; real unit conversion is not implemented.
- Temperature uses a simple 3N degree-of-freedom estimate and does not yet
  subtract center-of-mass or constraint degrees of freedom.
- Only orthorhombic boxes are supported.
- Periodic boundaries use minimum-image pair distances.
- Shifted LJ and Coulomb potentials are available; force shifting is not
  implemented.
- Thermostat support is limited to deterministic Berendsen velocity coupling.
  It is not a rigorous canonical ensemble sampler.
- There is no barostat, pressure control, virial reporting, NPT ensemble, Ewald,
  or Particle Mesh Ewald.

## Force Fields And Topology

- Topology support is project-local TOML, not a broad force-field format.
- There is no residue template assignment or automatic atom typing.
- There are no production biomolecular force-field compatibility claims.
- Coulomb is cutoff-based only.
- Bonded terms are intentionally basic: harmonic bonds, harmonic angles, and a
  periodic dihedral term.
- Constraints such as SHAKE or RATTLE are not implemented.

## Molecular File Support

- XYZ support is intentionally simple.
- PDB support is a small coordinate subset, not full PDB compliance.
- PDB unit cells, connectivity records, alternate locations, insertion codes,
  occupancies, and B-factors are not interpreted.
- mmCIF, GRO, DCD, XTC, TRR, and other trajectory formats are not supported.

## Performance Scope

- Parallelism is CPU-only and host-local through Rayon.
- There is no GPU backend.
- There is no SIMD-specific kernel selection.
- There is no domain decomposition.
- Benchmarks are useful for local comparisons, not portable performance claims.

## Analysis Assumptions

- Analysis assumes stable atom ordering and atom count.
- RMSD uses centered rigid alignment against a selected reference frame.
- Periodic unwrap assumes particles move less than half the box length between
  saved frames.
- RDF, RMSF, radius of gyration, and other richer analyses are future work.

## Reproducibility Boundaries

- Manifests record hashes, command line, engine version, git commit when
  available, platform metadata, and Rayon thread count.
- Reproducibility checks verify copied config/input/topology files against the
  manifest. They do not guarantee bit-identical floating-point behavior across
  every CPU, compiler, or OS.
- Binary checkpoints are project-local and not a stable cross-tool interchange
  format.

## Explicit Non-Goals

- Replacing established production MD engines or molecular conversion tools.
- Clinical, pharmaceutical, or drug-design validation.
- Full biomolecular preparation workflows.
- External runtime, container, VM, or proprietary toolchain dependency.
- Production-scale GPU or distributed MD.
