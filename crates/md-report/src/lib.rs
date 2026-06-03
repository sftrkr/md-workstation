//! Markdown report generation for completed MD run directories.

use std::fmt::Write;

use md_analysis::AnalysisSummary;
use md_io::{RunManifest, RunSummary};

#[derive(Debug, Clone, Copy)]
pub struct MarkdownReportInput<'a> {
    pub run_dir: &'a str,
    pub summary: &'a RunSummary,
    pub manifest: Option<&'a RunManifest>,
    pub analysis: Option<&'a AnalysisSummary>,
}

pub fn render_markdown_report(input: &MarkdownReportInput<'_>) -> String {
    let summary = input.summary;
    let mut report = String::new();

    writeln!(report, "# MD Run Report").unwrap();
    writeln!(report).unwrap();
    writeln!(report, "## Run").unwrap();
    writeln!(report).unwrap();
    writeln!(report, "- Run name: {}", summary.run_name).unwrap();
    writeln!(report, "- Run directory: {}", input.run_dir).unwrap();
    if let Some(project_name) = &summary.project_name {
        writeln!(report, "- Project: {project_name}").unwrap();
    }
    if let Some(description) = &summary.project_description {
        writeln!(report, "- Description: {description}").unwrap();
    }
    if let Some(manifest) = input.manifest {
        writeln!(report, "- Status: {}", manifest.status).unwrap();
        writeln!(report, "- Manifest schema: {}", manifest.schema_version).unwrap();
    }

    writeln!(report).unwrap();
    writeln!(report, "## System").unwrap();
    writeln!(report).unwrap();
    writeln!(report, "- Particles: {}", summary.particle_count).unwrap();
    writeln!(report, "- Steps: {}", summary.steps).unwrap();
    writeln!(report, "- Time step: {:.10}", summary.dt).unwrap();
    writeln!(report, "- Output interval: {}", summary.output_interval).unwrap();
    writeln!(report, "- Seed: {}", summary.seed).unwrap();
    writeln!(report, "- Workflow: {}", summary.workflow).unwrap();
    writeln!(report, "- Units: {}", summary.units).unwrap();
    writeln!(report, "- Boundary: {}", summary.boundary).unwrap();
    writeln!(report, "- Ensemble: {}", summary.ensemble).unwrap();
    writeln!(report, "- Thermostat: {}", summary.thermostat).unwrap();
    if let Some(target_temperature) = summary.thermostat_target_temperature {
        writeln!(
            report,
            "- Thermostat target temperature: {:.10}",
            target_temperature
        )
        .unwrap();
    }
    if let Some(tau) = summary.thermostat_tau {
        writeln!(report, "- Thermostat tau: {:.10}", tau).unwrap();
    }
    if let Some(step_size) = summary.minimization_step_size {
        writeln!(report, "- Minimization step size: {:.10}", step_size).unwrap();
    }
    if let Some(force_tolerance) = summary.minimization_force_tolerance {
        writeln!(
            report,
            "- Minimization force tolerance: {:.10}",
            force_tolerance
        )
        .unwrap();
    }
    if let Some(final_max_force) = summary.minimization_final_max_force {
        writeln!(
            report,
            "- Minimization final max force: {:.10}",
            final_max_force
        )
        .unwrap();
    }

    writeln!(report).unwrap();
    writeln!(report, "## Inputs").unwrap();
    writeln!(report).unwrap();
    match &summary.input_file {
        Some(input_file) => {
            writeln!(report, "- Coordinate input: {input_file}").unwrap();
            if let Some(format) = &summary.input_format {
                writeln!(report, "- Coordinate format: {format}").unwrap();
            }
            if let Some(frame) = summary.input_frame {
                writeln!(report, "- Coordinate frame: {frame}").unwrap();
            }
            if let Some(atom_metadata_file) = &summary.atom_metadata_file {
                writeln!(report, "- Atom metadata: {atom_metadata_file}").unwrap();
            }
        }
        None => {
            writeln!(report, "- Coordinate input: generated lattice").unwrap();
        }
    }
    match &summary.topology_file {
        Some(topology_file) => writeln!(report, "- Topology: {topology_file}").unwrap(),
        None => writeln!(report, "- Topology: none").unwrap(),
    }

    writeln!(report).unwrap();
    writeln!(report, "## Force Field").unwrap();
    writeln!(report).unwrap();
    writeln!(report, "- Non-bonded force: {}", summary.force).unwrap();
    writeln!(
        report,
        "- Shifted LJ potential: {}",
        summary.shifted_potential
    )
    .unwrap();
    writeln!(report, "- Harmonic bonds: {}", summary.bond_count).unwrap();
    writeln!(report, "- Harmonic angles: {}", summary.angle_count).unwrap();
    writeln!(report, "- Periodic dihedrals: {}", summary.dihedral_count).unwrap();
    writeln!(
        report,
        "- Non-bonded exclusions: {}",
        summary.excluded_pair_count
    )
    .unwrap();
    writeln!(report, "- Coulomb enabled: {}", summary.coulomb).unwrap();
    if let Some(cutoff) = summary.coulomb_cutoff {
        writeln!(report, "- Coulomb cutoff: {:.10}", cutoff).unwrap();
    }
    writeln!(report, "- Neighbor list: {}", summary.neighbor_list).unwrap();
    if let Some(skin) = summary.neighbor_skin {
        writeln!(report, "- Neighbor skin: {:.10}", skin).unwrap();
    }
    if let Some(interval) = summary.neighbor_rebuild_interval {
        writeln!(report, "- Neighbor rebuild interval: {}", interval).unwrap();
    }
    writeln!(report, "- Parallel force: {}", summary.parallel).unwrap();
    if let Some(threads) = summary.rayon_threads {
        writeln!(report, "- Rayon threads: {}", threads).unwrap();
    }

    writeln!(report).unwrap();
    writeln!(report, "## Final Energy").unwrap();
    writeln!(report).unwrap();
    writeln!(report, "- Step: {}", summary.final_energy.step).unwrap();
    writeln!(report, "- Time: {:.10}", summary.final_energy.time).unwrap();
    writeln!(report, "- Kinetic: {:.10}", summary.final_energy.kinetic).unwrap();
    writeln!(
        report,
        "- Potential: {:.10}",
        summary.final_energy.potential
    )
    .unwrap();
    writeln!(report, "- Total: {:.10}", summary.final_energy.total).unwrap();
    writeln!(
        report,
        "- Temperature: {:.10}",
        summary.final_energy.temperature
    )
    .unwrap();

    if let Some(analysis) = input.analysis {
        writeln!(report).unwrap();
        writeln!(report, "## Analysis").unwrap();
        writeln!(report).unwrap();
        writeln!(report, "- Frames: {}", analysis.frame_count).unwrap();
        writeln!(report, "- Atoms: {}", analysis.atom_count).unwrap();
        writeln!(report, "- RMSD alignment: {}", analysis.rmsd_alignment).unwrap();
        writeln!(report, "- Periodic unwrap: {}", analysis.periodic_unwrap).unwrap();
        if let Some([x, y, z]) = analysis.periodic_box {
            writeln!(report, "- Periodic box: {:.10} x {:.10} x {:.10}", x, y, z).unwrap();
        }
        writeln!(
            report,
            "- Energy drift: {:.10}",
            analysis.energy.total_drift
        )
        .unwrap();
        if let Some(drift_per_time) = analysis.energy.total_drift_per_time {
            writeln!(report, "- Energy drift per time: {:.10}", drift_per_time).unwrap();
        }
        writeln!(
            report,
            "- Mean total energy: {:.10}",
            analysis.energy.mean_total
        )
        .unwrap();
        writeln!(report, "- Final RMSD: {:.10}", analysis.rmsd.final_value).unwrap();
        if let Some(pair_distance) = &analysis.pair_distance {
            writeln!(
                report,
                "- Pair {}-{} final distance: {:.10}",
                pair_distance.atom_i, pair_distance.atom_j, pair_distance.final_value
            )
            .unwrap();
        }
    }

    writeln!(report).unwrap();
    writeln!(report, "## Outputs").unwrap();
    writeln!(report).unwrap();
    for output in &summary.outputs {
        writeln!(report, "- {output}").unwrap();
    }
    if let Some(manifest) = input.manifest {
        for output in &manifest.outputs {
            if !summary.outputs.contains(output) {
                writeln!(report, "- {output}").unwrap();
            }
        }
    }

    writeln!(report).unwrap();
    writeln!(report, "## Limitations").unwrap();
    writeln!(report).unwrap();
    writeln!(
        report,
        "This run was produced by an educational Rust-native MD engine. Results should not be used for clinical, pharmaceutical, or other high-stakes decisions without independent validation."
    )
    .unwrap();

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_core::EnergySample;

    #[test]
    fn markdown_report_includes_core_run_fields() {
        let summary = RunSummary {
            run_name: "phase9-test".to_string(),
            project_name: Some("demo".to_string()),
            project_description: Some("report smoke test".to_string()),
            particle_count: 2,
            steps: 10,
            dt: 0.001,
            output_interval: 5,
            seed: 7,
            workflow: "dynamics".to_string(),
            ensemble: "NVT".to_string(),
            thermostat: "berendsen".to_string(),
            thermostat_target_temperature: Some(0.2),
            thermostat_tau: Some(0.1),
            minimization_step_size: None,
            minimization_force_tolerance: None,
            minimization_final_max_force: None,
            input_file: Some("input.xyz".to_string()),
            input_format: Some("xyz".to_string()),
            input_frame: Some(0),
            atom_metadata_file: None,
            topology_file: Some("topology.toml".to_string()),
            bond_count: 1,
            angle_count: 1,
            dihedral_count: 1,
            excluded_pair_count: 2,
            coulomb: true,
            coulomb_cutoff: Some(3.0),
            force: "lennard-jones".to_string(),
            boundary: "open".to_string(),
            shifted_potential: false,
            neighbor_list: false,
            neighbor_skin: None,
            neighbor_rebuild_interval: None,
            parallel: false,
            rayon_threads: None,
            units: "reduced-lennard-jones".to_string(),
            final_energy: EnergySample::new(10, 0.01, 1.0, -0.5, 0.1),
            outputs: vec!["summary.json".to_string()],
        };
        let report = render_markdown_report(&MarkdownReportInput {
            run_dir: "runs/phase9-test-001",
            summary: &summary,
            manifest: None,
            analysis: None,
        });

        assert!(report.contains("# MD Run Report"));
        assert!(report.contains("- Project: demo"));
        assert!(report.contains("- Workflow: dynamics"));
        assert!(report.contains("- Ensemble: NVT"));
        assert!(report.contains("- Thermostat: berendsen"));
        assert!(report.contains("- Thermostat target temperature: 0.2000000000"));
        assert!(report.contains("- Harmonic bonds: 1"));
        assert!(report.contains("- Harmonic angles: 1"));
        assert!(report.contains("- Periodic dihedrals: 1"));
        assert!(report.contains("- Non-bonded exclusions: 2"));
        assert!(report.contains("- Coulomb enabled: true"));
    }
}
