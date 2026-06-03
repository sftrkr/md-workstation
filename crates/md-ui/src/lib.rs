use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use md_core::EnergySample;
use md_io::{
    read_energy_csv, read_manifest_json, read_summary_json, read_xyz_coordinate_frames,
    RunManifest, RunSummary,
};

#[cfg(feature = "desktop")]
pub mod desktop;

#[derive(Debug, Clone)]
pub struct EnergyPoint {
    pub step: usize,
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
}

#[derive(Debug, Clone)]
pub struct RunPreview {
    pub path: PathBuf,
    pub label: String,
    pub summary: Option<RunSummary>,
    pub manifest: Option<RunManifest>,
    pub energy: Vec<EnergyPoint>,
    pub trajectory: Option<TrajectoryPreview>,
    pub report: Option<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct TrajectoryPreview {
    pub frame_count: usize,
    pub particle_count: usize,
    pub first_step: Option<usize>,
    pub last_step: Option<usize>,
    pub first_atoms: Vec<String>,
}

pub fn discover_validation_configs(workspace_root: &Path) -> Result<Vec<PathBuf>> {
    let validation_dir = workspace_root.join("examples").join("validation");
    discover_files_with_extension(&validation_dir, "toml")
}

pub fn discover_run_dirs(workspace_root: &Path) -> Result<Vec<PathBuf>> {
    let runs_dir = workspace_root.join("runs");
    if !runs_dir.exists() {
        return Ok(Vec::new());
    }

    let mut dirs = Vec::new();
    for entry in
        fs::read_dir(&runs_dir).with_context(|| format!("failed to read {}", runs_dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir()
            && (path.join("summary.json").exists() || path.join("run-manifest.json").exists())
        {
            dirs.push(path);
        }
    }
    dirs.sort();
    Ok(dirs)
}

pub fn load_run_preview(run_dir: &Path) -> RunPreview {
    let mut warnings = Vec::new();
    let summary = read_optional_summary(run_dir, &mut warnings);
    let manifest = read_optional_manifest(run_dir, &mut warnings);
    let energy = read_optional_energy(run_dir, &mut warnings);
    let trajectory = read_optional_trajectory(run_dir, &mut warnings);
    let report = read_optional_report(run_dir, &mut warnings);
    let label = summary
        .as_ref()
        .map(|summary| summary.run_name.clone())
        .or_else(|| manifest.as_ref().map(|manifest| manifest.run_name.clone()))
        .unwrap_or_else(|| {
            run_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("run")
                .to_string()
        });

    RunPreview {
        path: run_dir.to_path_buf(),
        label,
        summary,
        manifest,
        energy,
        trajectory,
        report,
        warnings,
    }
}

pub fn load_config_text(config_path: &Path) -> Result<String> {
    fs::read_to_string(config_path)
        .with_context(|| format!("failed to read config {}", config_path.display()))
}

pub fn save_config_text(config_path: &Path, text: &str) -> Result<()> {
    fs::write(config_path, text)
        .with_context(|| format!("failed to write config {}", config_path.display()))
}

pub fn display_path(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn discover_files_with_extension(root: &Path, extension: &str) -> Result<Vec<PathBuf>> {
    if !root.exists() {
        return Ok(Vec::new());
    }

    let mut files = Vec::new();
    for entry in fs::read_dir(root).with_context(|| format!("failed to read {}", root.display()))? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() && path.extension().and_then(|value| value.to_str()) == Some(extension) {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn read_optional_summary(run_dir: &Path, warnings: &mut Vec<String>) -> Option<RunSummary> {
    let path = run_dir.join("summary.json");
    if !path.exists() {
        return None;
    }
    read_summary_json(&path)
        .map_err(|error| warnings.push(format!("summary.json: {error}")))
        .ok()
}

fn read_optional_manifest(run_dir: &Path, warnings: &mut Vec<String>) -> Option<RunManifest> {
    let path = run_dir.join("run-manifest.json");
    if !path.exists() {
        return None;
    }
    read_manifest_json(&path)
        .map_err(|error| warnings.push(format!("run-manifest.json: {error}")))
        .ok()
}

fn read_optional_energy(run_dir: &Path, warnings: &mut Vec<String>) -> Vec<EnergyPoint> {
    let path = run_dir.join("energy.csv");
    if !path.exists() {
        return Vec::new();
    }
    read_energy_csv(&path)
        .map(|samples| samples.into_iter().map(EnergyPoint::from).collect())
        .map_err(|error| warnings.push(format!("energy.csv: {error}")))
        .unwrap_or_default()
}

fn read_optional_trajectory(
    run_dir: &Path,
    warnings: &mut Vec<String>,
) -> Option<TrajectoryPreview> {
    let path = run_dir.join("trajectory.xyz");
    if !path.exists() {
        return None;
    }
    let frames = read_xyz_coordinate_frames(&path)
        .map_err(|error| warnings.push(format!("trajectory.xyz: {error}")))
        .ok()?;
    let first_frame = frames.first();
    let last_frame = frames.last();
    let first_atoms = first_frame
        .map(|frame| {
            frame
                .atoms
                .iter()
                .take(8)
                .map(|atom| format!("{} {:.6} {:.6} {:.6}", atom.element, atom.x, atom.y, atom.z))
                .collect()
        })
        .unwrap_or_default();

    Some(TrajectoryPreview {
        frame_count: frames.len(),
        particle_count: first_frame.map_or(0, |frame| frame.atoms.len()),
        first_step: first_frame.map(|frame| frame.step),
        last_step: last_frame.map(|frame| frame.step),
        first_atoms,
    })
}

fn read_optional_report(run_dir: &Path, warnings: &mut Vec<String>) -> Option<String> {
    let path = run_dir.join("run-report.md");
    if !path.exists() {
        return None;
    }
    fs::read_to_string(&path)
        .map_err(|error| warnings.push(format!("run-report.md: {error}")))
        .ok()
}

impl From<EnergySample> for EnergyPoint {
    fn from(sample: EnergySample) -> Self {
        Self {
            step: sample.step,
            time: sample.time,
            kinetic: sample.kinetic,
            potential: sample.potential,
            total: sample.total,
            temperature: sample.temperature,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn discovers_run_dirs_with_summary_or_manifest() {
        let root = unique_temp_dir("md-ui-runs");
        let runs = root.join("runs");
        fs::create_dir_all(runs.join("a-run")).unwrap();
        fs::create_dir_all(runs.join("b-run")).unwrap();
        fs::create_dir_all(runs.join("scratch")).unwrap();
        fs::write(runs.join("a-run/summary.json"), "{}").unwrap();
        fs::write(runs.join("b-run/run-manifest.json"), "{}").unwrap();

        let found = discover_run_dirs(&root).unwrap();

        assert_eq!(found.len(), 2);
        assert!(found.iter().any(|path| path.ends_with("a-run")));
        assert!(found.iter().any(|path| path.ends_with("b-run")));

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn load_run_preview_reads_energy_points_and_report() {
        let root = unique_temp_dir("md-ui-preview");
        let run_dir = root.join("runs/example-001");
        fs::create_dir_all(&run_dir).unwrap();
        fs::write(
            run_dir.join("energy.csv"),
            "step,time,kinetic,potential,total,temperature\n0,0.0000000000,1.0000000000,-2.0000000000,-1.0000000000,0.1000000000\n",
        )
        .unwrap();
        fs::write(
            run_dir.join("trajectory.xyz"),
            "2\nstep=0\nAr 0.0 0.0 0.0\nNe 1.0 0.0 0.0\n",
        )
        .unwrap();
        fs::write(run_dir.join("run-report.md"), "# Report\n").unwrap();

        let preview = load_run_preview(&run_dir);

        assert_eq!(preview.label, "example-001");
        assert_eq!(preview.energy.len(), 1);
        assert_eq!(preview.energy[0].total, -1.0);
        assert_eq!(preview.trajectory.as_ref().unwrap().frame_count, 1);
        assert_eq!(preview.trajectory.as_ref().unwrap().particle_count, 2);
        assert_eq!(preview.report.as_deref(), Some("# Report\n"));

        fs::remove_dir_all(root).unwrap();
    }

    fn unique_temp_dir(prefix: &str) -> PathBuf {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{}-{now}", std::process::id()))
    }
}
