#[cfg(feature = "desktop")]
use std::path::{Path, PathBuf};
#[cfg(feature = "desktop")]
use std::process::Command;
#[cfg(feature = "desktop")]
use std::sync::mpsc::{self, Receiver};
#[cfg(feature = "desktop")]
use std::thread;

#[cfg(feature = "desktop")]
use eframe::egui;
#[cfg(feature = "desktop")]
use egui_plot::{Line, Plot, PlotPoints};

#[cfg(feature = "desktop")]
use crate::{
    discover_run_dirs, discover_validation_configs, display_path, load_config_text,
    load_run_preview, save_config_text, RunPreview,
};

#[cfg(feature = "desktop")]
pub fn run() -> eframe::Result<()> {
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "MD Workstation",
        options,
        Box::new(|_cc| Ok(Box::new(MdUiApp::new()))),
    )
}

#[cfg(feature = "desktop")]
struct MdUiApp {
    workspace_root: PathBuf,
    validation_configs: Vec<PathBuf>,
    run_dirs: Vec<PathBuf>,
    config_path_text: String,
    config_text: String,
    selected_run: Option<PathBuf>,
    preview: Option<RunPreview>,
    status: String,
    run_output: String,
    run_receiver: Option<Receiver<RunResult>>,
}

#[cfg(feature = "desktop")]
struct RunResult {
    status: String,
    output: String,
}

#[cfg(feature = "desktop")]
impl MdUiApp {
    fn new() -> Self {
        let workspace_root = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let mut app = Self {
            workspace_root,
            validation_configs: Vec::new(),
            run_dirs: Vec::new(),
            config_path_text: String::new(),
            config_text: String::new(),
            selected_run: None,
            preview: None,
            status: "Ready".to_string(),
            run_output: String::new(),
            run_receiver: None,
        };
        app.refresh_workspace();
        if let Some(config_path) = app.validation_configs.first().cloned() {
            app.open_config(&config_path);
        }
        if let Some(run_dir) = app.run_dirs.last().cloned() {
            app.select_run(&run_dir);
        }
        app
    }

    fn refresh_workspace(&mut self) {
        match discover_validation_configs(&self.workspace_root) {
            Ok(configs) => self.validation_configs = configs,
            Err(error) => self.status = error.to_string(),
        }
        match discover_run_dirs(&self.workspace_root) {
            Ok(run_dirs) => self.run_dirs = run_dirs,
            Err(error) => self.status = error.to_string(),
        }
    }

    fn open_config(&mut self, config_path: &Path) {
        self.config_path_text = display_path(config_path);
        match load_config_text(config_path) {
            Ok(text) => {
                self.config_text = text;
                self.status = format!("Opened {}", config_path.display());
            }
            Err(error) => self.status = error.to_string(),
        }
    }

    fn save_config(&mut self) {
        let config_path = PathBuf::from(self.config_path_text.trim());
        match save_config_text(&config_path, &self.config_text) {
            Ok(()) => self.status = format!("Saved {}", config_path.display()),
            Err(error) => self.status = error.to_string(),
        }
    }

    fn select_run(&mut self, run_dir: &Path) {
        self.selected_run = Some(run_dir.to_path_buf());
        self.preview = Some(load_run_preview(run_dir));
    }

    fn start_run(&mut self) {
        if self.run_receiver.is_some() {
            return;
        }
        self.save_config();
        let config_path = PathBuf::from(self.config_path_text.trim());
        if config_path.as_os_str().is_empty() {
            self.status = "Config path is empty".to_string();
            return;
        }

        let workspace_root = self.workspace_root.clone();
        let (sender, receiver) = mpsc::channel();
        self.run_receiver = Some(receiver);
        self.status = format!("Running {}", config_path.display());
        self.run_output.clear();

        thread::spawn(move || {
            let output = Command::new("cargo")
                .args(["run", "-p", "md-cli", "--", "run"])
                .arg(&config_path)
                .current_dir(&workspace_root)
                .output();
            let result = match output {
                Ok(output) => {
                    let mut text = String::new();
                    text.push_str(&String::from_utf8_lossy(&output.stdout));
                    text.push_str(&String::from_utf8_lossy(&output.stderr));
                    let status = if output.status.success() {
                        "Run completed".to_string()
                    } else {
                        format!("Run failed: {}", output.status)
                    };
                    RunResult {
                        status,
                        output: text,
                    }
                }
                Err(error) => RunResult {
                    status: format!("Run failed: {error}"),
                    output: String::new(),
                },
            };
            let _ = sender.send(result);
        });
    }

    fn poll_run_receiver(&mut self) {
        let Some(receiver) = &self.run_receiver else {
            return;
        };
        match receiver.try_recv() {
            Ok(result) => {
                self.status = result.status;
                self.run_output = result.output;
                self.run_receiver = None;
                self.refresh_workspace();
                if let Some(run_dir) = self.run_dirs.last().cloned() {
                    self.select_run(&run_dir);
                }
            }
            Err(mpsc::TryRecvError::Empty) => {}
            Err(mpsc::TryRecvError::Disconnected) => {
                self.status = "Run worker disconnected".to_string();
                self.run_receiver = None;
            }
        }
    }
}

#[cfg(feature = "desktop")]
impl eframe::App for MdUiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_run_receiver();

        egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Workspace");
                ui.monospace(self.workspace_root.display().to_string());
                if ui.button("Refresh").clicked() {
                    self.refresh_workspace();
                }
                ui.separator();
                ui.label(&self.status);
            });
        });

        egui::SidePanel::left("selection_panel")
            .resizable(true)
            .default_width(260.0)
            .show(ctx, |ui| {
                ui.heading("Configs");
                egui::ScrollArea::vertical()
                    .id_source("configs")
                    .max_height(220.0)
                    .show(ui, |ui| {
                        for config_path in self.validation_configs.clone() {
                            let label = config_path
                                .file_name()
                                .and_then(|name| name.to_str())
                                .unwrap_or("config");
                            if ui.button(label).clicked() {
                                self.open_config(&config_path);
                            }
                        }
                    });

                ui.separator();
                ui.heading("Runs");
                egui::ScrollArea::vertical()
                    .id_source("runs")
                    .show(ui, |ui| {
                        for run_dir in self.run_dirs.clone() {
                            let selected = self.selected_run.as_ref() == Some(&run_dir);
                            let label = run_dir
                                .file_name()
                                .and_then(|name| name.to_str())
                                .unwrap_or("run");
                            if ui.selectable_label(selected, label).clicked() {
                                self.select_run(&run_dir);
                            }
                        }
                    });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.columns(2, |columns| {
                self.config_editor(&mut columns[0]);
                self.run_preview(&mut columns[1]);
            });
        });
    }
}

#[cfg(feature = "desktop")]
impl MdUiApp {
    fn config_editor(&mut self, ui: &mut egui::Ui) {
        ui.heading("Config");
        ui.horizontal(|ui| {
            ui.label("Path");
            ui.text_edit_singleline(&mut self.config_path_text);
        });
        ui.horizontal(|ui| {
            if ui.button("Open").clicked() {
                let path = PathBuf::from(self.config_path_text.trim());
                self.open_config(&path);
            }
            if ui.button("Save").clicked() {
                self.save_config();
            }
            let running = self.run_receiver.is_some();
            if ui.add_enabled(!running, egui::Button::new("Run")).clicked() {
                self.start_run();
            }
        });

        egui::ScrollArea::vertical()
            .id_source("config_editor")
            .show(ui, |ui| {
                ui.add(
                    egui::TextEdit::multiline(&mut self.config_text)
                        .desired_rows(32)
                        .code_editor(),
                );
            });

        if !self.run_output.is_empty() {
            ui.separator();
            ui.heading("Run Log");
            egui::ScrollArea::vertical()
                .id_source("run_log")
                .max_height(180.0)
                .show(ui, |ui| {
                    ui.monospace(&self.run_output);
                });
        }
    }

    fn run_preview(&mut self, ui: &mut egui::Ui) {
        ui.heading("Run Preview");
        let Some(preview) = &self.preview else {
            ui.label("No run selected");
            return;
        };

        ui.monospace(preview.path.display().to_string());
        if let Some(summary) = &preview.summary {
            ui.horizontal(|ui| {
                ui.label("Particles");
                ui.monospace(summary.particle_count.to_string());
                ui.label("Steps");
                ui.monospace(summary.steps.to_string());
                ui.label("Workflow");
                ui.monospace(&summary.workflow);
            });
            ui.horizontal(|ui| {
                ui.label("Final total");
                ui.monospace(format!("{:.10}", summary.final_energy.total));
                ui.label("Temperature");
                ui.monospace(format!("{:.6}", summary.final_energy.temperature));
            });
        }
        if let Some(manifest) = &preview.manifest {
            ui.horizontal(|ui| {
                ui.label("Status");
                ui.monospace(&manifest.status);
                ui.label("Engine");
                ui.monospace(manifest.engine_version.as_deref().unwrap_or("unknown"));
            });
        }
        for warning in &preview.warnings {
            ui.colored_label(egui::Color32::YELLOW, warning);
        }

        ui.separator();
        ui.heading("Trajectory");
        if let Some(trajectory) = &preview.trajectory {
            ui.horizontal(|ui| {
                ui.label("Frames");
                ui.monospace(trajectory.frame_count.to_string());
                ui.label("Particles");
                ui.monospace(trajectory.particle_count.to_string());
                ui.label("Steps");
                ui.monospace(format!(
                    "{}..{}",
                    trajectory
                        .first_step
                        .map(|step| step.to_string())
                        .unwrap_or_else(|| "?".to_string()),
                    trajectory
                        .last_step
                        .map(|step| step.to_string())
                        .unwrap_or_else(|| "?".to_string())
                ));
            });
            egui::ScrollArea::vertical()
                .id_source("trajectory_preview")
                .max_height(120.0)
                .show(ui, |ui| {
                    for atom in &trajectory.first_atoms {
                        ui.monospace(atom);
                    }
                });
        } else {
            ui.label("No trajectory.xyz");
        }

        ui.separator();
        ui.heading("Energy");
        if preview.energy.is_empty() {
            ui.label("No energy.csv");
        } else {
            let points =
                PlotPoints::from_iter(preview.energy.iter().map(|point| [point.time, point.total]));
            Plot::new("energy_plot")
                .height(220.0)
                .allow_scroll(false)
                .show(ui, |plot_ui| {
                    plot_ui.line(Line::new(points));
                });
        }

        ui.separator();
        ui.heading("Report");
        if let Some(report) = &preview.report {
            egui::ScrollArea::vertical()
                .id_source("report_view")
                .show(ui, |ui| {
                    ui.monospace(report);
                });
        } else {
            ui.label("No run-report.md");
        }
    }
}
