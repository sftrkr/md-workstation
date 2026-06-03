#[cfg(feature = "desktop")]
fn main() -> eframe::Result<()> {
    md_ui::desktop::run()
}

#[cfg(not(feature = "desktop"))]
fn main() {
    eprintln!("md-ui desktop support is optional.");
    eprintln!("Run with: cargo run -p md-ui --features desktop");
}
