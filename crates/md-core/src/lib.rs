//! Core data structures for MD Workstation.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors produced by core state validation and state manipulation.
#[derive(Debug, Error)]
pub enum CoreError {
    #[error("state vectors must all have {expected} entries but {field} has {actual}")]
    LengthMismatch {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("{field}[{index}] is not finite: {value}")]
    NonFiniteValue {
        field: &'static str,
        index: usize,
        value: f64,
    },
    #[error("mass[{index}] must be positive and finite: {value}")]
    InvalidMass { index: usize, value: f64 },
    #[error("target temperature must be finite and non-negative: {0}")]
    InvalidTemperature(f64),
    #[error("cannot rescale zero kinetic energy to a positive temperature")]
    ZeroKineticEnergy,
    #[error("box dimensions must be positive and finite: ({x}, {y}, {z})")]
    InvalidBox { x: f64, y: f64, z: f64 },
}

/// Orthorhombic simulation box dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SimulationBox {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl SimulationBox {
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, CoreError> {
        let simulation_box = Self { x, y, z };
        simulation_box.validate()?;
        Ok(simulation_box)
    }

    pub fn validate(&self) -> Result<(), CoreError> {
        if self.x.is_finite()
            && self.y.is_finite()
            && self.z.is_finite()
            && self.x > 0.0
            && self.y > 0.0
            && self.z > 0.0
        {
            Ok(())
        } else {
            Err(CoreError::InvalidBox {
                x: self.x,
                y: self.y,
                z: self.z,
            })
        }
    }

    pub fn min_dimension(&self) -> f64 {
        self.x.min(self.y).min(self.z)
    }

    pub fn minimum_image_delta(&self, dx: f64, dy: f64, dz: f64) -> (f64, f64, f64) {
        (
            minimum_image_component(dx, self.x),
            minimum_image_component(dy, self.y),
            minimum_image_component(dz, self.z),
        )
    }

    pub fn wrap_position(&self, x: &mut f64, y: &mut f64, z: &mut f64) {
        *x = wrap_component(*x, self.x);
        *y = wrap_component(*y, self.y);
        *z = wrap_component(*z, self.z);
    }
}

fn minimum_image_component(delta: f64, length: f64) -> f64 {
    delta - length * (delta / length).round()
}

fn wrap_component(value: f64, length: f64) -> f64 {
    value.rem_euclid(length)
}

/// Structure-of-arrays simulation state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SystemState {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub vz: Vec<f64>,
    pub fx: Vec<f64>,
    pub fy: Vec<f64>,
    pub fz: Vec<f64>,
    pub mass: Vec<f64>,
    pub charge: Vec<f64>,
    pub element: Vec<String>,
    #[serde(default)]
    pub atom_name: Vec<Option<String>>,
    #[serde(default)]
    pub residue_name: Vec<Option<String>>,
    #[serde(default)]
    pub residue_id: Vec<Option<i32>>,
    #[serde(default)]
    pub chain_id: Vec<Option<String>>,
}

impl SystemState {
    pub fn new(particle_count: usize) -> Self {
        Self {
            x: vec![0.0; particle_count],
            y: vec![0.0; particle_count],
            z: vec![0.0; particle_count],
            vx: vec![0.0; particle_count],
            vy: vec![0.0; particle_count],
            vz: vec![0.0; particle_count],
            fx: vec![0.0; particle_count],
            fy: vec![0.0; particle_count],
            fz: vec![0.0; particle_count],
            mass: vec![1.0; particle_count],
            charge: vec![0.0; particle_count],
            element: vec!["X".to_string(); particle_count],
            atom_name: vec![None; particle_count],
            residue_name: vec![None; particle_count],
            residue_id: vec![None; particle_count],
            chain_id: vec![None; particle_count],
        }
    }

    pub fn particle_count(&self) -> usize {
        self.x.len()
    }

    pub fn clear_forces(&mut self) {
        self.fx.fill(0.0);
        self.fy.fill(0.0);
        self.fz.fill(0.0);
    }

    pub fn has_molecular_metadata(&self) -> bool {
        has_optional_string_metadata(&self.atom_name)
            || has_optional_string_metadata(&self.residue_name)
            || self.residue_id.iter().any(Option::is_some)
            || has_optional_string_metadata(&self.chain_id)
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.mass
            .iter()
            .zip(&self.vx)
            .zip(&self.vy)
            .zip(&self.vz)
            .map(|(((mass, vx), vy), vz)| 0.5 * mass * (vx * vx + vy * vy + vz * vz))
            .sum()
    }

    /// Temperature estimate in reduced units with k_B = 1.
    ///
    /// This simple Phase 2 estimate uses 3N degrees of freedom and does not
    /// subtract center-of-mass or constraint degrees of freedom.
    pub fn temperature(&self) -> f64 {
        let n = self.particle_count();
        if n == 0 {
            return 0.0;
        }

        2.0 * self.kinetic_energy() / (3.0 * n as f64)
    }

    pub fn remove_center_of_mass_velocity(&mut self) {
        let total_mass: f64 = self.mass.iter().sum();
        if total_mass == 0.0 {
            return;
        }

        let mut px = 0.0;
        let mut py = 0.0;
        let mut pz = 0.0;

        for i in 0..self.particle_count() {
            px += self.mass[i] * self.vx[i];
            py += self.mass[i] * self.vy[i];
            pz += self.mass[i] * self.vz[i];
        }

        let cm_vx = px / total_mass;
        let cm_vy = py / total_mass;
        let cm_vz = pz / total_mass;

        for i in 0..self.particle_count() {
            self.vx[i] -= cm_vx;
            self.vy[i] -= cm_vy;
            self.vz[i] -= cm_vz;
        }
    }

    pub fn rescale_temperature(&mut self, target_temperature: f64) -> Result<(), CoreError> {
        if !target_temperature.is_finite() || target_temperature < 0.0 {
            return Err(CoreError::InvalidTemperature(target_temperature));
        }

        if target_temperature == 0.0 {
            self.vx.fill(0.0);
            self.vy.fill(0.0);
            self.vz.fill(0.0);
            return Ok(());
        }

        let current_temperature = self.temperature();
        if current_temperature <= 0.0 || !current_temperature.is_finite() {
            return Err(CoreError::ZeroKineticEnergy);
        }

        let scale = (target_temperature / current_temperature).sqrt();
        for i in 0..self.particle_count() {
            self.vx[i] *= scale;
            self.vy[i] *= scale;
            self.vz[i] *= scale;
        }

        Ok(())
    }

    pub fn validate(&self) -> Result<(), CoreError> {
        let expected = self.x.len();
        for (field, actual) in [
            ("y", self.y.len()),
            ("z", self.z.len()),
            ("vx", self.vx.len()),
            ("vy", self.vy.len()),
            ("vz", self.vz.len()),
            ("fx", self.fx.len()),
            ("fy", self.fy.len()),
            ("fz", self.fz.len()),
            ("mass", self.mass.len()),
            ("charge", self.charge.len()),
            ("element", self.element.len()),
        ] {
            if actual != expected {
                return Err(CoreError::LengthMismatch {
                    field,
                    expected,
                    actual,
                });
            }
        }
        for (field, actual) in [
            ("atom_name", self.atom_name.len()),
            ("residue_name", self.residue_name.len()),
            ("residue_id", self.residue_id.len()),
            ("chain_id", self.chain_id.len()),
        ] {
            if actual != 0 && actual != expected {
                return Err(CoreError::LengthMismatch {
                    field,
                    expected,
                    actual,
                });
            }
        }

        validate_finite("x", &self.x)?;
        validate_finite("y", &self.y)?;
        validate_finite("z", &self.z)?;
        validate_finite("vx", &self.vx)?;
        validate_finite("vy", &self.vy)?;
        validate_finite("vz", &self.vz)?;
        validate_finite("fx", &self.fx)?;
        validate_finite("fy", &self.fy)?;
        validate_finite("fz", &self.fz)?;
        validate_finite("charge", &self.charge)?;

        for (index, mass) in self.mass.iter().copied().enumerate() {
            if !mass.is_finite() || mass <= 0.0 {
                return Err(CoreError::InvalidMass { index, value: mass });
            }
        }

        Ok(())
    }

    pub fn wrap_positions(&mut self, simulation_box: &SimulationBox) {
        for i in 0..self.particle_count() {
            simulation_box.wrap_position(&mut self.x[i], &mut self.y[i], &mut self.z[i]);
        }
    }
}

fn has_optional_string_metadata(values: &[Option<String>]) -> bool {
    values.iter().any(|value| {
        value
            .as_deref()
            .is_some_and(|value| !value.trim().is_empty())
    })
}

fn validate_finite(field: &'static str, values: &[f64]) -> Result<(), CoreError> {
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(CoreError::NonFiniteValue {
                field,
                index,
                value,
            });
        }
    }
    Ok(())
}

/// Energies sampled at a simulation step.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EnergySample {
    pub step: usize,
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
}

impl EnergySample {
    pub fn new(step: usize, time: f64, kinetic: f64, potential: f64, temperature: f64) -> Self {
        Self {
            step,
            time,
            kinetic,
            potential,
            total: kinetic + potential,
            temperature,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn system_state_validates_vector_lengths() {
        let mut state = SystemState::new(2);
        state.y.pop();

        let error = state.validate().unwrap_err();
        assert!(matches!(
            error,
            CoreError::LengthMismatch { field: "y", .. }
        ));
    }

    #[test]
    fn temperature_rescale_reaches_target() {
        let mut state = SystemState::new(2);
        state.vx[0] = 1.0;
        state.vx[1] = -1.0;

        state.rescale_temperature(0.5).unwrap();

        assert!((state.temperature() - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn center_of_mass_velocity_is_removed() {
        let mut state = SystemState::new(2);
        state.mass[0] = 1.0;
        state.mass[1] = 3.0;
        state.vx[0] = 4.0;
        state.vx[1] = 0.0;

        state.remove_center_of_mass_velocity();

        let total_px = state.mass[0] * state.vx[0] + state.mass[1] * state.vx[1];
        assert!(total_px.abs() < 1.0e-12);
    }

    #[test]
    fn invalid_box_dimensions_are_rejected() {
        let error = SimulationBox::new(1.0, 0.0, 1.0).unwrap_err();

        assert!(matches!(error, CoreError::InvalidBox { .. }));
    }

    #[test]
    fn non_finite_positions_are_rejected() {
        let mut state = SystemState::new(1);
        state.x[0] = f64::NAN;

        let error = state.validate().unwrap_err();

        assert!(matches!(
            error,
            CoreError::NonFiniteValue { field: "x", .. }
        ));
    }

    #[test]
    fn minimum_image_delta_uses_shortest_periodic_displacement() {
        let simulation_box = SimulationBox::new(10.0, 8.0, 6.0).unwrap();

        let (dx, dy, dz) = simulation_box.minimum_image_delta(9.0, -5.0, 2.5);

        assert!((dx + 1.0).abs() < 1.0e-12);
        assert!((dy - 3.0).abs() < 1.0e-12);
        assert!((dz - 2.5).abs() < 1.0e-12);
    }

    #[test]
    fn wrap_positions_keeps_particles_inside_box() {
        let simulation_box = SimulationBox::new(10.0, 8.0, 6.0).unwrap();
        let mut state = SystemState::new(2);
        state.x[0] = -0.25;
        state.y[0] = 8.25;
        state.z[0] = 12.1;
        state.x[1] = 10.0;
        state.y[1] = -8.0;
        state.z[1] = -0.1;

        state.wrap_positions(&simulation_box);

        assert!((state.x[0] - 9.75).abs() < 1.0e-12);
        assert!((state.y[0] - 0.25).abs() < 1.0e-12);
        assert!((state.z[0] - 0.1).abs() < 1.0e-12);
        assert!(state.x[1].abs() < 1.0e-12);
        assert!(state.y[1].abs() < 1.0e-12);
        assert!((state.z[1] - 5.9).abs() < 1.0e-12);
    }
}
