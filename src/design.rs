// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, a United States Government Work. See the
// LICENSE file at the crate root.

//! File-neutral regression design matrices.
//!
//! A design matrix has one stored row per retained observation and one column
//! per regressor.  The optional observation map relates those rows to the full
//! pre-censoring timeline.  File syntax such as AFNI's comment-wrapped
//! `X.xmat.1D` header belongs in `afni-io`; this module owns the validated
//! meaning that other formats and generated matrices can share.

use crate::{Error, Result};

/// The broad purpose of a design-matrix regressor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegressorRole {
    /// A baseline term such as a run polynomial.
    Baseline,
    /// A nuisance term such as motion or physiological noise.
    Nuisance,
    /// A regressor of interest, commonly a task stimulus.
    Interest,
    /// The producer did not describe the regressor's purpose.
    Unknown,
}

/// A named design-matrix column and its grouping metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Regressor {
    label: String,
    role: RegressorRole,
    group: Option<i32>,
}

impl Regressor {
    /// Create a regressor.
    pub fn new(label: impl Into<String>, role: RegressorRole, group: Option<i32>) -> Result<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::Empty("design-matrix regressor label".into()));
        }
        Ok(Self { label, role, group })
    }

    /// The column label.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// The regressor's broad purpose.
    pub fn role(&self) -> RegressorRole {
        self.role
    }

    /// An optional producer-defined group number.
    pub fn group(&self) -> Option<i32> {
        self.group
    }
}

/// A named inclusive range of design-matrix columns belonging to one stimulus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stimulus {
    label: String,
    first_regressor: usize,
    last_regressor: usize,
}

impl Stimulus {
    /// Create a stimulus covering `first_regressor..=last_regressor`.
    pub fn new(
        label: impl Into<String>,
        first_regressor: usize,
        last_regressor: usize,
    ) -> Result<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::Empty("stimulus label".into()));
        }
        if first_regressor > last_regressor {
            return Err(Error::InvalidParameter {
                name: "stimulus regressor range".into(),
                reason: format!("{first_regressor}..={last_regressor} is reversed"),
            });
        }
        Ok(Self {
            label,
            first_regressor,
            last_regressor,
        })
    }

    /// The stimulus label.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// First column belonging to the stimulus.
    pub fn first_regressor(&self) -> usize {
        self.first_regressor
    }

    /// Last column belonging to the stimulus, inclusive.
    pub fn last_regressor(&self) -> usize {
        self.last_regressor
    }
}

/// A validated, row-major regression design matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignMatrix {
    rows: usize,
    cols: usize,
    values: Vec<f64>,
    regressors: Vec<Regressor>,
    full_observations: usize,
    retained_observations: Vec<usize>,
    run_starts: Vec<usize>,
    sample_period_seconds: Option<f64>,
    stimuli: Vec<Stimulus>,
}

impl DesignMatrix {
    /// Build and validate a design matrix.
    ///
    /// `retained_observations[row]` is the corresponding row in the full,
    /// uncensored timeline. It must be strictly increasing. `run_starts` uses
    /// that same full-timeline coordinate system and must begin at zero.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rows: usize,
        cols: usize,
        values: Vec<f64>,
        regressors: Vec<Regressor>,
        full_observations: usize,
        retained_observations: Vec<usize>,
        run_starts: Vec<usize>,
        sample_period_seconds: Option<f64>,
        stimuli: Vec<Stimulus>,
    ) -> Result<Self> {
        if rows == 0 {
            return Err(Error::Empty("design-matrix observations".into()));
        }
        if cols == 0 {
            return Err(Error::Empty("design-matrix regressors".into()));
        }
        let expected = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::InvalidParameter {
                name: "design-matrix dimensions".into(),
                reason: format!("{rows} * {cols} overflows usize"),
            })?;
        if values.len() != expected {
            return Err(Error::LengthMismatch {
                what: "design-matrix values".into(),
                expected,
                found: values.len(),
            });
        }
        if let Some(&value) = values.iter().find(|value| !value.is_finite()) {
            return Err(Error::NonFinite {
                what: "design-matrix value".into(),
                value,
            });
        }
        if regressors.len() != cols {
            return Err(Error::LengthMismatch {
                what: "design-matrix regressors".into(),
                expected: cols,
                found: regressors.len(),
            });
        }
        if full_observations == 0 {
            return Err(Error::Empty("full observation timeline".into()));
        }
        if retained_observations.len() != rows {
            return Err(Error::LengthMismatch {
                what: "retained-observation map".into(),
                expected: rows,
                found: retained_observations.len(),
            });
        }
        validate_sorted_indices(
            "retained-observation map",
            &retained_observations,
            full_observations,
        )?;
        if run_starts.first() != Some(&0) {
            return Err(Error::InvalidParameter {
                name: "run starts".into(),
                reason: "the first run must start at observation 0".into(),
            });
        }
        validate_sorted_indices("run starts", &run_starts, full_observations)?;
        if let Some(period) = sample_period_seconds {
            if !period.is_finite() || period <= 0.0 {
                return Err(Error::InvalidParameter {
                    name: "sample period".into(),
                    reason: format!("must be finite and positive, got {period}"),
                });
            }
        }
        for stimulus in &stimuli {
            if stimulus.last_regressor >= cols {
                return Err(Error::ColumnIndexOutOfRange {
                    index: stimulus.last_regressor,
                    len: cols,
                });
            }
        }
        Ok(Self {
            rows,
            cols,
            values,
            regressors,
            full_observations,
            retained_observations,
            run_starts,
            sample_period_seconds,
            stimuli,
        })
    }

    /// Number of retained observations (stored rows).
    pub fn observation_count(&self) -> usize {
        self.rows
    }

    /// Number of regressors (stored columns).
    pub fn regressor_count(&self) -> usize {
        self.cols
    }

    /// Number of observations before censoring.
    pub fn full_observation_count(&self) -> usize {
        self.full_observations
    }

    /// All values in row-major order.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Value at a retained observation and regressor, if both are in range.
    pub fn get(&self, observation: usize, regressor: usize) -> Option<f64> {
        if observation >= self.rows || regressor >= self.cols {
            return None;
        }
        self.values
            .get(observation * self.cols + regressor)
            .copied()
    }

    /// Regressor descriptions in column order.
    pub fn regressors(&self) -> &[Regressor] {
        &self.regressors
    }

    /// Map each stored row to the full, uncensored observation timeline.
    pub fn retained_observations(&self) -> &[usize] {
        &self.retained_observations
    }

    /// Whether any full-timeline observations were censored.
    pub fn is_censored(&self) -> bool {
        self.rows != self.full_observations
            || self
                .retained_observations
                .iter()
                .enumerate()
                .any(|(row, &full)| row != full)
    }

    /// Full-timeline indices at which runs begin.
    pub fn run_starts(&self) -> &[usize] {
        &self.run_starts
    }

    /// Sample period in seconds, when known.
    pub fn sample_period_seconds(&self) -> Option<f64> {
        self.sample_period_seconds
    }

    /// Named stimulus-to-regressor ranges.
    pub fn stimuli(&self) -> &[Stimulus] {
        &self.stimuli
    }
}

fn validate_sorted_indices(what: &str, indices: &[usize], len: usize) -> Result<()> {
    if indices.is_empty() {
        return Err(Error::Empty(what.into()));
    }
    for (position, &index) in indices.iter().enumerate() {
        if index >= len {
            return Err(Error::IndexOutOfRange {
                index: i64::try_from(index).unwrap_or(i64::MAX),
                len,
            });
        }
        if position > 0 && indices[position - 1] >= index {
            return Err(Error::InvalidParameter {
                name: what.into(),
                reason: "indices must be strictly increasing".into(),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn regressor(label: &str, role: RegressorRole, group: i32) -> Regressor {
        Regressor::new(label, role, Some(group)).unwrap()
    }

    #[test]
    fn validates_and_exposes_a_censored_multirun_design() {
        let design = DesignMatrix::new(
            3,
            2,
            vec![1.0, 10.0, 1.0, 20.0, 1.0, 40.0],
            vec![
                regressor("Run#1Pol#0", RegressorRole::Baseline, -1),
                regressor("visual#0", RegressorRole::Interest, 1),
            ],
            5,
            vec![0, 1, 3],
            vec![0, 3],
            Some(2.0),
            vec![Stimulus::new("visual", 1, 1).unwrap()],
        )
        .unwrap();

        assert_eq!(design.get(2, 1), Some(40.0));
        assert_eq!(design.full_observation_count(), 5);
        assert_eq!(design.retained_observations(), [0, 1, 3]);
        assert!(design.is_censored());
        assert_eq!(design.run_starts(), [0, 3]);
        assert_eq!(design.stimuli()[0].label(), "visual");
    }

    #[test]
    fn rejects_inconsistent_semantics() {
        let regs = vec![regressor("constant", RegressorRole::Baseline, -1)];
        assert!(DesignMatrix::new(
            2,
            1,
            vec![1.0, 1.0],
            regs.clone(),
            3,
            vec![0, 0],
            vec![0],
            Some(1.0),
            vec![],
        )
        .is_err());
        assert!(DesignMatrix::new(
            2,
            1,
            vec![1.0, 1.0],
            regs,
            2,
            vec![0, 1],
            vec![1],
            Some(1.0),
            vec![],
        )
        .is_err());
    }
}
