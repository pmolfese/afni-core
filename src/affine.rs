// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, a United States Government Work. See the
// crate LICENSE for the full public-domain and CC0 notice.

//! Validated spatial affine transformations.
//!
//! AFNI commonly stores a spatial affine as twelve row-major numbers: a 3x3
//! linear transform followed by a translation column. This module owns the
//! file-neutral meaning and mathematics of those transforms. Reading and
//! writing `.aff12.1D` files belongs in `afni-io`.
//!
//! A matrix maps an input point to an output point:
//!
//! ```text
//! output = linear * input + translation
//! ```
//!
//! The numbers do not encode names such as "base" or "source". Callers must
//! retain that direction explicitly in their own variables or metadata. The
//! coordinate-axis convention *is* recorded here, because silently mixing
//! AFNI DICOM/RAI coordinates with NIfTI RAS coordinates changes the result.

use crate::{Error, Result};

/// A row-major homogeneous 4x4 affine matrix.
pub type AffineMatrix = [[f64; 4]; 4];

/// World-coordinate axis convention used by an affine transform.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CoordinateConvention {
    /// AFNI's native DICOM/RAI convention.
    ///
    /// Converting its coordinates to NIfTI RAS negates x and y.
    AfniDICOM,
    /// NIfTI's common x=Right, y=Anterior, z=Superior (RAS) convention.
    Ras,
}

/// One validated, invertible spatial affine transformation.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineTransform {
    matrix: AffineMatrix,
    convention: CoordinateConvention,
}

impl AffineTransform {
    /// Construct an affine from a homogeneous 4x4 matrix.
    ///
    /// Every entry must be finite, the final row must be `[0, 0, 0, 1]`, and
    /// the 3x3 linear part must be invertible. A nearly singular matrix is
    /// rejected rather than producing unstable inverse coordinates.
    pub fn new(matrix: AffineMatrix, convention: CoordinateConvention) -> Result<Self> {
        validate_matrix(&matrix)?;
        Ok(Self { matrix, convention })
    }

    /// Construct an AFNI DICOM/RAI transform from one `aff12` row.
    pub fn from_aff12_row(values: [f64; 12]) -> Result<Self> {
        Self::from_aff12_row_in(values, CoordinateConvention::AfniDICOM)
    }

    /// Construct a transform from twelve row-major values in `convention`.
    pub fn from_aff12_row_in(values: [f64; 12], convention: CoordinateConvention) -> Result<Self> {
        Self::new(
            [
                [values[0], values[1], values[2], values[3]],
                [values[4], values[5], values[6], values[7]],
                [values[8], values[9], values[10], values[11]],
                [0.0, 0.0, 0.0, 1.0],
            ],
            convention,
        )
    }

    /// Identity transform in the requested coordinate convention.
    pub fn identity(convention: CoordinateConvention) -> Self {
        Self {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            convention,
        }
    }

    /// The row-major homogeneous matrix.
    pub fn matrix(&self) -> &AffineMatrix {
        &self.matrix
    }

    /// The coordinate-axis convention in which this matrix is expressed.
    pub fn convention(&self) -> CoordinateConvention {
        self.convention
    }

    /// Return the twelve values used by AFNI's 3x4 textual representation.
    pub fn to_aff12_row(&self) -> [f64; 12] {
        [
            self.matrix[0][0],
            self.matrix[0][1],
            self.matrix[0][2],
            self.matrix[0][3],
            self.matrix[1][0],
            self.matrix[1][1],
            self.matrix[1][2],
            self.matrix[1][3],
            self.matrix[2][0],
            self.matrix[2][1],
            self.matrix[2][2],
            self.matrix[2][3],
        ]
    }

    /// Apply the affine to a point, including its translation.
    pub fn apply_point(&self, point: [f64; 3]) -> Result<[f64; 3]> {
        require_finite_triplet("affine input point", point)?;
        let result = std::array::from_fn(|row| {
            self.matrix[row][0] * point[0]
                + self.matrix[row][1] * point[1]
                + self.matrix[row][2] * point[2]
                + self.matrix[row][3]
        });
        require_finite_triplet("affine output point", result)?;
        Ok(result)
    }

    /// Apply only the 3x3 linear part to a vector, ignoring translation.
    pub fn apply_vector(&self, vector: [f64; 3]) -> Result<[f64; 3]> {
        require_finite_triplet("affine input vector", vector)?;
        let result = std::array::from_fn(|row| {
            self.matrix[row][0] * vector[0]
                + self.matrix[row][1] * vector[1]
                + self.matrix[row][2] * vector[2]
        });
        require_finite_triplet("affine output vector", result)?;
        Ok(result)
    }

    /// Return the inverse transformation.
    pub fn inverse(&self) -> Result<Self> {
        let linear_inverse = inverse_linear_part(&self.matrix)?;
        let translation = [self.matrix[0][3], self.matrix[1][3], self.matrix[2][3]];
        let inverse_translation: [f64; 3] = std::array::from_fn(|row| {
            -(linear_inverse[row][0] * translation[0]
                + linear_inverse[row][1] * translation[1]
                + linear_inverse[row][2] * translation[2])
        });
        Self::new(
            [
                [
                    linear_inverse[0][0],
                    linear_inverse[0][1],
                    linear_inverse[0][2],
                    inverse_translation[0],
                ],
                [
                    linear_inverse[1][0],
                    linear_inverse[1][1],
                    linear_inverse[1][2],
                    inverse_translation[1],
                ],
                [
                    linear_inverse[2][0],
                    linear_inverse[2][1],
                    linear_inverse[2][2],
                    inverse_translation[2],
                ],
                [0.0, 0.0, 0.0, 1.0],
            ],
            self.convention,
        )
    }

    /// Compose transforms in application order: apply `self`, then `next`.
    ///
    /// This matches `cat_matvec A B`: the transform represented by `B` follows
    /// the transform represented by `A`, so the matrix product is `B * A`.
    pub fn then(&self, next: &Self) -> Result<Self> {
        require_same_convention(self.convention, next.convention)?;
        Self::new(
            multiply_matrices(&next.matrix, &self.matrix),
            self.convention,
        )
    }

    /// Express the same physical mapping in another axis convention.
    ///
    /// AFNI DICOM/RAI and NIfTI RAS differ by negating x and y. Both the input
    /// and output coordinate
    /// bases are converted, so the result is `F * M * F`, where
    /// `F = diag(-1, -1, 1, 1)`.
    pub fn convert_coordinates(&self, target: CoordinateConvention) -> Result<Self> {
        if target == self.convention {
            return Ok(self.clone());
        }
        let flip = [
            [-1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let converted = multiply_matrices(&flip, &multiply_matrices(&self.matrix, &flip));
        Self::new(converted, target)
    }
}

/// A non-empty sequence of affine transforms, such as one motion transform
/// per volume in a time-series dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineTransformSeries {
    transforms: Vec<AffineTransform>,
    convention: CoordinateConvention,
}

impl AffineTransformSeries {
    /// Build a series, requiring at least one transform and one convention.
    pub fn new(transforms: Vec<AffineTransform>) -> Result<Self> {
        let convention = transforms
            .first()
            .ok_or_else(|| Error::Empty("affine transform series".into()))?
            .convention();
        for transform in &transforms {
            require_same_convention(convention, transform.convention())?;
        }
        Ok(Self {
            transforms,
            convention,
        })
    }

    /// Number of transforms in the series.
    pub fn len(&self) -> usize {
        self.transforms.len()
    }

    /// Whether this series has no transforms.
    ///
    /// A constructed series is never empty; this method is provided alongside
    /// [`len`](Self::len) for ordinary collection-style code.
    pub fn is_empty(&self) -> bool {
        self.transforms.is_empty()
    }

    /// Coordinate-axis convention shared by every transform.
    pub fn convention(&self) -> CoordinateConvention {
        self.convention
    }

    /// All transforms in file or time-point order.
    pub fn as_slice(&self) -> &[AffineTransform] {
        &self.transforms
    }

    /// Retrieve one transform by zero-based index.
    pub fn get(&self, index: usize) -> Option<&AffineTransform> {
        self.transforms.get(index)
    }

    /// Apply one indexed transform to a point.
    pub fn apply_point(&self, index: usize, point: [f64; 3]) -> Result<[f64; 3]> {
        let transform = self.transforms.get(index).ok_or(Error::IndexOutOfRange {
            index: i64::try_from(index).unwrap_or(i64::MAX),
            len: self.len(),
        })?;
        transform.apply_point(point)
    }

    /// Invert every transform without changing series order.
    pub fn inverse(&self) -> Result<Self> {
        Self::new(
            self.transforms
                .iter()
                .map(AffineTransform::inverse)
                .collect::<Result<_>>()?,
        )
    }

    /// Apply one static transform after every transform in the series.
    pub fn then(&self, next: &AffineTransform) -> Result<Self> {
        Self::new(
            self.transforms
                .iter()
                .map(|transform| transform.then(next))
                .collect::<Result<_>>()?,
        )
    }

    /// Pairwise composition in application order.
    ///
    /// Both series must have the same length. For each index, the transform in
    /// `self` is applied first and the corresponding transform in `next`
    /// follows it.
    pub fn then_series(&self, next: &Self) -> Result<Self> {
        if self.len() != next.len() {
            return Err(Error::LengthMismatch {
                what: "affine transform series composition".into(),
                expected: self.len(),
                found: next.len(),
            });
        }
        Self::new(
            self.transforms
                .iter()
                .zip(&next.transforms)
                .map(|(first, second)| first.then(second))
                .collect::<Result<_>>()?,
        )
    }

    /// Express every transform in another coordinate-axis convention.
    pub fn convert_coordinates(&self, target: CoordinateConvention) -> Result<Self> {
        Self::new(
            self.transforms
                .iter()
                .map(|transform| transform.convert_coordinates(target))
                .collect::<Result<_>>()?,
        )
    }
}

fn require_same_convention(
    first: CoordinateConvention,
    second: CoordinateConvention,
) -> Result<()> {
    if first == second {
        Ok(())
    } else {
        Err(Error::InvalidParameter {
            name: "affine coordinate convention".into(),
            reason: format!("cannot combine {first:?} with {second:?}"),
        })
    }
}

fn require_finite_triplet(what: &str, values: [f64; 3]) -> Result<()> {
    for value in values {
        if !value.is_finite() {
            return Err(Error::NonFinite {
                what: what.into(),
                value,
            });
        }
    }
    Ok(())
}

fn validate_matrix(matrix: &AffineMatrix) -> Result<()> {
    for (row, values) in matrix.iter().enumerate() {
        for (column, &value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(Error::NonFinite {
                    what: format!("affine matrix entry [{row}][{column}]"),
                    value,
                });
            }
        }
    }
    if matrix[3] != [0.0, 0.0, 0.0, 1.0] {
        return Err(Error::InvalidParameter {
            name: "affine matrix".into(),
            reason: "last row must be exactly [0, 0, 0, 1]".into(),
        });
    }
    inverse_linear_part(matrix).map(|_| ())
}

fn inverse_linear_part(matrix: &AffineMatrix) -> Result<[[f64; 3]; 3]> {
    let scale = matrix
        .iter()
        .take(3)
        .flat_map(|row| row.iter().take(3))
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Err(singular_error());
    }

    // Work on a scaled matrix so the determinant test is independent of the
    // transform's units and cannot overflow merely because entries are large.
    let a = matrix[0][0] / scale;
    let b = matrix[0][1] / scale;
    let c = matrix[0][2] / scale;
    let d = matrix[1][0] / scale;
    let e = matrix[1][1] / scale;
    let f = matrix[1][2] / scale;
    let g = matrix[2][0] / scale;
    let h = matrix[2][1] / scale;
    let i = matrix[2][2] / scale;
    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if !det.is_finite() || det.abs() <= 64.0 * f64::EPSILON {
        return Err(singular_error());
    }

    let inverse_scale = 1.0 / (det * scale);
    let inverse = [
        [
            (e * i - f * h) * inverse_scale,
            (c * h - b * i) * inverse_scale,
            (b * f - c * e) * inverse_scale,
        ],
        [
            (f * g - d * i) * inverse_scale,
            (a * i - c * g) * inverse_scale,
            (c * d - a * f) * inverse_scale,
        ],
        [
            (d * h - e * g) * inverse_scale,
            (b * g - a * h) * inverse_scale,
            (a * e - b * d) * inverse_scale,
        ],
    ];
    if inverse.iter().flatten().all(|value| value.is_finite()) {
        Ok(inverse)
    } else {
        Err(singular_error())
    }
}

fn singular_error() -> Error {
    Error::InvalidParameter {
        name: "affine matrix".into(),
        reason: "3x3 linear part is singular or numerically degenerate".into(),
    }
}

fn multiply_matrices(left: &AffineMatrix, right: &AffineMatrix) -> AffineMatrix {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            (0..4)
                .map(|inner| left[row][inner] * right[inner][column])
                .sum()
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: [f64; 3], expected: [f64; 3]) {
        for axis in 0..3 {
            assert!((actual[axis] - expected[axis]).abs() < 1.0e-12);
        }
    }

    fn translation(x: f64, y: f64, z: f64) -> AffineTransform {
        AffineTransform::from_aff12_row([1.0, 0.0, 0.0, x, 0.0, 1.0, 0.0, y, 0.0, 0.0, 1.0, z])
            .unwrap()
    }

    #[test]
    fn aff12_round_trip_apply_and_inverse() {
        let row = [
            2.0, 0.5, 0.0, 10.0, 0.0, 3.0, 0.25, -4.0, 0.0, 0.0, 4.0, 7.0,
        ];
        let transform = AffineTransform::from_aff12_row(row).unwrap();
        assert_eq!(transform.to_aff12_row(), row);
        let point = [1.0, 2.0, 3.0];
        let moved = transform.apply_point(point).unwrap();
        close(moved, [13.0, 2.75, 19.0]);
        close(
            transform.inverse().unwrap().apply_point(moved).unwrap(),
            point,
        );
    }

    #[test]
    fn composition_is_in_application_order() {
        let shift = translation(10.0, 0.0, 0.0);
        let scale = AffineTransform::from_aff12_row([
            2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0,
        ])
        .unwrap();
        close(
            shift
                .then(&scale)
                .unwrap()
                .apply_point([1.0, 0.0, 0.0])
                .unwrap(),
            [22.0, 0.0, 0.0],
        );
    }

    #[test]
    fn rai_ras_conversion_changes_both_coordinate_bases() {
        let rai = AffineTransform::from_aff12_row([
            1.0, 0.0, 2.0, 3.0, 0.0, 1.0, 0.0, 4.0, 5.0, 0.0, 1.0, 6.0,
        ])
        .unwrap();
        let ras = rai.convert_coordinates(CoordinateConvention::Ras).unwrap();
        assert_eq!(
            ras.to_aff12_row(),
            [1.0, 0.0, -2.0, -3.0, 0.0, 1.0, 0.0, -4.0, -5.0, 0.0, 1.0, 6.0]
        );
        assert_eq!(
            ras.convert_coordinates(CoordinateConvention::AfniDICOM)
                .unwrap(),
            rai
        );
    }

    #[test]
    fn invalid_matrices_are_rejected() {
        assert!(AffineTransform::from_aff12_row([0.0; 12]).is_err());
        let mut non_affine = AffineTransform::identity(CoordinateConvention::Ras)
            .matrix()
            .to_owned();
        non_affine[3][0] = 1.0;
        assert!(AffineTransform::new(non_affine, CoordinateConvention::Ras).is_err());
        let mut non_finite = translation(0.0, 0.0, 0.0).to_aff12_row();
        non_finite[3] = f64::NAN;
        assert!(AffineTransform::from_aff12_row(non_finite).is_err());
    }

    #[test]
    fn series_supports_indexing_inversion_and_composition() {
        let series = AffineTransformSeries::new(vec![
            translation(1.0, 0.0, 0.0),
            translation(2.0, 0.0, 0.0),
        ])
        .unwrap();
        assert_eq!(series.len(), 2);
        close(
            series.apply_point(1, [3.0, 0.0, 0.0]).unwrap(),
            [5.0, 0.0, 0.0],
        );
        close(
            series
                .inverse()
                .unwrap()
                .apply_point(1, [5.0, 0.0, 0.0])
                .unwrap(),
            [3.0, 0.0, 0.0],
        );
        let shifted = series.then(&translation(10.0, 0.0, 0.0)).unwrap();
        close(
            shifted.apply_point(0, [0.0, 0.0, 0.0]).unwrap(),
            [11.0, 0.0, 0.0],
        );
        assert!(series.apply_point(2, [0.0; 3]).is_err());
        assert!(AffineTransformSeries::new(Vec::new()).is_err());
    }
}
