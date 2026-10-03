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
// A network: nodes (an index, a position, a label) joined by edges that carry one or
// more measures (correlation, fractional anisotropy, ...). This is the file-neutral
// model behind AFNI/FATCAT `Graph_Bucket` datasets (what `3dNetCorr`, `3dTrackID` and
// `ConvertDset -graphize` produce) and the matrix views in sumaru.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `afni-io` parses and writes the NIML `Graph_Bucket` file; its `adapt` module
//   converts to and from `Graph`. Nothing here knows about NIML.
// * `column::ColumnRange` is the range type (so a viewer's color scaling is the same
//   as for any other column) and `threshold::Threshold` selects edges, the same
//   thresholds that select surface nodes.
// * `tract.rs` holds the tracts that connect the nodes' regions; a network file may
//   link to one (that link is file business and stays in `afni-io`).
//
// HOW EDGES ARE STORED (matching `suma_datasets.c`)
//
// Edges are rows of a table, one value per measure, in one of four layouts:
//
//   * `Full`: an n x n matrix in COLUMN-major order. Row `e` is the matrix element
//     `(row = e % n, column = e / n)`. Direction is kept: `(r, c)` and `(c, r)` are
//     different rows.
//   * `LowerTriangle`: the part below the diagonal, column by column
//     (`(1,0) (2,0) .. (n-1,0) (2,1) ..`), `n (n - 1) / 2` rows. Symmetric: the
//     element `(c, r)` is the same edge.
//   * `LowerTriangleWithDiagonal`: the same including the diagonal, `n (n + 1) / 2`.
//   * `Sparse`: an explicit list of edges, each naming its two end nodes by NODE
//     INDEX (the first column of the node table), not by position, and carrying the
//     edge's own id. Direction is as listed.
//
// Positions are stored as the file has them, in AFNI's DICOM frame ("RAI"); use
// `positions_ras` or `domain::flip_dicom_ras` for the frame surfaces use.
// ---------------------------------------------------------------------------

//! A validated network with matrix, triangle and sparse edge layouts.

use std::collections::HashMap;

use crate::column::ColumnRange;
use crate::domain::flip_dicom_ras;
use crate::error::{Error, Result};
use crate::threshold::Threshold;

/// One node of the network.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphNode {
    /// The node's index as written in the file (usually 0 to n-1, but any distinct
    /// integers are allowed; sparse edges refer to these).
    pub index: i32,
    /// Where the node is, in AFNI's DICOM frame (+x left, +y posterior, +z up).
    pub position: [f32; 3],
    /// The node's name (an ROI label, for example).
    pub label: String,
}

/// One edge of a sparse graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseEdge {
    /// The edge's own id (the first column of the file's index list).
    pub id: i32,
    /// The index (not the position) of the first end node: the matrix row.
    pub row_node: i32,
    /// The index of the second end node: the matrix column.
    pub column_node: i32,
}

/// How the edge table is laid out.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeLayout {
    /// The whole `n x n` matrix, column-major.
    Full,
    /// The strict lower triangle, column by column.
    LowerTriangle,
    /// The lower triangle with the diagonal, column by column.
    LowerTriangleWithDiagonal,
    /// An explicit edge list.
    Sparse(Vec<SparseEdge>),
}

/// One measure: a name and one value per edge.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeMeasure {
    /// Name of the measure (the file's column label).
    pub label: String,
    /// One value per edge, in table order. `NaN` means missing.
    pub values: Vec<f32>,
}

/// A network, validated on construction.
#[derive(Debug, Clone, PartialEq)]
pub struct Graph {
    nodes: Vec<GraphNode>,
    layout: EdgeLayout,
    measures: Vec<EdgeMeasure>,
    /// For a sparse layout: each edge's end nodes as positions in `nodes`.
    resolved: Vec<(usize, usize)>,
}

/// Number of edges a layout holds for `n` nodes (`None` for a sparse layout, which
/// says for itself, or on overflow).
fn implied_edge_count(layout: &EdgeLayout, n: usize) -> Option<usize> {
    match layout {
        EdgeLayout::Full => n.checked_mul(n),
        EdgeLayout::LowerTriangle => n.checked_mul(n.saturating_sub(1)).map(|v| v / 2),
        EdgeLayout::LowerTriangleWithDiagonal => n.checked_mul(n + 1).map(|v| v / 2),
        EdgeLayout::Sparse(edges) => Some(edges.len()),
    }
}

impl Graph {
    /// Build a graph. Checks: at least one node and one measure; distinct node
    /// indices and finite positions; every measure has one value per edge; a sparse
    /// layout names only existing nodes.
    pub fn new(
        nodes: Vec<GraphNode>,
        layout: EdgeLayout,
        measures: Vec<EdgeMeasure>,
    ) -> Result<Self> {
        if nodes.is_empty() {
            return Err(Error::Empty("graph nodes".into()));
        }
        if measures.is_empty() {
            return Err(Error::Empty("graph measures".into()));
        }
        let mut by_index: HashMap<i32, usize> = HashMap::with_capacity(nodes.len());
        for (position, node) in nodes.iter().enumerate() {
            for c in node.position {
                crate::numeric::ensure_finite("node position", f64::from(c))?;
            }
            if by_index.insert(node.index, position).is_some() {
                return Err(Error::InvalidParameter {
                    name: "node index".into(),
                    reason: format!("index {} appears more than once", node.index),
                });
            }
        }
        let expected =
            implied_edge_count(&layout, nodes.len()).ok_or_else(|| Error::InvalidParameter {
                name: "node count".into(),
                reason: format!("{} nodes overflow the edge count", nodes.len()),
            })?;
        for m in &measures {
            if m.values.len() != expected {
                return Err(Error::LengthMismatch {
                    what: format!("measure '{}' values", m.label),
                    expected,
                    found: m.values.len(),
                });
            }
        }
        let resolved = match &layout {
            EdgeLayout::Sparse(edges) => edges
                .iter()
                .map(|e| {
                    let find = |index: i32| {
                        by_index
                            .get(&index)
                            .copied()
                            .ok_or_else(|| Error::InvalidParameter {
                                name: "edge end node".into(),
                                reason: format!(
                                    "edge {} names node {index}, which does not exist",
                                    e.id
                                ),
                            })
                    };
                    Ok((find(e.row_node)?, find(e.column_node)?))
                })
                .collect::<Result<Vec<_>>>()?,
            _ => Vec::new(),
        };
        Ok(Self {
            nodes,
            layout,
            measures,
            resolved,
        })
    }

    /// The nodes, in file order.
    pub fn nodes(&self) -> &[GraphNode] {
        &self.nodes
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// How the edges are laid out.
    pub fn layout(&self) -> &EdgeLayout {
        &self.layout
    }

    /// The measures, in file order.
    pub fn measures(&self) -> &[EdgeMeasure] {
        &self.measures
    }

    /// Number of edge rows.
    pub fn edge_count(&self) -> usize {
        self.measures[0].values.len()
    }

    /// Node positions in the RAS frame (see [`flip_dicom_ras`]).
    pub fn positions_ras(&self) -> Vec<[f32; 3]> {
        self.nodes
            .iter()
            .map(|n| flip_dicom_ras(n.position))
            .collect()
    }

    /// Edges in a lower-triangular layout before column `j`.
    fn triangle_prior(&self, j: usize) -> usize {
        let n = self.node_count();
        let diag = matches!(self.layout, EdgeLayout::LowerTriangleWithDiagonal);
        // sum over columns c < j of (n - 1 - c) (or n - c with the diagonal)
        let per = if diag { n } else { n - 1 };
        j * per - j * j.saturating_sub(1) / 2
    }

    /// The `(row, column)` node positions an edge row connects.
    pub fn edge_endpoints(&self, edge: usize) -> Option<(usize, usize)> {
        let n = self.node_count();
        if edge >= self.edge_count() {
            return None;
        }
        match &self.layout {
            EdgeLayout::Full => Some((edge % n, edge / n)),
            EdgeLayout::LowerTriangle | EdgeLayout::LowerTriangleWithDiagonal => {
                let diag = matches!(self.layout, EdgeLayout::LowerTriangleWithDiagonal);
                // Find the column: the last j whose first edge is at or before `edge`.
                let (mut lo, mut hi) = (0, n);
                while lo + 1 < hi {
                    let mid = (lo + hi) / 2;
                    if self.triangle_prior(mid) <= edge {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let offset = edge - self.triangle_prior(lo);
                let row = lo + offset + usize::from(!diag);
                Some((row, lo))
            }
            EdgeLayout::Sparse(_) => self.resolved.get(edge).copied(),
        }
    }

    /// The edge row for the matrix element `(row, column)` (positions in
    /// [`nodes`](Self::nodes)). Full layouts keep direction; triangular layouts treat
    /// `(r, c)` and `(c, r)` as one edge and have no diagonal unless stored; a sparse
    /// layout finds the first edge listed with exactly that order. `None` when the
    /// element is not stored. A sparse lookup scans the list; use
    /// [`matrix`](Self::matrix) for many.
    pub fn edge_index(&self, row: usize, column: usize) -> Option<usize> {
        let n = self.node_count();
        if row >= n || column >= n {
            return None;
        }
        match &self.layout {
            EdgeLayout::Full => Some(column * n + row),
            EdgeLayout::LowerTriangle | EdgeLayout::LowerTriangleWithDiagonal => {
                let diag = matches!(self.layout, EdgeLayout::LowerTriangleWithDiagonal);
                let (i, j) = (row.max(column), row.min(column));
                if !diag && i == j {
                    return None;
                }
                Some(self.triangle_prior(j) + (i - j) - usize::from(!diag))
            }
            EdgeLayout::Sparse(_) => self.resolved.iter().position(|&e| e == (row, column)),
        }
    }

    /// One value, or `None` if the element is not stored or the measure does not
    /// exist. A stored `NaN` is returned as is.
    pub fn value(&self, row: usize, column: usize, measure: usize) -> Option<f32> {
        let edge = self.edge_index(row, column)?;
        self.measures.get(measure)?.values.get(edge).copied()
    }

    /// One measure as a dense `n x n` matrix in display order (`[row * n + column]`).
    /// Triangular layouts fill both halves; a full layout is copied as it is; a
    /// sparse layout fills only the listed `(row, column)` cells (the first edge wins
    /// if one is listed twice). Cells with no edge are `None`.
    pub fn matrix(&self, measure: usize) -> Option<Vec<Option<f32>>> {
        let m = self.measures.get(measure)?;
        let n = self.node_count();
        let mut cells = vec![None; n * n];
        for (edge, &value) in m.values.iter().enumerate() {
            let (row, column) = self.edge_endpoints(edge)?;
            let mut put = |r: usize, c: usize| {
                if cells[r * n + c].is_none() {
                    cells[r * n + c] = Some(value);
                }
            };
            put(row, column);
            if !matches!(self.layout, EdgeLayout::Full | EdgeLayout::Sparse(_)) {
                put(column, row); // triangular layouts are symmetric
            }
        }
        Some(cells)
    }

    /// The finite extent of one measure (`NaN` and infinities are ignored). `None`
    /// if the measure does not exist or has no finite value.
    pub fn measure_range(&self, measure: usize) -> Option<ColumnRange> {
        let values = &self.measures.get(measure)?.values;
        let finite = values.iter().copied().filter(|v| v.is_finite());
        let (lo, hi) = finite.fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        });
        (lo <= hi)
            .then(|| ColumnRange::new(f64::from(lo), f64::from(hi)).ok())
            .flatten()
    }

    /// The edge rows whose value in `measure` passes `threshold` (a missing value
    /// never passes; see [`Threshold::passes`]). The same thresholds that select
    /// surface nodes select edges.
    pub fn edges_passing(&self, measure: usize, threshold: &Threshold) -> Result<Vec<usize>> {
        threshold.validate()?;
        let m = self.measures.get(measure).ok_or(Error::IndexOutOfRange {
            index: measure as i64,
            len: self.measures.len(),
        })?;
        Ok(m.values
            .iter()
            .enumerate()
            .filter(|(_, &v)| threshold.passes(f64::from(v)))
            .map(|(i, _)| i)
            .collect())
    }

    /// Sum of the finite values of `measure` over the edges touching each node (a
    /// node's "strength"). An edge counts for both of its ends; a diagonal element
    /// counts once.
    pub fn node_strength(&self, measure: usize) -> Option<Vec<f64>> {
        let m = self.measures.get(measure)?;
        let mut strength = vec![0.0; self.node_count()];
        for (edge, &v) in m.values.iter().enumerate() {
            if !v.is_finite() {
                continue;
            }
            let (a, b) = self.edge_endpoints(edge)?;
            strength[a] += f64::from(v);
            if a != b {
                strength[b] += f64::from(v);
            }
        }
        Some(strength)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nodes(n: usize) -> Vec<GraphNode> {
        (0..n)
            .map(|i| GraphNode {
                index: i as i32,
                position: [i as f32, 2.0 * i as f32, -(i as f32)],
                label: format!("n{i}"),
            })
            .collect()
    }

    fn measure(label: &str, values: Vec<f32>) -> EdgeMeasure {
        EdgeMeasure {
            label: label.into(),
            values,
        }
    }

    #[test]
    fn full_layout_is_column_major_and_directed() {
        // values chosen as 10*row + column of the matrix element they represent
        let n = 3;
        let mut v = vec![0.0; 9];
        for r in 0..n {
            for c in 0..n {
                v[c * n + r] = (10 * r + c) as f32;
            }
        }
        let g = Graph::new(nodes(n), EdgeLayout::Full, vec![measure("m", v)]).unwrap();
        assert_eq!(g.edge_count(), 9);
        assert_eq!(g.edge_endpoints(1), Some((1, 0)));
        assert_eq!(g.edge_endpoints(3), Some((0, 1)));
        assert_eq!(g.value(2, 1, 0), Some(21.0));
        assert_eq!(g.value(1, 2, 0), Some(12.0)); // direction kept
        let dense = g.matrix(0).unwrap();
        assert_eq!(dense[2 * n + 1], Some(21.0));
        assert_eq!(g.edge_endpoints(9), None);
    }

    #[test]
    fn triangular_layouts_follow_suma_packing() {
        // 4 nodes: lower triangle order (1,0) (2,0) (3,0) (2,1) (3,1) (3,2)
        let g = Graph::new(
            nodes(4),
            EdgeLayout::LowerTriangle,
            vec![measure("m", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0])],
        )
        .unwrap();
        let want = [(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2)];
        for (e, &(r, c)) in want.iter().enumerate() {
            assert_eq!(g.edge_endpoints(e), Some((r, c)), "edge {e}");
            assert_eq!(g.edge_index(r, c), Some(e));
            assert_eq!(g.edge_index(c, r), Some(e), "mirrored");
        }
        assert_eq!(g.edge_index(2, 2), None); // no diagonal
        assert_eq!(g.value(0, 3, 0), Some(3.0));
        let d = g.matrix(0).unwrap();
        assert_eq!(d[3], Some(3.0));
        assert_eq!(d[12], Some(3.0));
        assert_eq!(d[0], None);

        // With the diagonal: (0,0) (1,0) (2,0) (1,1) (2,1) (2,2) for 3 nodes.
        let g = Graph::new(
            nodes(3),
            EdgeLayout::LowerTriangleWithDiagonal,
            vec![measure("m", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0])],
        )
        .unwrap();
        let want = [(0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (2, 2)];
        for (e, &(r, c)) in want.iter().enumerate() {
            assert_eq!(g.edge_endpoints(e), Some((r, c)), "diag edge {e}");
            assert_eq!(g.edge_index(r, c), Some(e));
        }
        assert_eq!(g.value(1, 1, 0), Some(4.0));
    }

    #[test]
    fn triangular_endpoint_search_agrees_with_index_for_a_large_graph() {
        let n = 37;
        for layout in [
            EdgeLayout::LowerTriangle,
            EdgeLayout::LowerTriangleWithDiagonal,
        ] {
            let count = implied_edge_count(&layout, n).unwrap();
            let g = Graph::new(
                nodes(n),
                layout,
                vec![measure("m", (0..count).map(|i| i as f32).collect())],
            )
            .unwrap();
            for e in 0..count {
                let (r, c) = g.edge_endpoints(e).unwrap();
                assert_eq!(g.edge_index(r, c), Some(e));
            }
        }
    }

    #[test]
    fn sparse_edges_name_nodes_by_index_not_position() {
        // Node indices 5..8: an edge "5 -> 6" joins positions 0 and 1.
        let ns: Vec<GraphNode> = (0..4)
            .map(|i| GraphNode {
                index: 5 + i,
                position: [0.0; 3],
                label: format!("{}", 5 + i),
            })
            .collect();
        let edges = vec![
            SparseEdge {
                id: 0,
                row_node: 5,
                column_node: 6,
            },
            SparseEdge {
                id: 1,
                row_node: 6,
                column_node: 7,
            },
            SparseEdge {
                id: 2,
                row_node: 8,
                column_node: 5,
            },
        ];
        let g = Graph::new(
            ns,
            EdgeLayout::Sparse(edges),
            vec![
                measure("a", vec![1.0, 2.0, 3.0]),
                measure("b", vec![10.0, 20.0, 30.0]),
            ],
        )
        .unwrap();
        assert_eq!(g.edge_endpoints(2), Some((3, 0)));
        assert_eq!(g.value(0, 1, 1), Some(10.0));
        assert_eq!(g.value(1, 0, 0), None, "sparse keeps listed direction only");
        let d = g.matrix(0).unwrap();
        assert_eq!(d[12], Some(3.0));
        assert_eq!(d.iter().filter(|c| c.is_some()).count(), 3);
        // An edge naming a missing node is an error.
        let bad = Graph::new(
            nodes(2),
            EdgeLayout::Sparse(vec![SparseEdge {
                id: 0,
                row_node: 0,
                column_node: 9,
            }]),
            vec![measure("a", vec![1.0])],
        );
        assert!(bad.is_err());
    }

    #[test]
    fn construction_checks() {
        assert!(Graph::new(vec![], EdgeLayout::Full, vec![measure("m", vec![])]).is_err());
        assert!(Graph::new(nodes(2), EdgeLayout::Full, vec![]).is_err());
        // wrong number of values for the layout
        assert!(Graph::new(nodes(2), EdgeLayout::Full, vec![measure("m", vec![0.0; 3])]).is_err());
        assert!(Graph::new(
            nodes(3),
            EdgeLayout::LowerTriangle,
            vec![measure("m", vec![0.0; 4])]
        )
        .is_err());
        // duplicate node index, non-finite position
        let mut dup = nodes(2);
        dup[1].index = 0;
        assert!(Graph::new(dup, EdgeLayout::Full, vec![measure("m", vec![0.0; 4])]).is_err());
        let mut nan = nodes(2);
        nan[0].position[1] = f32::NAN;
        assert!(Graph::new(nan, EdgeLayout::Full, vec![measure("m", vec![0.0; 4])]).is_err());
    }

    #[test]
    fn ranges_thresholds_and_strength_reuse_core_types() {
        let g = Graph::new(
            nodes(3),
            EdgeLayout::LowerTriangle,
            vec![
                measure("r", vec![0.5, -0.8, f32::NAN]),
                measure("fa", vec![1.0, 2.0, 3.0]),
            ],
        )
        .unwrap();
        let range = g.measure_range(0).unwrap();
        assert_eq!((range.min, range.max), (-0.8f32 as f64, 0.5));
        assert!(g.measure_range(5).is_none());
        // |r| >= 0.6 passes only the -0.8 edge; NaN never passes.
        assert_eq!(
            g.edges_passing(0, &Threshold::AbsoluteAbove(0.6)).unwrap(),
            vec![1]
        );
        assert_eq!(
            g.edges_passing(1, &Threshold::Above(2.0)).unwrap(),
            vec![1, 2]
        );
        assert!(g.edges_passing(0, &Threshold::AbsoluteAbove(-1.0)).is_err());
        assert!(g.edges_passing(9, &Threshold::Off).is_err());
        // strength of fa: edges (1,0)=1, (2,0)=2, (2,1)=3
        assert_eq!(g.node_strength(1).unwrap(), vec![3.0, 4.0, 5.0]);
        // NaN is skipped: edges (1,0)=0.5 and (2,0)=-0.8 count, (2,1)=NaN does not.
        let s = g.node_strength(0).unwrap();
        assert!((s[0] - (0.5 - 0.8)).abs() < 1e-6);
        assert!((s[1] - 0.5).abs() < 1e-6 && (s[2] + 0.8).abs() < 1e-6);
    }

    #[test]
    fn positions_can_be_flipped_to_ras() {
        let g = Graph::new(nodes(2), EdgeLayout::Full, vec![measure("m", vec![0.0; 4])]).unwrap();
        assert_eq!(g.positions_ras()[1], [-1.0, -2.0, -1.0]);
        assert_eq!(
            g.nodes()[1].position,
            [1.0, 2.0, -1.0],
            "stored as in the file"
        );
    }
}
