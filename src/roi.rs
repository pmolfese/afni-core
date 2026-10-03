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
// The semantic model of a region of interest (ROI) on a surface: who it belongs
// to, what it is called, how it is drawn on screen, how it was made (an open path,
// a closed path, a filled area, or just a collection of nodes), and the ordered
// strokes that make it up. It is the file-neutral version of SUMA's "drawn ROI"
// and of sumaru's `roi.rs`.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `afni-io` reads and writes `.niml.roi` files (`NodeRoi`); its `adapt` module
//   converts them to and from the types here. Nothing in this file touches a file.
// * `NodeSet` is the canonical form of "which nodes": sorted, no repeats. The
//   operations in `roi_ops.rs` (boundary, grow, shrink, fill, shortest path,
//   components) work on node sets and a `SurfaceTopology`/`SurfaceMesh`.
// * `roi_edit.rs` has the edit commands (append a stroke, join the ends, fill,
//   relabel) and an undo/redo editor. Mouse gestures stay in the viewer.
// * `Roi::to_dataset` and `rois_to_dataset` produce a `Dataset` (kind `Roi`,
//   sparse rows, one integer label column), what SUMA's `ROI2dataset` writes.
// * `labels::LabelEntry` is what `Roi::label_entry` returns, so a set of ROIs can
//   become a label table.
//
// LOSSLESSNESS
//
// SUMA writes its codes as integers (`Type`, the element type and action of each
// stroke) and the hemisphere as text. Files in the wild contain codes SUMA itself
// does not define (sumaru writes `Type="4"`). Every enum here has an `Other(..)`
// variant that keeps an unknown code, so reading and writing a file never loses it.
//
// WHICH NODES ARE IN AN ROI
//
// Like SUMA (`SUMA_NodesInROI`), an ROI's nodes are the nodes of all its strokes,
// whatever their element kind. `ordered_nodes` keeps the order they were drawn in,
// dropping a node that repeats the last node of the previous stroke (the junction
// of two strokes); `node_set` sorts them and removes repeats. `drawn_nodes` and
// `drawn_nodes_unique` are the two lists `ROI2dataset -nodelist[.nodups]` writes
// (note that `-nodelist` does NOT drop junction repeats).
// ---------------------------------------------------------------------------

//! The ROI model: [`Roi`], its strokes, and canonical node sets.

use crate::color::Rgba;
use crate::column::{ColumnData, ColumnRole, DataColumn};
use crate::dataset::{Dataset, DatasetKind, ParentIds};
use crate::domain::{DomainId, SurfaceDomain};
use crate::error::{Error, Result};
use crate::labels::LabelEntry;
use crate::topology::SurfaceTopology;

// ---------------------------------------------------------------------------
// Small enums (each keeps unknown codes)
// ---------------------------------------------------------------------------

/// How an ROI was drawn (SUMA's `SUMA_ROI_DRAWING_TYPE`, the file's `Type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RoiDrawingType {
    /// 0: an open path of connected nodes.
    OpenPath,
    /// 1: a closed path.
    ClosedPath,
    /// 2: a filled closed path.
    FilledArea,
    /// 3: a plain collection of nodes.
    Collection,
    /// A code SUMA does not define; kept as written.
    Other(i32),
}

impl RoiDrawingType {
    /// The variant for a file's `Type` code (unknown codes are kept).
    pub fn from_code(code: i32) -> Self {
        match code {
            0 => Self::OpenPath,
            1 => Self::ClosedPath,
            2 => Self::FilledArea,
            3 => Self::Collection,
            other => Self::Other(other),
        }
    }

    /// The code to write.
    pub fn code(self) -> i32 {
        match self {
            Self::OpenPath => 0,
            Self::ClosedPath => 1,
            Self::FilledArea => 2,
            Self::Collection => 3,
            Self::Other(c) => c,
        }
    }
}

/// What a stroke's list describes (SUMA's `SUMA_ROI_TYPE`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RoiElementKind {
    /// 1: a set of nodes.
    NodeGroup,
    /// 2: a set of edges.
    EdgeGroup,
    /// 3: a set of faces.
    FaceGroup,
    /// 4: a series of connected nodes.
    NodeSegment,
    /// 0 (undefined) or a code SUMA does not define; kept as written.
    Other(i32),
}

impl RoiElementKind {
    /// The variant for an element-type code.
    pub fn from_code(code: i32) -> Self {
        match code {
            1 => Self::NodeGroup,
            2 => Self::EdgeGroup,
            3 => Self::FaceGroup,
            4 => Self::NodeSegment,
            other => Self::Other(other),
        }
    }

    /// The code to write.
    pub fn code(self) -> i32 {
        match self {
            Self::NodeGroup => 1,
            Self::EdgeGroup => 2,
            Self::FaceGroup => 3,
            Self::NodeSegment => 4,
            Self::Other(c) => c,
        }
    }
}

/// The drawing action that produced a stroke (SUMA's `SUMA_BRUSH_STROKE_ACTION`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RoiBrushAction {
    /// 1: add the stroke.
    AppendStroke,
    /// 2: add the stroke, or fill if it closes an area.
    AppendStrokeOrFill,
    /// 3: join the path's ends.
    JoinEnds,
    /// 4: fill the enclosed area.
    FillArea,
    /// 0 (undefined) or a code SUMA does not define; kept as written.
    Other(i32),
}

impl RoiBrushAction {
    /// The variant for an action code.
    pub fn from_code(code: i32) -> Self {
        match code {
            1 => Self::AppendStroke,
            2 => Self::AppendStrokeOrFill,
            3 => Self::JoinEnds,
            4 => Self::FillArea,
            other => Self::Other(other),
        }
    }

    /// The code to write.
    pub fn code(self) -> i32 {
        match self {
            Self::AppendStroke => 1,
            Self::AppendStrokeOrFill => 2,
            Self::JoinEnds => 3,
            Self::FillArea => 4,
            Self::Other(c) => c,
        }
    }
}

/// A surface's hemisphere, as SUMA writes it (`Parent_side`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RoiSide {
    /// `L`.
    Left,
    /// `R`.
    Right,
    /// `LR`: both hemispheres.
    Both,
    /// `no_side`.
    NoSide,
    /// Any other text; kept as written.
    Other(String),
}

impl RoiSide {
    /// The variant for a `Parent_side` string.
    pub fn from_name(name: &str) -> Self {
        match name {
            "L" => Self::Left,
            "R" => Self::Right,
            "LR" => Self::Both,
            "no_side" => Self::NoSide,
            other => Self::Other(other.to_owned()),
        }
    }

    /// The text to write.
    pub fn name(&self) -> &str {
        match self {
            Self::Left => "L",
            Self::Right => "R",
            Self::Both => "LR",
            Self::NoSide => "no_side",
            Self::Other(s) => s,
        }
    }
}

/// Where an ROI is in its life.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RoiDrawStatus {
    /// Being drawn.
    InCreation,
    /// Done (every ROI read from a file is finished).
    Finished,
    /// Being edited again.
    InEdit,
}

/// How an ROI came to exist (viewer-side provenance; not in a `.niml.roi` file).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RoiSource {
    /// Typed in or built by code.
    Manual,
    /// Drawn on a surface.
    Drawn,
    /// Read from a NIML ROI file.
    NimlRoi,
    /// Taken from a dataset (for example a label dataset).
    Dataset,
    /// Made from a thresholded overlay (clusters, for example).
    ThresholdedOverlay,
    /// Imported from another format.
    Imported,
    /// Anything else; the text is kept.
    Other(String),
}

/// An ROI's own identity (the file's `self_idcode`).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RoiId(String);

impl RoiId {
    /// A non-blank id (surrounding whitespace is trimmed).
    pub fn new(id: impl Into<String>) -> Result<Self> {
        let id = id.into();
        let id = id.trim();
        if id.is_empty() {
            return Err(Error::Empty("ROI id".into()));
        }
        Ok(Self(id.to_owned()))
    }

    /// The text of the id.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

// ---------------------------------------------------------------------------
// NodeSet
// ---------------------------------------------------------------------------

/// A set of node indices in canonical form: sorted ascending, no repeats. Two sets
/// with the same members are always equal, whatever order they were built in.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct NodeSet(Vec<u32>);

impl NodeSet {
    /// The empty set.
    pub fn empty() -> Self {
        Self(Vec::new())
    }

    /// The set of the given nodes, in any order and with any repeats.
    pub fn new(nodes: impl IntoIterator<Item = u32>) -> Self {
        let mut v: Vec<u32> = nodes.into_iter().collect();
        v.sort_unstable();
        v.dedup();
        Self(v)
    }

    /// The nodes whose entry in `mask` is `true`.
    pub fn from_mask(mask: &[bool]) -> Self {
        // Already ascending and unique, so no sort is needed.
        Self(
            mask.iter()
                .enumerate()
                .filter(|(_, &m)| m)
                .map(|(i, _)| i as u32)
                .collect(),
        )
    }

    /// The members as a sorted slice.
    pub fn as_slice(&self) -> &[u32] {
        &self.0
    }

    /// Number of members.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether the set has no members.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Whether `node` is a member (binary search).
    pub fn contains(&self, node: u32) -> bool {
        self.0.binary_search(&node).is_ok()
    }

    /// Smallest and largest member, if any.
    pub fn range(&self) -> Option<(u32, u32)> {
        Some((*self.0.first()?, *self.0.last()?))
    }

    /// Iterate the members in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.0.iter().copied()
    }

    /// Members of either set.
    pub fn union(&self, other: &Self) -> Self {
        let (mut i, mut j) = (0, 0);
        let mut out = Vec::with_capacity(self.len() + other.len());
        // Merge two sorted lists, keeping one copy of a node present in both.
        while i < self.0.len() && j < other.0.len() {
            match self.0[i].cmp(&other.0[j]) {
                std::cmp::Ordering::Less => {
                    out.push(self.0[i]);
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    out.push(other.0[j]);
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    out.push(self.0[i]);
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&self.0[i..]);
        out.extend_from_slice(&other.0[j..]);
        Self(out)
    }

    /// Members of both sets.
    pub fn intersection(&self, other: &Self) -> Self {
        Self(
            self.0
                .iter()
                .copied()
                .filter(|&n| other.contains(n))
                .collect(),
        )
    }

    /// Members of `self` that are not in `other`.
    pub fn difference(&self, other: &Self) -> Self {
        Self(
            self.0
                .iter()
                .copied()
                .filter(|&n| !other.contains(n))
                .collect(),
        )
    }

    /// Error unless every member is a valid node of a surface with `node_count` nodes.
    pub fn validate_within(&self, node_count: usize) -> Result<()> {
        match self.0.last() {
            Some(&max) if max as usize >= node_count => Err(Error::IndexOutOfRange {
                index: i64::from(max),
                len: node_count,
            }),
            _ => Ok(()),
        }
    }

    /// One flag per node: `true` for members. Errors if a member is out of range.
    pub fn to_mask(&self, node_count: usize) -> Result<Vec<bool>> {
        self.validate_within(node_count)?;
        let mut mask = vec![false; node_count];
        for &n in &self.0 {
            mask[n as usize] = true;
        }
        Ok(mask)
    }
}

impl FromIterator<u32> for NodeSet {
    fn from_iter<T: IntoIterator<Item = u32>>(iter: T) -> Self {
        Self::new(iter)
    }
}

// ---------------------------------------------------------------------------
// Strokes and ROIs
// ---------------------------------------------------------------------------

/// One stroke of an ROI (one record of a `.niml.roi` body): an ordered list of
/// nodes produced by one drawing action.
#[derive(Debug, Clone, PartialEq)]
pub struct RoiStroke {
    /// What the list describes.
    pub kind: RoiElementKind,
    /// The action that made it.
    pub action: RoiBrushAction,
    /// The nodes, in drawing order.
    pub nodes: Vec<u32>,
    /// Triangles, for a face-group stroke (not stored in NIML files).
    pub triangles: Vec<u32>,
    /// Straight-line length of the stroke, if known (not stored in NIML files).
    pub node_distance: Option<f32>,
    /// Length along the surface, if known (not stored in NIML files).
    pub surface_distance: Option<f32>,
}

impl RoiStroke {
    /// A stroke of the given kind and action over `nodes`.
    pub fn new(kind: RoiElementKind, action: RoiBrushAction, nodes: Vec<u32>) -> Self {
        Self {
            kind,
            action,
            nodes,
            triangles: Vec::new(),
            node_distance: None,
            surface_distance: None,
        }
    }

    /// An unordered group of nodes with no drawing action.
    pub fn node_group(nodes: Vec<u32>) -> Self {
        Self::new(RoiElementKind::NodeGroup, RoiBrushAction::Other(0), nodes)
    }

    /// A drawn path (a series of connected nodes).
    pub fn node_segment(nodes: Vec<u32>, action: RoiBrushAction) -> Self {
        Self::new(RoiElementKind::NodeSegment, action, nodes)
    }

    /// Check the stroke has something in it and its distances are sane.
    pub fn validate(&self) -> Result<()> {
        let empty = match self.kind {
            RoiElementKind::FaceGroup => self.triangles.is_empty(),
            RoiElementKind::Other(_) => self.nodes.is_empty() && self.triangles.is_empty(),
            _ => self.nodes.is_empty(),
        };
        if empty {
            return Err(Error::Empty("ROI stroke".into()));
        }
        for (what, d) in [
            ("node distance", self.node_distance),
            ("surface distance", self.surface_distance),
        ] {
            if let Some(d) = d {
                crate::numeric::ensure_finite(what, f64::from(d))?;
                if d < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: what.into(),
                        reason: format!("{d} is negative"),
                    });
                }
            }
        }
        Ok(())
    }
}

/// A region of interest on a surface.
#[derive(Debug, Clone, PartialEq)]
pub struct Roi {
    /// The ROI's own id (the file's `self_idcode`), if it has one.
    pub id: Option<RoiId>,
    /// The surface domain it was drawn on (the file's `domain_parent_idcode`).
    pub parent_domain: Option<DomainId>,
    /// The id of the surface it was drawn on, if that differs from the domain.
    pub parent_surface: Option<String>,
    /// The hemisphere.
    pub parent_side: Option<RoiSide>,
    /// Human-readable name.
    pub label: String,
    /// The integer written in a dataset made from this ROI.
    pub integer_label: i32,
    /// Fill color.
    pub fill_color: Rgba,
    /// Outline color.
    pub edge_color: Rgba,
    /// Outline thickness in pixels.
    pub edge_thickness: u32,
    /// Name of the color plane the viewer draws it in (the file's `ColPlaneName`).
    pub color_plane: Option<String>,
    /// Whether to color by the label's table color instead of `fill_color`
    /// (viewer-side; not in a NIML file).
    pub color_by_label: bool,
    /// How it was drawn.
    pub drawing_type: RoiDrawingType,
    /// Where it is in its life.
    pub draw_status: RoiDrawStatus,
    /// Where it came from.
    pub source: RoiSource,
    /// The source's own identifier, if any.
    pub source_id: Option<String>,
    /// The ordered strokes.
    pub strokes: Vec<RoiStroke>,
}

impl Roi {
    /// An empty, finished ROI with SUMA's default look (translucent red fill, black
    /// 1-pixel outline). `label` must not be blank.
    pub fn new(label: impl Into<String>, integer_label: i32) -> Result<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::Empty("ROI label".into()));
        }
        Ok(Self {
            id: None,
            parent_domain: None,
            parent_surface: None,
            parent_side: None,
            label,
            integer_label,
            fill_color: Rgba::from_u8(255, 0, 0, 180),
            edge_color: Rgba::from_u8(0, 0, 0, 255),
            edge_thickness: 1,
            color_plane: None,
            color_by_label: false,
            drawing_type: RoiDrawingType::Collection,
            draw_status: RoiDrawStatus::Finished,
            source: RoiSource::Manual,
            source_id: None,
            strokes: Vec::new(),
        })
    }

    /// A collection ROI made of the given nodes (one node-group stroke).
    pub fn from_nodes(
        label: impl Into<String>,
        integer_label: i32,
        nodes: impl IntoIterator<Item = u32>,
    ) -> Result<Self> {
        let mut roi = Self::new(label, integer_label)?;
        // Keep the nodes in the order given: this is a group, but order is free to keep.
        let stroke = RoiStroke::node_group(nodes.into_iter().collect());
        stroke.validate()?;
        roi.strokes.push(stroke);
        Ok(roi)
    }

    /// The nodes in drawing order. A node equal to the last node of the previous
    /// stroke is dropped (the junction of two strokes), as SUMA does; repeats
    /// elsewhere are kept.
    pub fn ordered_nodes(&self) -> Vec<u32> {
        let mut out = Vec::new();
        let mut last_of_previous: Option<u32> = None;
        for stroke in &self.strokes {
            for &n in &stroke.nodes {
                if Some(n) != last_of_previous {
                    out.push(n);
                }
            }
            last_of_previous = stroke.nodes.last().copied();
        }
        out
    }

    /// Every node of every stroke exactly as drawn, repeats and junction nodes
    /// included (what `ROI2dataset -nodelist` writes).
    pub fn drawn_nodes(&self) -> Vec<u32> {
        self.strokes
            .iter()
            .flat_map(|s| s.nodes.iter().copied())
            .collect()
    }

    /// The drawn order with each node kept only the first time it appears (what
    /// `ROI2dataset -nodelist.nodups` writes): the path a stroke traced, without
    /// revisits.
    pub fn drawn_nodes_unique(&self) -> Vec<u32> {
        let mut seen = std::collections::BTreeSet::new();
        self.drawn_nodes()
            .into_iter()
            .filter(|&n| seen.insert(n))
            .collect()
    }

    /// The nodes as a canonical set.
    pub fn node_set(&self) -> NodeSet {
        self.strokes
            .iter()
            .flat_map(|s| s.nodes.iter().copied())
            .collect()
    }

    /// Whether `node` is in any stroke.
    pub fn contains_node(&self, node: u32) -> bool {
        self.strokes.iter().any(|s| s.nodes.contains(&node))
    }

    /// Smallest and largest node, if there are any.
    pub fn node_range(&self) -> Option<(u32, u32)> {
        self.node_set().range()
    }

    /// Check the ROI's own consistency: a label, non-empty strokes, a color with
    /// finite channels.
    pub fn validate(&self) -> Result<()> {
        if self.label.trim().is_empty() {
            return Err(Error::Empty("ROI label".into()));
        }
        for stroke in &self.strokes {
            stroke.validate()?;
        }
        Ok(())
    }

    /// Check the ROI against a surface domain: the parent domain (when both sides
    /// know one) must match, and every node must exist. Triangle numbers cannot be
    /// checked without the faces; see [`validate_for_topology`](Self::validate_for_topology).
    pub fn validate_for_domain(&self, domain: &SurfaceDomain) -> Result<()> {
        self.validate()?;
        if let (Some(mine), Some(theirs)) = (&self.parent_domain, domain.id()) {
            if mine != theirs {
                return Err(Error::InvalidParameter {
                    name: "ROI parent domain".into(),
                    reason: format!(
                        "the ROI belongs to domain '{}', not '{}'",
                        mine.as_str(),
                        theirs.as_str()
                    ),
                });
            }
        }
        self.node_set().validate_within(domain.node_count())
    }

    /// Like [`validate_for_domain`](Self::validate_for_domain) but also checks the
    /// triangle numbers against a surface's faces.
    pub fn validate_for_topology(&self, topology: &SurfaceTopology) -> Result<()> {
        self.validate()?;
        self.node_set().validate_within(topology.node_count())?;
        for stroke in &self.strokes {
            for &t in &stroke.triangles {
                if t as usize >= topology.face_count() {
                    return Err(Error::IndexOutOfRange {
                        index: i64::from(t),
                        len: topology.face_count(),
                    });
                }
            }
        }
        Ok(())
    }

    /// The ROI as a label-table entry (key = integer label, name = label, color =
    /// fill color).
    pub fn label_entry(&self) -> LabelEntry {
        LabelEntry {
            key: i64::from(self.integer_label),
            name: self.label.clone(),
            rgba: Some([
                self.fill_color.r,
                self.fill_color.g,
                self.fill_color.b,
                self.fill_color.a,
            ]),
        }
    }

    /// The ROI as a sparse dataset on `domain`: one row per node, every row holding
    /// the integer label. Errors if the ROI has no nodes or does not fit the domain.
    pub fn to_dataset(&self, domain: &SurfaceDomain) -> Result<Dataset> {
        rois_to_dataset(
            std::slice::from_ref(self),
            domain,
            &RoiDatasetOptions::default(),
        )
    }
}

// ---------------------------------------------------------------------------
// ROIs to a dataset (SUMA's ROI2dataset)
// ---------------------------------------------------------------------------

/// Who keeps a node that two ROIs with DIFFERENT labels both hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OverlapPolicy {
    /// The first ROI (in the order given) keeps it. Deterministic, and the usual
    /// reading of SUMA's "duplicate entries were eliminated". SUMA itself does not
    /// promise this: it sorts with the C library's `qsort`, so which label survives
    /// depends on the platform (see the roadmap log).
    #[default]
    FirstWins,
    /// The last ROI keeps it.
    LastWins,
    /// Refuse: return an error naming a contested node.
    Error,
}

/// Options for [`rois_to_dataset`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RoiDatasetOptions {
    /// What to do when ROIs with different labels share a node.
    pub overlap: OverlapPolicy,
    /// Write a row for every node from 0 to this index (inclusive), giving nodes
    /// that are in no ROI `pad_label`. `None` writes only the ROI nodes.
    /// (`ROI2dataset -pad_to_node`.)
    pub pad_to: Option<u32>,
    /// The label of a padded node (`-pad_label`; SUMA's default is 0).
    pub pad_label: i32,
}

/// Combine ROIs into one dataset on `domain`: a sparse dataset of kind `Roi` with a
/// node-index mapping and one `Int32` label column, as SUMA's `ROI2dataset` writes.
///
/// A node held by several ROIs with the same label is simply listed once. If the
/// labels differ, `options.overlap` decides (by default the first ROI keeps it; SUMA
/// leaves this to the C library's sort). Use [`overlapping_nodes`] to find out
/// whether ROIs overlap. Errors if no ROI has a
/// node, if a node is outside the domain, or if `pad_to` is smaller than the
/// largest node (SUMA refuses that case too).
pub fn rois_to_dataset(
    rois: &[Roi],
    domain: &SurfaceDomain,
    options: &RoiDatasetOptions,
) -> Result<Dataset> {
    for roi in rois {
        roi.validate_for_domain(domain)?;
    }
    // node -> label, resolving contested nodes by the overlap policy.
    let mut claimed: std::collections::BTreeMap<u32, i32> = std::collections::BTreeMap::new();
    for roi in rois {
        for node in roi.node_set().iter() {
            match claimed.get(&node).copied() {
                None => {
                    claimed.insert(node, roi.integer_label);
                }
                Some(label) if label == roi.integer_label => {}
                Some(label) => match options.overlap {
                    OverlapPolicy::FirstWins => {}
                    OverlapPolicy::LastWins => {
                        claimed.insert(node, roi.integer_label);
                    }
                    OverlapPolicy::Error => {
                        return Err(Error::InvalidParameter {
                            name: "ROI overlap".into(),
                            reason: format!(
                                "node {node} is in ROIs labelled {label} and {}",
                                roi.integer_label
                            ),
                        })
                    }
                },
            }
        }
    }
    let max_node = claimed.keys().next_back().copied();
    let max_node = max_node.ok_or_else(|| Error::Empty("ROI nodes".into()))?;

    let (nodes, labels): (Vec<u32>, Vec<i32>) = match options.pad_to {
        Some(pad_to) => {
            if pad_to < max_node {
                return Err(Error::InvalidParameter {
                    name: "pad_to".into(),
                    reason: format!(
                        "an ROI contains node {max_node}, beyond the padding limit {pad_to}"
                    ),
                });
            }
            // Every node from 0 to pad_to; ROI nodes keep their label.
            (0..=pad_to)
                .map(|n| (n, claimed.get(&n).copied().unwrap_or(options.pad_label)))
                .unzip()
        }
        None => claimed.into_iter().unzip(),
    };
    // Padding past the domain would make an invalid row map; say so plainly.
    if let Some(&last) = nodes.last() {
        if last as usize >= domain.node_count() {
            return Err(Error::IndexOutOfRange {
                index: i64::from(last),
                len: domain.node_count(),
            });
        }
    }
    let column = DataColumn::new(
        "integer label",
        ColumnRole::Label,
        ColumnData::Int32(labels),
    )?;
    let dataset = Dataset::indexed(
        DatasetKind::Roi,
        crate::domain::Domain::Surface(domain.clone()),
        nodes,
        vec![column],
    )?;
    let ids = ParentIds {
        self_id: None,
        domain_parent: domain.id().map(|d| d.as_str().to_owned()),
        geometry_parent: rois.iter().find_map(|r| r.parent_surface.clone()),
    };
    Ok(dataset.with_parent_ids(ids))
}

/// The nodes that belong to more than one of `rois` (what SUMA warns about when it
/// "eliminates duplicate entries").
pub fn overlapping_nodes(rois: &[Roi]) -> NodeSet {
    let mut seen: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
    let mut dup = Vec::new();
    for roi in rois {
        for node in roi.node_set().iter() {
            if !seen.insert(node) {
                dup.push(node);
            }
        }
    }
    NodeSet::new(dup)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn domain(n: usize) -> SurfaceDomain {
        SurfaceDomain::new(Some(DomainId::new("surf-1").unwrap()), n).unwrap()
    }

    #[test]
    fn codes_round_trip_including_unknown_ones() {
        for code in -1..8 {
            assert_eq!(RoiDrawingType::from_code(code).code(), code);
            assert_eq!(RoiElementKind::from_code(code).code(), code);
            assert_eq!(RoiBrushAction::from_code(code).code(), code);
        }
        // SUMA defines 0..=3; sumaru's 4 is kept, not lost.
        assert_eq!(RoiDrawingType::from_code(4), RoiDrawingType::Other(4));
        assert_eq!(RoiDrawingType::from_code(2), RoiDrawingType::FilledArea);
        assert_eq!(RoiElementKind::from_code(0), RoiElementKind::Other(0));
        for name in ["L", "R", "LR", "no_side", "weird"] {
            assert_eq!(RoiSide::from_name(name).name(), name);
        }
        assert_eq!(RoiSide::from_name("weird"), RoiSide::Other("weird".into()));
    }

    #[test]
    fn node_sets_are_canonical_and_support_set_algebra() {
        let a = NodeSet::new([5, 1, 3, 3, 1]);
        assert_eq!(a.as_slice(), &[1, 3, 5]);
        let b = NodeSet::new([3, 4, 5, 6]);
        assert_eq!(a.union(&b).as_slice(), &[1, 3, 4, 5, 6]);
        assert_eq!(a.intersection(&b).as_slice(), &[3, 5]);
        assert_eq!(a.difference(&b).as_slice(), &[1]);
        assert_eq!(b.difference(&a).as_slice(), &[4, 6]);
        assert!(a.contains(3) && !a.contains(2));
        assert_eq!(a.range(), Some((1, 5)));
        assert_eq!(NodeSet::empty().range(), None);
        assert!(a.validate_within(6).is_ok());
        assert!(a.validate_within(5).is_err());
        assert_eq!(
            a.to_mask(6).unwrap(),
            [false, true, false, true, false, true]
        );
        assert_eq!(NodeSet::from_mask(&a.to_mask(6).unwrap()), a);
    }

    #[test]
    fn ordered_nodes_drop_only_stroke_junction_repeats() {
        let mut roi = Roi::new("path", 1).unwrap();
        roi.strokes.push(RoiStroke::node_segment(
            vec![4, 5, 6],
            RoiBrushAction::AppendStroke,
        ));
        // Starts at the previous stroke's last node (6): that one is dropped.
        roi.strokes.push(RoiStroke::node_segment(
            vec![6, 7, 5],
            RoiBrushAction::AppendStroke,
        ));
        // 5 appears again (not a junction), so it stays.
        assert_eq!(roi.ordered_nodes(), vec![4, 5, 6, 7, 5]);
        assert_eq!(roi.node_set().as_slice(), &[4, 5, 6, 7]);
        // The raw lists: every node as drawn, then first occurrences only.
        assert_eq!(roi.drawn_nodes(), vec![4, 5, 6, 6, 7, 5]);
        assert_eq!(roi.drawn_nodes_unique(), vec![4, 5, 6, 7]);
        assert_eq!(roi.node_range(), Some((4, 7)));
        assert!(roi.contains_node(7) && !roi.contains_node(8));
    }

    #[test]
    fn validation_against_a_domain() {
        let roi = Roi::from_nodes("a", 1, [0, 3, 9]).unwrap();
        assert!(roi.validate_for_domain(&domain(10)).is_ok());
        assert!(roi.validate_for_domain(&domain(9)).is_err());
        // A different parent domain is refused; a matching one is fine.
        let mut other = roi.clone();
        other.parent_domain = Some(DomainId::new("surf-2").unwrap());
        assert!(other.validate_for_domain(&domain(10)).is_err());
        other.parent_domain = Some(DomainId::new("surf-1").unwrap());
        assert!(other.validate_for_domain(&domain(10)).is_ok());
        assert!(Roi::new("  ", 1).is_err());
        assert!(RoiStroke::node_group(vec![]).validate().is_err());
    }

    #[test]
    fn one_roi_becomes_a_sparse_label_dataset() {
        let roi = Roi::from_nodes("blob", 7, [9, 2, 2, 5]).unwrap();
        let d = roi.to_dataset(&domain(10)).unwrap();
        assert_eq!(d.kind(), &DatasetKind::Roi);
        assert!(d.is_sparse());
        assert_eq!(d.row_count(), 3);
        assert_eq!(d.sample_for_row(0), Some(2));
        assert_eq!(d.parent_ids().domain_parent.as_deref(), Some("surf-1"));
        let entry = roi.label_entry();
        assert_eq!((entry.key, entry.name.as_str()), (7, "blob"));
    }

    #[test]
    fn several_rois_first_one_wins_and_padding_works() {
        let a = Roi::from_nodes("a", 1, [1, 2, 3]).unwrap();
        let b = Roi::from_nodes("b", 2, [3, 4]).unwrap();
        assert_eq!(overlapping_nodes(&[a.clone(), b.clone()]).as_slice(), &[3]);
        let dom = domain(8);
        let d =
            rois_to_dataset(&[a.clone(), b.clone()], &dom, &RoiDatasetOptions::default()).unwrap();
        let labels = match d.columns()[0].values() {
            ColumnData::Int32(v) => v.clone(),
            _ => panic!("labels are Int32"),
        };
        // Nodes 1,2,3,4 -> 1,1,1(first ROI wins),2.
        assert_eq!(labels, vec![1, 1, 1, 2]);
        // Padding to node 6 gives 7 rows; unclaimed nodes get the pad label.
        let padded = rois_to_dataset(
            &[a.clone(), b.clone()],
            &dom,
            &RoiDatasetOptions {
                pad_to: Some(6),
                pad_label: -1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(padded.row_count(), 7);
        let labels = match padded.columns()[0].values() {
            ColumnData::Int32(v) => v.clone(),
            _ => panic!(),
        };
        assert_eq!(labels, vec![-1, 1, 1, 1, 2, -1, -1]);
        // Padding below the largest node, or past the domain, is an error.
        let c = Roi::from_nodes("c", 3, [5]).unwrap();
        assert!(rois_to_dataset(
            std::slice::from_ref(&c),
            &dom,
            &RoiDatasetOptions {
                pad_to: Some(4),
                ..Default::default()
            }
        )
        .is_err());
        assert!(rois_to_dataset(
            &[c],
            &dom,
            &RoiDatasetOptions {
                pad_to: Some(20),
                ..Default::default()
            }
        )
        .is_err());
        // The overlap policy: last wins, or refuse.
        let last = rois_to_dataset(
            &[a.clone(), b.clone()],
            &dom,
            &RoiDatasetOptions {
                overlap: OverlapPolicy::LastWins,
                ..Default::default()
            },
        )
        .unwrap();
        let labels = match last.columns()[0].values() {
            ColumnData::Int32(v) => v.clone(),
            _ => panic!(),
        };
        assert_eq!(labels, vec![1, 1, 2, 2]);
        assert!(rois_to_dataset(
            &[a, b],
            &dom,
            &RoiDatasetOptions {
                overlap: OverlapPolicy::Error,
                ..Default::default()
            }
        )
        .is_err());
        // No nodes at all is an error.
        assert!(rois_to_dataset(&[Roi::new("e", 1).unwrap()], &dom, &Default::default()).is_err());
    }
}
