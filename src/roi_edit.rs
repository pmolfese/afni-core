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
// Editing an ROI as a list of explicit COMMANDS that can be undone and redone:
// add a stroke, remove the last stroke, rename, recolor, change the drawing type
// or status, or a `Batch` of those. `RoiEditor` holds an ROI together with its
// undo and redo history and offers the three drawing actions SUMA has: draw a
// path through picked nodes, join the ends, and fill the enclosed area.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `roi.rs` is the ROI being edited. `roi_ops.rs` supplies the geometry the
//   drawing actions use (shortest paths, filling).
// * The viewer (sumaru, later afniru) owns the MOUSE: it turns clicks into picked
//   node numbers and calls the editor. Keeping picking out of here is what lets a
//   tool, a test or a script make exactly the same edits.
//
// HOW UNDO WORKS
//
// Applying a command returns its INVERSE, another command. The editor keeps pairs of
// (command, inverse): undo applies the inverse, redo applies the command again.
// Because commands are plain data they can also be logged, saved or sent elsewhere.
// A new edit clears the redo list, like every editor.
// ---------------------------------------------------------------------------

//! Undoable ROI edit commands and the [`RoiEditor`].

use crate::color::Rgba;
use crate::error::{Error, Result};
use crate::mesh::SurfaceMesh;
use crate::roi::{
    NodeSet, Roi, RoiBrushAction, RoiDrawStatus, RoiDrawingType, RoiElementKind, RoiStroke,
};
use crate::roi_ops::{fill_enclosed, join_ends, shortest_path};
use crate::topology::SurfaceTopology;

/// One edit of an ROI.
#[derive(Debug, Clone, PartialEq)]
pub enum RoiCommand {
    /// Add a stroke at the end.
    AddStroke(RoiStroke),
    /// Remove the last stroke (an error if there is none).
    RemoveLastStroke,
    /// Change the name and the integer label.
    SetLabel {
        /// New name (not blank).
        label: String,
        /// New integer label.
        integer_label: i32,
    },
    /// Change the fill color, outline color and outline thickness.
    SetColors {
        /// Fill color.
        fill: Rgba,
        /// Outline color.
        edge: Rgba,
        /// Outline thickness in pixels.
        edge_thickness: u32,
    },
    /// Change how the ROI is classified (open path, closed path, filled, ...).
    SetDrawingType(RoiDrawingType),
    /// Change whether it is being drawn, finished or being edited.
    SetDrawStatus(RoiDrawStatus),
    /// Several commands applied in order; undone in reverse order.
    Batch(Vec<RoiCommand>),
}

impl RoiCommand {
    /// Apply the command to `roi` and return the command that undoes it. If a
    /// command inside a `Batch` fails, the ones already applied are undone, so a
    /// failed command leaves the ROI as it was.
    pub fn apply(&self, roi: &mut Roi) -> Result<RoiCommand> {
        match self {
            Self::AddStroke(stroke) => {
                stroke.validate()?;
                roi.strokes.push(stroke.clone());
                Ok(Self::RemoveLastStroke)
            }
            Self::RemoveLastStroke => roi
                .strokes
                .pop()
                .map(Self::AddStroke)
                .ok_or_else(|| Error::Empty("ROI strokes (nothing to remove)".into())),
            Self::SetLabel {
                label,
                integer_label,
            } => {
                if label.trim().is_empty() {
                    return Err(Error::Empty("ROI label".into()));
                }
                let inverse = Self::SetLabel {
                    label: std::mem::replace(&mut roi.label, label.clone()),
                    integer_label: std::mem::replace(&mut roi.integer_label, *integer_label),
                };
                Ok(inverse)
            }
            Self::SetColors {
                fill,
                edge,
                edge_thickness,
            } => Ok(Self::SetColors {
                fill: std::mem::replace(&mut roi.fill_color, *fill),
                edge: std::mem::replace(&mut roi.edge_color, *edge),
                edge_thickness: std::mem::replace(&mut roi.edge_thickness, *edge_thickness),
            }),
            Self::SetDrawingType(t) => Ok(Self::SetDrawingType(std::mem::replace(
                &mut roi.drawing_type,
                *t,
            ))),
            Self::SetDrawStatus(s) => Ok(Self::SetDrawStatus(std::mem::replace(
                &mut roi.draw_status,
                *s,
            ))),
            Self::Batch(commands) => {
                let mut inverses: Vec<RoiCommand> = Vec::with_capacity(commands.len());
                for command in commands {
                    match command.apply(roi) {
                        Ok(inverse) => inverses.push(inverse),
                        Err(e) => {
                            // Roll back what already happened, newest first.
                            for undo in inverses.iter().rev() {
                                // An inverse cannot fail right after its command ran.
                                let _ = undo.apply(roi);
                            }
                            return Err(e);
                        }
                    }
                }
                inverses.reverse();
                Ok(Self::Batch(inverses))
            }
        }
    }

    /// A stroke that appends `nodes` as a drawn path (action "append stroke").
    pub fn append_path(nodes: Vec<u32>) -> Self {
        Self::AddStroke(RoiStroke::node_segment(nodes, RoiBrushAction::AppendStroke))
    }
}

/// An ROI plus its undo and redo history.
#[derive(Debug, Clone, PartialEq)]
pub struct RoiEditor {
    roi: Roi,
    /// (command, its inverse), oldest first.
    done: Vec<(RoiCommand, RoiCommand)>,
    /// Undone pairs, most recently undone last.
    undone: Vec<(RoiCommand, RoiCommand)>,
}

impl RoiEditor {
    /// Start editing `roi` with an empty history.
    pub fn new(roi: Roi) -> Self {
        Self {
            roi,
            done: Vec::new(),
            undone: Vec::new(),
        }
    }

    /// The ROI as it is now.
    pub fn roi(&self) -> &Roi {
        &self.roi
    }

    /// Stop editing and take the ROI.
    pub fn into_roi(self) -> Roi {
        self.roi
    }

    /// Whether there is something to undo.
    pub fn can_undo(&self) -> bool {
        !self.done.is_empty()
    }

    /// Whether there is something to redo.
    pub fn can_redo(&self) -> bool {
        !self.undone.is_empty()
    }

    /// The commands applied so far (oldest first), e.g. to save a log.
    pub fn history(&self) -> Vec<&RoiCommand> {
        self.done.iter().map(|(c, _)| c).collect()
    }

    /// Apply a command and record it. Clears the redo list. On error the ROI is
    /// unchanged and nothing is recorded.
    pub fn apply(&mut self, command: RoiCommand) -> Result<()> {
        let inverse = command.apply(&mut self.roi)?;
        self.done.push((command, inverse));
        self.undone.clear();
        Ok(())
    }

    /// Undo the last command. Returns `false` if there was nothing to undo.
    pub fn undo(&mut self) -> bool {
        match self.done.pop() {
            Some((command, inverse)) => {
                // An inverse cannot fail on the state its command produced.
                let _ = inverse.apply(&mut self.roi);
                self.undone.push((command, inverse));
                true
            }
            None => false,
        }
    }

    /// Redo the last undone command. Returns `false` if there was nothing to redo.
    pub fn redo(&mut self) -> bool {
        match self.undone.pop() {
            Some((command, inverse)) => {
                let _ = command.apply(&mut self.roi);
                self.done.push((command, inverse));
                true
            }
            None => false,
        }
    }

    // ----- the three drawing actions -----

    /// Draw a path through the picked nodes. Consecutive picks are joined by the
    /// shortest path along the surface (what SUMA does as you click), and if the
    /// ROI already has strokes the new path starts from its last node. The first
    /// path drawn on an empty ROI makes it an open path.
    pub fn draw_path(&mut self, mesh: &SurfaceMesh, picks: &[u32]) -> Result<()> {
        if picks.is_empty() {
            return Err(Error::Empty("picked nodes".into()));
        }
        let mut path: Vec<u32> = Vec::new();
        // The path continues from the last node of the ROI, if it has one.
        let mut previous = self.roi.ordered_nodes().last().copied();
        for &pick in picks {
            match previous {
                Some(prev) if prev != pick => {
                    let leg = shortest_path(mesh, prev, pick)?;
                    // The leg starts at `prev`, which is already in the path.
                    let skip = usize::from(!path.is_empty() || self.roi.strokes_have_nodes());
                    path.extend_from_slice(&leg[skip.min(leg.len())..]);
                }
                Some(_) => {}
                None => path.push(pick),
            }
            previous = Some(pick);
        }
        if path.is_empty() {
            path.push(picks[0]);
        }
        let mut commands = vec![RoiCommand::append_path(path)];
        if self.roi.strokes.is_empty() {
            commands.push(RoiCommand::SetDrawingType(RoiDrawingType::OpenPath));
        }
        self.apply(RoiCommand::Batch(commands))
    }

    /// Close the path: add the shortest path from its last node back to its first,
    /// and mark the ROI a closed path. Errors for an ROI with no nodes.
    pub fn join_ends(&mut self, mesh: &SurfaceMesh) -> Result<()> {
        let nodes = self.roi.ordered_nodes();
        let segment = join_ends(mesh, &nodes)?;
        self.apply(RoiCommand::Batch(vec![
            RoiCommand::AddStroke(RoiStroke::node_segment(segment, RoiBrushAction::JoinEnds)),
            RoiCommand::SetDrawingType(RoiDrawingType::ClosedPath),
        ]))
    }

    /// Fill the area around `seed` enclosed by the ROI's current nodes, and mark the
    /// ROI a filled area. The new stroke holds the filled nodes that were not
    /// already in the ROI. Errors if the fill leaks off the rim of an open surface
    /// (the path does not enclose the seed) or if the seed is on the path.
    pub fn fill_area(&mut self, topology: &SurfaceTopology, seed: u32) -> Result<()> {
        let boundary = self.roi.node_set();
        let fill = fill_enclosed(topology, &boundary, seed)?;
        if fill.touches_surface_rim {
            return Err(Error::InvalidParameter {
                name: "fill".into(),
                reason: format!(
                    "the path does not enclose node {seed}: the fill reached the rim of the surface"
                ),
            });
        }
        let added: NodeSet = fill.nodes.difference(&boundary);
        let stroke = RoiStroke::new(
            RoiElementKind::NodeGroup,
            RoiBrushAction::FillArea,
            added.iter().collect(),
        );
        self.apply(RoiCommand::Batch(vec![
            RoiCommand::AddStroke(stroke),
            RoiCommand::SetDrawingType(RoiDrawingType::FilledArea),
        ]))
    }

    /// Mark the ROI finished.
    pub fn finish(&mut self) -> Result<()> {
        self.apply(RoiCommand::SetDrawStatus(RoiDrawStatus::Finished))
    }
}

impl Roi {
    /// Whether any stroke has a node (a helper for the editor).
    fn strokes_have_nodes(&self) -> bool {
        self.strokes.iter().any(|s| !s.nodes.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::roi_ops::check_path;

    fn grid(w: u32, h: u32) -> SurfaceMesh {
        let stride = w + 1;
        let mut vertices = Vec::new();
        for j in 0..=h {
            for i in 0..=w {
                vertices.push([i as f32, j as f32, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for j in 0..h {
            for i in 0..w {
                let a = j * stride + i;
                let (b, c, d) = (a + 1, a + stride, a + stride + 1);
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            }
        }
        SurfaceMesh::from_triangles(vertices, faces).unwrap()
    }

    fn id(w: u32, i: u32, j: u32) -> u32 {
        j * (w + 1) + i
    }

    #[test]
    fn commands_return_their_inverse() {
        let mut roi = Roi::new("a", 1).unwrap();
        let before = roi.clone();
        let add = RoiCommand::append_path(vec![1, 2, 3]);
        let inverse = add.apply(&mut roi).unwrap();
        assert_eq!(roi.strokes.len(), 1);
        inverse.apply(&mut roi).unwrap();
        assert_eq!(roi, before);
        // Removing from an empty ROI is an error.
        assert!(RoiCommand::RemoveLastStroke.apply(&mut roi).is_err());
        // Setters give back the old values.
        let inv = RoiCommand::SetLabel {
            label: "b".into(),
            integer_label: 9,
        }
        .apply(&mut roi)
        .unwrap();
        assert_eq!((roi.label.as_str(), roi.integer_label), ("b", 9));
        inv.apply(&mut roi).unwrap();
        assert_eq!((roi.label.as_str(), roi.integer_label), ("a", 1));
        assert!(RoiCommand::SetLabel {
            label: " ".into(),
            integer_label: 1
        }
        .apply(&mut roi)
        .is_err());
    }

    #[test]
    fn a_failed_batch_changes_nothing() {
        let mut roi = Roi::new("a", 1).unwrap();
        let before = roi.clone();
        let batch = RoiCommand::Batch(vec![
            RoiCommand::append_path(vec![4, 5]),
            RoiCommand::SetDrawingType(RoiDrawingType::OpenPath),
            RoiCommand::AddStroke(RoiStroke::node_group(vec![])), // invalid: empty
        ]);
        assert!(batch.apply(&mut roi).is_err());
        assert_eq!(roi, before);
    }

    #[test]
    fn undo_and_redo_walk_the_history() {
        let m = grid(8, 8);
        let mut ed = RoiEditor::new(Roi::new("path", 3).unwrap());
        assert!(!ed.can_undo() && !ed.undo() && !ed.redo());
        ed.draw_path(&m, &[id(8, 1, 1), id(8, 4, 1)]).unwrap();
        let after_draw = ed.roi().clone();
        assert_eq!(after_draw.drawing_type, RoiDrawingType::OpenPath);
        ed.finish().unwrap();
        ed.apply(RoiCommand::SetDrawingType(RoiDrawingType::Collection))
            .unwrap();
        assert_eq!(ed.history().len(), 3);
        // Undo twice: back to just after the draw.
        assert!(ed.undo() && ed.undo());
        assert_eq!(ed.roi(), &after_draw);
        assert!(ed.can_redo());
        // Redo one, then a NEW edit clears the redo list.
        assert!(ed.redo());
        ed.apply(RoiCommand::SetColors {
            fill: Rgba::from_u8(0, 255, 0, 255),
            edge: Rgba::from_u8(0, 0, 0, 255),
            edge_thickness: 3,
        })
        .unwrap();
        assert!(!ed.can_redo());
        // Undo everything: the empty ROI comes back exactly.
        while ed.undo() {}
        assert_eq!(ed.roi(), &Roi::new("path", 3).unwrap());
    }

    #[test]
    fn drawing_joining_and_filling_a_region() {
        let m = grid(10, 10);
        let t = m.topology();
        let mut ed = RoiEditor::new(Roi::new("loop", 5).unwrap());
        // Click the corners of a square: legs between clicks are shortest paths.
        let corners = [id(10, 3, 3), id(10, 7, 3), id(10, 7, 7), id(10, 3, 7)];
        ed.draw_path(&m, &corners).unwrap();
        let path = ed.roi().ordered_nodes();
        assert!(check_path(t, &path).unwrap().is_connected());
        assert_eq!(path[0], corners[0]);
        assert_eq!(*path.last().unwrap(), corners[3]);
        assert_eq!(ed.roi().drawing_type, RoiDrawingType::OpenPath);
        // Join the ends: the path is now a closed loop.
        ed.join_ends(&m).unwrap();
        assert_eq!(ed.roi().drawing_type, RoiDrawingType::ClosedPath);
        let closed = ed.roi().ordered_nodes();
        assert!(check_path(t, &closed).unwrap().is_closed);
        let loop_nodes = ed.roi().node_set();
        // Filling from the middle gives everything inside; from outside, an error
        // (the fill reaches the rim of the open surface).
        assert!(ed.fill_area(t, id(10, 0, 0)).is_err());
        assert_eq!(
            ed.roi().drawing_type,
            RoiDrawingType::ClosedPath,
            "failed fill changes nothing"
        );
        ed.fill_area(t, id(10, 5, 5)).unwrap();
        assert_eq!(ed.roi().drawing_type, RoiDrawingType::FilledArea);
        let filled = ed.roi().node_set();
        assert!(filled.contains(id(10, 5, 5)) && !filled.contains(id(10, 1, 1)));
        assert!(filled.len() > loop_nodes.len());
        // The last stroke is a fill action holding only the new nodes.
        let last = ed.roi().strokes.last().unwrap();
        assert_eq!(last.action, RoiBrushAction::FillArea);
        assert_eq!(last.nodes.len(), filled.len() - loop_nodes.len());
        // Undo the fill and the join: back to the open path.
        assert!(ed.undo() && ed.undo());
        assert_eq!(ed.roi().drawing_type, RoiDrawingType::OpenPath);
        assert_eq!(ed.roi().ordered_nodes(), path);
    }

    #[test]
    fn drawing_continues_from_the_last_node() {
        let m = grid(10, 10);
        let mut ed = RoiEditor::new(Roi::new("p", 1).unwrap());
        ed.draw_path(&m, &[id(10, 1, 1), id(10, 3, 1)]).unwrap();
        ed.draw_path(&m, &[id(10, 3, 4)]).unwrap();
        let nodes = ed.roi().ordered_nodes();
        // One continuous path: no duplicated junction node, every step an edge.
        assert!(check_path(m.topology(), &nodes).unwrap().is_connected());
        assert!(!check_path(m.topology(), &nodes).unwrap().has_repeats);
        assert_eq!(*nodes.last().unwrap(), id(10, 3, 4));
        assert!(ed.draw_path(&m, &[]).is_err());
        assert!(ed.join_ends(&m).is_ok());
        assert!(RoiEditor::new(Roi::new("x", 1).unwrap())
            .join_ends(&m)
            .is_err());
    }
}
