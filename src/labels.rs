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
// A minimal semantic label table: integer keys with names and optional
// colors, in file order. Used as column metadata for label datasets
// (for example a FreeSurfer parcellation, where key 1007 means "fusiform").
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `column.rs` attaches an optional `LabelTable` to a column.
// * Phase 4 (colors) will extend this with fallback colors, key-zero policy and
//   lookup helpers. Phase 1 only needs identity: keys stay integers (never
//   squeezed through `f32`) and order is preserved.
// * `afni-io` keeps the raw table (with every attribute, for round trips); an
//   adapter there converts to and from this type.
// ---------------------------------------------------------------------------

//! Label tables: integer keys, names, optional colors.

use crate::error::{Error, Result};

/// One label.
#[derive(Debug, Clone, PartialEq)]
pub struct LabelEntry {
    /// The integer stored in the data. Kept as `i64`, never as a float.
    pub key: i64,
    /// Human-readable name.
    pub name: String,
    /// Display color as RGBA in `[0, 1]`, if the source has one.
    pub rgba: Option<[f32; 4]>,
}

/// A label table with unique keys, in the order the source listed them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LabelTable {
    entries: Vec<LabelEntry>,
}

impl LabelTable {
    /// Build a table. Keys must be unique (a repeated key would make "what is
    /// label 7?" ambiguous) and any color component must be finite.
    pub fn new(entries: Vec<LabelEntry>) -> Result<Self> {
        for (i, e) in entries.iter().enumerate() {
            if entries[..i].iter().any(|prev| prev.key == e.key) {
                return Err(Error::InvalidParameter {
                    name: "label key".into(),
                    reason: format!("key {} appears more than once", e.key),
                });
            }
            if let Some(rgba) = e.rgba {
                for c in rgba {
                    crate::numeric::ensure_finite("label color", f64::from(c))?;
                }
            }
        }
        Ok(Self { entries })
    }

    /// The entries, in source order.
    pub fn entries(&self) -> &[LabelEntry] {
        &self.entries
    }

    /// Number of labels.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table has no labels.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The entry for `key`, if any.
    pub fn get(&self, key: i64) -> Option<&LabelEntry> {
        self.entries.iter().find(|e| e.key == key)
    }

    /// The key whose name is exactly `name`, if any.
    pub fn key_for(&self, name: &str) -> Option<i64> {
        self.entries.iter().find(|e| e.name == name).map(|e| e.key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(key: i64, name: &str) -> LabelEntry {
        LabelEntry {
            key,
            name: name.into(),
            rgba: None,
        }
    }

    #[test]
    fn lookup_both_directions_and_order_kept() {
        let t = LabelTable::new(vec![entry(1007, "fusiform"), entry(3, "other")]).unwrap();
        assert_eq!(t.get(1007).unwrap().name, "fusiform");
        assert_eq!(t.key_for("other"), Some(3));
        assert_eq!(t.get(5), None);
        assert_eq!(t.entries()[0].key, 1007);
    }

    #[test]
    fn duplicate_keys_and_bad_colors_are_rejected() {
        assert!(LabelTable::new(vec![entry(1, "a"), entry(1, "b")]).is_err());
        let bad = LabelEntry {
            key: 1,
            name: "a".into(),
            rgba: Some([f32::NAN, 0.0, 0.0, 1.0]),
        };
        assert!(LabelTable::new(vec![bad]).is_err());
    }

    #[test]
    fn large_keys_are_not_squeezed_through_float() {
        // 2^53 + 1 is not representable as f64; i64 keeps it exactly.
        let key = (1_i64 << 53) + 1;
        assert_eq!(
            LabelTable::new(vec![entry(key, "big")])
                .unwrap()
                .get(key)
                .unwrap()
                .key,
            key
        );
    }
}
