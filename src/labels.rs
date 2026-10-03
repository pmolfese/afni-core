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
// * Phase 4 adds the color side: a `LabelColorPolicy` (what color an unlabeled
//   key, an uncolored entry, or key 0 gets), SUMA-style stable fallback colors,
//   and `LabelColorMap`, a table plus a policy. Keys stay integers (never
//   squeezed through `f32`), order is preserved, and a table's own colors are
//   returned exactly as stored (no quantization to bytes).
// * `afni-io` keeps the raw table (with every attribute, for round trips); an
//   adapter there converts to and from this type.
// ---------------------------------------------------------------------------

//! Label tables: integer keys, names, optional colors, and label-to-color lookup.

use crate::color::Rgba;
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

// ---------------------------------------------------------------------------
// Colors for labels
// ---------------------------------------------------------------------------

/// Ten saturated, well-separated colors for labels that have none of their own.
/// This palette is a sumaru design choice, not an AFNI definition; it exists so a
/// label dataset with no color table still draws distinguishable regions, and so
/// the same key gets the same color in every viewer.
const STABLE_PALETTE: [[u8; 3]; 10] = [
    [0, 194, 255],
    [255, 242, 0],
    [57, 255, 20],
    [255, 117, 24],
    [255, 0, 255],
    [0, 255, 255],
    [255, 49, 49],
    [157, 0, 255],
    [180, 255, 0],
    [255, 0, 170],
];

/// The fallback color for key 0, mid gray.
const ZERO_KEY_GRAY: [u8; 3] = [128, 128, 128];

/// A fallback color that depends only on the key: key 0 is gray; any other key
/// `k` uses palette entry `(|k| - 1) mod 10`, so `k` and `-k` share a color and
/// the assignment never changes between runs or viewers.
///
/// ```
/// use afni_core::labels::stable_label_rgb;
/// assert_eq!(stable_label_rgb(1), [0, 194, 255]);
/// assert_eq!(stable_label_rgb(11), stable_label_rgb(1)); // the palette repeats
/// assert_eq!(stable_label_rgb(0), [128, 128, 128]);
/// ```
pub fn stable_label_rgb(key: i64) -> [u8; 3] {
    if key == 0 {
        return ZERO_KEY_GRAY;
    }
    // `unsigned_abs` is defined even for i64::MIN, unlike `abs`.
    let index = (key.unsigned_abs() - 1) % STABLE_PALETTE.len() as u64;
    STABLE_PALETTE[index as usize]
}

/// Where a label color comes from when the table does not supply one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KeyColor {
    /// One fixed color (use [`Rgba::TRANSPARENT`] for "draw nothing").
    Fixed(Rgba),
    /// The key-based [`stable_label_rgb`] color, opaque.
    StablePalette,
}

impl KeyColor {
    fn resolve(self, key: i64) -> Rgba {
        match self {
            KeyColor::Fixed(c) => c,
            KeyColor::StablePalette => {
                let [r, g, b] = stable_label_rgb(key);
                Rgba::from_u8(r, g, b, 255)
            }
        }
    }
}

/// What color key 0 gets. In many label datasets 0 means "no label".
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ZeroKeyPolicy {
    /// Treat 0 like any other key: use the table's entry if there is one,
    /// otherwise the unlabeled rule (the default; faithful to the file).
    #[default]
    FromTable,
    /// Always use the unlabeled color for 0, even if the table colors it.
    Unlabeled,
    /// Always use this color for 0.
    Fixed(Rgba),
}

/// Rules for choosing a color for a label key.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LabelColorPolicy {
    /// For a key that is not in the table. Default: transparent.
    pub unlabeled: KeyColor,
    /// For a table entry that carries no color. Default: the stable palette.
    pub uncolored: KeyColor,
    /// For key 0.
    pub zero: ZeroKeyPolicy,
}

impl Default for LabelColorPolicy {
    fn default() -> Self {
        Self {
            unlabeled: KeyColor::Fixed(Rgba::TRANSPARENT),
            uncolored: KeyColor::StablePalette,
            zero: ZeroKeyPolicy::FromTable,
        }
    }
}

impl LabelTable {
    /// The color for `key` under `policy`.
    ///
    /// A table's own color is returned exactly as stored (an `f32` RGBA; a color
    /// read from 8-bit data stays `value / 255`, never re-quantized).
    pub fn color_for_key(&self, key: i64, policy: &LabelColorPolicy) -> Rgba {
        if key == 0 {
            match policy.zero {
                ZeroKeyPolicy::Unlabeled => return policy.unlabeled.resolve(key),
                ZeroKeyPolicy::Fixed(color) => return color,
                ZeroKeyPolicy::FromTable => {}
            }
        }
        match self.get(key) {
            Some(entry) => match entry.rgba {
                Some([r, g, b, a]) => Rgba { r, g, b, a },
                None => policy.uncolored.resolve(key),
            },
            None => policy.unlabeled.resolve(key),
        }
    }
}

/// A label table together with the policy that colors keys it does not color.
/// This is the "color map" for integer label data.
#[derive(Debug, Clone, PartialEq)]
pub struct LabelColorMap {
    table: LabelTable,
    policy: LabelColorPolicy,
}

impl LabelColorMap {
    /// A color map for `table` under `policy`.
    pub fn new(table: LabelTable, policy: LabelColorPolicy) -> Self {
        Self { table, policy }
    }

    /// A color map with the default policy.
    pub fn with_defaults(table: LabelTable) -> Self {
        Self::new(table, LabelColorPolicy::default())
    }

    /// The table (keys, names, and their own colors).
    pub fn table(&self) -> &LabelTable {
        &self.table
    }

    /// The fallback policy.
    pub fn policy(&self) -> &LabelColorPolicy {
        &self.policy
    }

    /// The color for one key.
    pub fn color_for_key(&self, key: i64) -> Rgba {
        self.table.color_for_key(key, &self.policy)
    }

    /// Colors for a list of keys (for example every row of a label column).
    pub fn colors_for_keys(&self, keys: &[i64]) -> Vec<Rgba> {
        keys.iter().map(|&k| self.color_for_key(k)).collect()
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

    fn colored(key: i64, name: &str, c: [u8; 3]) -> LabelEntry {
        let [r, g, b] = c.map(|v| f32::from(v) / 255.0);
        LabelEntry {
            key,
            name: name.into(),
            rgba: Some([r, g, b, 1.0]),
        }
    }

    #[test]
    fn stable_palette_is_periodic_symmetric_and_total() {
        assert_eq!(stable_label_rgb(0), [128, 128, 128]);
        assert_eq!(stable_label_rgb(1), [0, 194, 255]);
        assert_eq!(stable_label_rgb(10), [255, 0, 170]);
        assert_eq!(stable_label_rgb(11), stable_label_rgb(1));
        assert_eq!(stable_label_rgb(-3), stable_label_rgb(3));
        // Never panics, even at the extremes.
        let _ = stable_label_rgb(i64::MIN);
        let _ = stable_label_rgb(i64::MAX);
    }

    #[test]
    fn table_colors_are_returned_exactly_as_stored() {
        let t = LabelTable::new(vec![colored(1, "a", [76, 25, 204])]).unwrap();
        let c = t.color_for_key(1, &LabelColorPolicy::default());
        // Exactly value/255: no re-quantization.
        assert_eq!(
            c.to_array(),
            [76.0 / 255.0, 25.0 / 255.0, 204.0 / 255.0, 1.0]
        );
        // A color that is NOT a multiple of 1/255 also survives untouched.
        let odd = LabelEntry {
            key: 2,
            name: "b".into(),
            rgba: Some([0.1234, 0.5, 0.7777, 0.4]),
        };
        let t2 = LabelTable::new(vec![odd]).unwrap();
        assert_eq!(
            t2.color_for_key(2, &LabelColorPolicy::default()).to_array(),
            [0.1234, 0.5, 0.7777, 0.4]
        );
    }

    #[test]
    fn unlabeled_and_uncolored_follow_the_policy() {
        let uncolored = LabelEntry {
            key: 5,
            name: "no color".into(),
            rgba: None,
        };
        let t = LabelTable::new(vec![colored(1, "a", [255, 0, 0]), uncolored]).unwrap();
        let default = LabelColorPolicy::default();
        assert_eq!(
            t.color_for_key(99, &default),
            Rgba::TRANSPARENT,
            "unlabeled: transparent"
        );
        let [r, g, b] = stable_label_rgb(5);
        assert_eq!(
            t.color_for_key(5, &default),
            Rgba::from_u8(r, g, b, 255),
            "uncolored: palette"
        );
        let gray = Rgba::from_u8(9, 9, 9, 255);
        let custom = LabelColorPolicy {
            unlabeled: KeyColor::Fixed(gray),
            uncolored: KeyColor::Fixed(Rgba::WHITE),
            zero: ZeroKeyPolicy::FromTable,
        };
        assert_eq!(t.color_for_key(99, &custom), gray);
        assert_eq!(t.color_for_key(5, &custom), Rgba::WHITE);
        // Unlabeled keys can use the stable palette too (a dataset with no table).
        let palette = LabelColorPolicy {
            unlabeled: KeyColor::StablePalette,
            ..default
        };
        let [r, g, b] = stable_label_rgb(99);
        assert_eq!(t.color_for_key(99, &palette), Rgba::from_u8(r, g, b, 255));
    }

    #[test]
    fn key_zero_policies() {
        let t = LabelTable::new(vec![
            colored(0, "undefined", [10, 20, 30]),
            colored(1, "a", [255, 0, 0]),
        ])
        .unwrap();
        let from_table = LabelColorPolicy::default();
        assert_eq!(
            t.color_for_key(0, &from_table),
            Rgba::from_u8(10, 20, 30, 255)
        );
        let hide = LabelColorPolicy {
            zero: ZeroKeyPolicy::Unlabeled,
            ..from_table
        };
        assert_eq!(t.color_for_key(0, &hide), Rgba::TRANSPARENT);
        let fixed = LabelColorPolicy {
            zero: ZeroKeyPolicy::Fixed(Rgba::BLACK),
            ..from_table
        };
        assert_eq!(t.color_for_key(0, &fixed), Rgba::BLACK);
        assert_eq!(
            t.color_for_key(1, &hide),
            Rgba::from_u8(255, 0, 0, 255),
            "other keys unaffected"
        );
        // With no entry for 0 at all, FromTable falls through to the unlabeled rule.
        let no_zero = LabelTable::new(vec![colored(1, "a", [255, 0, 0])]).unwrap();
        assert_eq!(no_zero.color_for_key(0, &from_table), Rgba::TRANSPARENT);
    }

    #[test]
    fn label_color_map_keeps_keys_and_colors_a_column() {
        let t = LabelTable::new(vec![
            colored(1007, "fusiform", [1, 2, 3]),
            colored(3, "other", [4, 5, 6]),
        ])
        .unwrap();
        let m = LabelColorMap::with_defaults(t);
        assert_eq!(m.table().get(1007).unwrap().name, "fusiform");
        let colors = m.colors_for_keys(&[1007, 3, 7, 1007]);
        assert_eq!(colors[0], Rgba::from_u8(1, 2, 3, 255));
        assert_eq!(colors[1], Rgba::from_u8(4, 5, 6, 255));
        assert_eq!(colors[2], Rgba::TRANSPARENT);
        assert_eq!(colors[3], colors[0]);
        // Large keys are not squeezed through float on the way.
        let big = (1_i64 << 53) + 1;
        let tb = LabelTable::new(vec![colored(big, "big", [9, 8, 7])]).unwrap();
        let mb = LabelColorMap::with_defaults(tb);
        assert_eq!(mb.color_for_key(big), Rgba::from_u8(9, 8, 7, 255));
        assert_eq!(mb.color_for_key(big - 1), Rgba::TRANSPARENT);
    }
}
