//! Parses secondary structure, and other categories that aren't atoms, from mmCIF files.

use std::{collections::HashMap, io};

use crate::{BackboneSS, SecondaryStructure};

// todo: Save SS to CIF.

#[allow(unused)]
#[derive(Clone, Copy, PartialEq, Debug)]
enum LoopKind {
    None,
    StructConf,
    AtomSite,
    SheetRange,
}

// pub fn load_ss<R: Read + Seek>(mut data: R) -> io::Result<Vec<BackboneSS>> {
// pub fn load_ss<R: Read + Seek>(text: &str) -> io::Result<Vec<BackboneSS>> {
pub fn load_ss(text: &str) -> io::Result<Vec<BackboneSS>> {
    // data.seek(SeekFrom::Start(0))?;
    // let rdr = BufReader::new(data);

    // Caches
    let mut ca_xyz: HashMap<(String, i32), u32> = HashMap::new();
    let mut helix_rows: Vec<(Vec<String>, Vec<String>)> = Vec::new();
    let mut sheet_rows: Vec<(Vec<String>, Vec<String>)> = Vec::new();

    let mut kind = LoopKind::None;
    let mut head: Vec<String> = Vec::new();

    // atom-site column indices (filled on first row)
    let mut a_idx = (None, None, None, None, None, None, None); // asym, seq, atom, x,y,z, id

    let lines: Vec<&str> = text.lines().collect();
    for line in lines {
        // let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }

        if t == "loop_" {
            kind = LoopKind::None;
            head.clear();
            a_idx = (None, None, None, None, None, None, None);
            continue;
        }
        if t == "#" {
            kind = LoopKind::None;
            continue;
        }

        match kind {
            LoopKind::None => {
                if t.starts_with("_struct_conf.") {
                    kind = LoopKind::StructConf;
                    head.push(t.to_owned());
                } else if t.starts_with("_atom_site.") {
                    kind = LoopKind::AtomSite;
                    head.push(t.to_owned());
                } else if t.starts_with("_struct_sheet_range.") {
                    kind = LoopKind::SheetRange;
                    head.push(t.to_owned());
                }
            }

            // ───────────── _struct_conf  (helices/turns/strands) ─────────────
            LoopKind::StructConf => {
                if t.starts_with('_') {
                    head.push(t.to_owned());
                    continue;
                }
                let cols: Vec<String> = t.split_whitespace().map(str::to_owned).collect();
                helix_rows.push((head.clone(), cols));
            }

            // ───────────── _struct_sheet_range (β-strands) ─────────────
            LoopKind::SheetRange => {
                if t.starts_with('_') {
                    head.push(t.to_owned());
                    continue;
                }
                let cols: Vec<String> = t.split_whitespace().map(str::to_owned).collect();
                sheet_rows.push((head.clone(), cols));
            }

            // ───────────── _atom_site (coordinates) ─────────────
            LoopKind::AtomSite => {
                if t.starts_with('_') {
                    head.push(t.to_owned());
                    continue;
                }

                // first data row of atom_site → locate columns
                if a_idx.0.is_none() {
                    for (i, h) in head.iter().enumerate() {
                        match &h[h.rfind('.').unwrap() + 1..] {
                            "label_asym_id" => a_idx.0 = Some(i),
                            "label_seq_id" => a_idx.1 = Some(i),
                            "label_atom_id" => a_idx.2 = Some(i),
                            "id" => a_idx.6 = Some(i),
                            "Cartn_x" => a_idx.3 = Some(i),
                            "Cartn_y" => a_idx.4 = Some(i),
                            "Cartn_z" => a_idx.5 = Some(i),
                            _ => {}
                        }
                    }
                }

                let (ia, isq, iat, _ix, _iy, iz, id) = match a_idx {
                    (Some(a), Some(s), Some(at), Some(x), Some(y), Some(z), Some(id)) => {
                        (a, s, at, x, y, z, id)
                    }
                    _ => continue,
                };

                let c: Vec<&str> = t.split_whitespace().collect();
                if c.len() <= iz || c[iat] != "CA" {
                    continue;
                }

                if let (Ok(seq), Ok(serial)) = (c[isq].parse::<i32>(), c[id].parse::<u32>()) {
                    ca_xyz.insert((c[ia].to_owned(), seq), serial);
                }
            }
        }
    }

    let mut ss = Vec::new();

    // Helices from _struct_conf -----
    for (h, c) in helix_rows {
        // resolve indices once per header set
        fn find(h: &[String], tag: &str) -> Option<usize> {
            h.iter().position(|s| s.ends_with(tag))
        }
        let i_type = find(&h, "conf_type_id");
        let i_ba = find(&h, "beg_label_asym_id");
        let i_bs = find(&h, "beg_label_seq_id");
        let i_ea = find(&h, "end_label_asym_id");
        let i_es = find(&h, "end_label_seq_id");
        let (i_type, i_ba, i_bs, i_ea, i_es) = match (i_type, i_ba, i_bs, i_ea, i_es) {
            (Some(a), Some(b), Some(c), Some(d), Some(e)) => (a, b, c, d, e),
            _ => continue,
        };

        if !c[i_type].starts_with("HELX") {
            continue;
        }

        let beg_seq = c[i_bs].parse().ok();
        let end_seq = c[i_es].parse().ok();
        if beg_seq.is_none() || end_seq.is_none() {
            continue;
        }

        let start_sn = match ca_xyz.get(&(c[i_ba].clone(), beg_seq.unwrap())) {
            Some(v) => *v,
            None => continue,
        };
        let end_sn = match ca_xyz.get(&(c[i_ea].clone(), end_seq.unwrap())) {
            Some(v) => *v,
            None => continue,
        };

        ss.push(BackboneSS {
            start_sn,
            end_sn,
            sec_struct: SecondaryStructure::Helix,
        });
    }

    // ----- β-strands from _struct_sheet_range -----
    for (h, c) in sheet_rows {
        fn idx(h: &[String], tag: &str) -> Option<usize> {
            h.iter().position(|s| s.ends_with(tag))
        }
        let ib_a = idx(&h, "beg_label_asym_id");
        let ib_s = idx(&h, "beg_label_seq_id");
        let ie_a = idx(&h, "end_label_asym_id");
        let ie_s = idx(&h, "end_label_seq_id");
        let (ib_a, ib_s, ie_a, ie_s) = match (ib_a, ib_s, ie_a, ie_s) {
            (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
            _ => continue,
        };

        let beg_seq = c[ib_s].parse().ok();
        let end_seq = c[ie_s].parse().ok();
        if beg_seq.is_none() || end_seq.is_none() {
            continue;
        }

        let start_sn = match ca_xyz.get(&(c[ib_a].clone(), beg_seq.unwrap())) {
            Some(v) => *v,
            None => continue,
        };
        let end_sn = match ca_xyz.get(&(c[ie_a].clone(), end_seq.unwrap())) {
            Some(v) => *v,
            None => continue,
        };

        ss.push(BackboneSS {
            start_sn,
            end_sn,
            sec_struct: SecondaryStructure::Sheet,
        });
    }

    Ok(ss)
}

/// One row of a category: item names (e.g. `db_name`, lowercase), to values.
pub(crate) type CifRow = HashMap<String, String>;

/// Rows of the named categories (e.g. `_struct_ref`), whether written as a loop or as key-value
/// items. Unknown (`?`) and inapplicable (`.`) values are omitted.
pub(crate) fn load_categories(text: &str, names: &[&str]) -> HashMap<String, Vec<CifRow>> {
    let mut result: HashMap<String, Vec<CifRow>> = HashMap::new();
    // Key-value items are a category's single row.
    let mut kv_rows: HashMap<String, CifRow> = HashMap::new();

    // The category and item name of a tag, if it's one we're reading. E.g. `_struct_ref.db_name`
    // -> (`_struct_ref`, `db_name`).
    let split = |tag: &str| -> Option<(&str, String)> {
        let (cat, item) = tag.split_once('.')?;
        let cat = names.iter().find(|n| n.eq_ignore_ascii_case(cat))?;
        Some((*cat, item.to_ascii_lowercase()))
    };

    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;

    while i < lines.len() {
        // Skip over text fields, so we don't mistake their contents for tags.
        if lines[i].starts_with(';') {
            i = read_text_field(&lines, i).1;
            continue;
        }

        let line = lines[i].trim();

        if line.eq_ignore_ascii_case("loop_") {
            let mut j = i + 1;
            let mut tags = Vec::new();
            while j < lines.len() && lines[j].trim_start().starts_with('_') {
                tags.push(lines[j].trim());
                j += 1;
            }

            let Some(cat) = tags.first().and_then(|t| split(t)).map(|(cat, _)| cat) else {
                // Not a category we're reading; its rows are skipped over line by line.
                i = j;
                continue;
            };
            let items: Vec<String> = tags
                .iter()
                .filter_map(|t| t.split_once('.'))
                .map(|(_, item)| item.to_ascii_lowercase())
                .collect();

            let (values, end) = read_values(&lines, j);
            let rows = result.entry(cat.to_owned()).or_default();
            for chunk in values.chunks_exact(items.len()) {
                rows.push(
                    items
                        .iter()
                        .zip(chunk)
                        .filter(|(_, v)| !is_null(v))
                        .map(|(item, v)| (item.clone(), v.clone()))
                        .collect(),
                );
            }

            i = end;
            continue;
        }

        if line.starts_with('_')
            && let Some(tag) = line.split_whitespace().next()
            && let Some((cat, item)) = split(tag)
        {
            // The value is either on the same line, or the next: e.g. a long quoted one, or a
            // text field.
            let (value, next) = match tokenize_line(&line[tag.len()..]).into_iter().next() {
                Some(v) => (v, i + 1),
                None if lines.get(i + 1).is_some_and(|l| l.starts_with(';')) => {
                    read_text_field(&lines, i + 1)
                }
                None => (
                    lines
                        .get(i + 1)
                        .and_then(|l| tokenize_line(l).into_iter().next())
                        .unwrap_or_default(),
                    i + 2,
                ),
            };

            if !is_null(&value) {
                kv_rows
                    .entry(cat.to_owned())
                    .or_default()
                    .insert(item, value);
            }

            i = next;
            continue;
        }

        i += 1;
    }

    for (cat, row) in kv_rows {
        result.entry(cat).or_default().push(row);
    }

    result
}

/// `?` and `.` are mmCIF's markers for unknown and inapplicable values.
fn is_null(value: &str) -> bool {
    matches!(value, "" | "?" | ".")
}

/// Whether a line starts with a reserved word, e.g. `loop_` or `data_`. These, and tags, end a
/// loop's values.
fn is_reserved(line: &str) -> bool {
    ["loop_", "data_", "save_", "global_", "stop_"]
        .iter()
        .any(|w| {
            line.as_bytes()
                .get(..w.len())
                .is_some_and(|b| b.eq_ignore_ascii_case(w.as_bytes()))
        })
}

/// A loop's values, from line `i` up to the next tag or reserved word. Returns them, and the
/// index of the line after.
fn read_values(lines: &[&str], mut i: usize) -> (Vec<String>, usize) {
    let mut result = Vec::new();

    while i < lines.len() {
        if lines[i].starts_with(';') {
            let (value, next) = read_text_field(lines, i);
            result.push(value);
            i = next;
            continue;
        }

        let line = lines[i].trim_start();
        if line.starts_with('_') || is_reserved(line) {
            break;
        }

        result.extend(tokenize_line(line));
        i += 1;
    }

    (result, i)
}

/// A text field: a value spanning lines, from one starting with `;` to the next. Returns it
/// without its delimiters, and the index of the line after it.
fn read_text_field(lines: &[&str], i: usize) -> (String, usize) {
    let mut result = lines[i][1..].to_owned();
    let mut j = i + 1;

    while j < lines.len() && !lines[j].starts_with(';') {
        result.push('\n');
        result.push_str(lines[j]);
        j += 1;
    }

    (result.trim().to_owned(), j + 1)
}

/// Split a line into values, removing quotes, and stopping at a comment. A quoted value ends at
/// its quote only when followed by whitespace, so it may contain that quote, e.g. `'N-(2'-OH)'`.
fn tokenize_line(line: &str) -> Vec<String> {
    let mut result = Vec::new();
    // Slicing at these byte indices is safe: the delimiters we split on are all ASCII.
    let b = line.as_bytes();
    let mut i = 0;

    while i < b.len() {
        match b[i] {
            c if c.is_ascii_whitespace() => i += 1,
            b'#' => break,
            q @ (b'\'' | b'"') => {
                let start = i + 1;
                let mut end = start;
                while end < b.len()
                    && !(b[end] == q && b.get(end + 1).is_none_or(|c| c.is_ascii_whitespace()))
                {
                    end += 1;
                }
                result.push(line[start..end].to_owned());
                i = end + 1;
            }
            _ => {
                let start = i;
                while i < b.len() && !b[i].is_ascii_whitespace() {
                    i += 1;
                }
                result.push(line[start..i].to_owned());
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenizes_quoted_values() {
        assert_eq!(
            tokenize_line("EMDB 'A spike, one RBD up' EMD-21375 'N-(2'-OH)' # comment"),
            vec!["EMDB", "A spike, one RBD up", "EMD-21375", "N-(2'-OH)"]
        );
    }

    /// A loop with text fields and rows on one line, as in 6GO7, and key-value items with a value
    /// on the next line, as in 7K3G.
    #[test]
    fn reads_loops_and_key_value_items() {
        let text = "data_TEST
loop_
_struct_ref.id
_struct_ref.db_name
_struct_ref.db_code
_struct_ref.pdbx_db_accession
_struct_ref.pdbx_db_isoform
_struct_ref.entity_id
_struct_ref.pdbx_seq_one_letter_code
_struct_ref.pdbx_align_begin
1 UNP TDT_MOUSE   P09838 ?        1
;SPSPVPGSQNVPAPAVKKISQYACQRRTTLNNYNQLFTDALDILAENDELRENEGSCLAFMRASSVLKSLPFPITSMKDT
_EGIPCLGDKVKSIIEGIIEDGESSEAKAV
;
132
2 UNP DPOLM_MOUSE Q9JIW4 ?        1 HQYHRSHLADSAHNLRQRSSTMDAFERSFC 363
4 PDB 6GO7        6GO7   ?        2 ? 1
#
_pdbx_database_related.db_name        BMRB
_pdbx_database_related.details
'SARS-CoV-2 Envelope Protein Transmembrane Domain'
_pdbx_database_related.db_id          30795
_pdbx_database_related.content_type   unspecified
#
";
        let cats = load_categories(text, &["_struct_ref", "_pdbx_database_related"]);

        let refs = &cats["_struct_ref"];
        assert_eq!(refs.len(), 3);
        assert_eq!(refs[0]["pdbx_db_accession"], "P09838");
        assert!(refs[0]["pdbx_seq_one_letter_code"].contains("\n_EGIPCLG"));
        assert!(!refs[0].contains_key("pdbx_db_isoform"));
        assert_eq!(refs[1]["pdbx_db_accession"], "Q9JIW4");
        assert_eq!(refs[2]["db_name"], "PDB");

        let related = &cats["_pdbx_database_related"];
        assert_eq!(related.len(), 1);
        assert_eq!(related[0]["db_id"], "30795");
        assert_eq!(
            related[0]["details"],
            "SARS-CoV-2 Envelope Protein Transmembrane Domain"
        );
    }
}
