//! For reading and writing GenBank (and GenPept) flat files: one or more records, each with a
//! header, a feature table, and the sequence, and ending with `//`.
//!
//! We read the sequence and the main header fields into the sequence's name and metadata. The
//! feature table and references are skipped for now. Writing produces a minimal record from the
//! same fields, with an empty feature table.
//!
//! [Format reference](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
//!
//! todo: See also `plascad`'s genbank implementaiton

use std::{
    collections::HashMap,
    fs,
    fs::File,
    io::{self, BufWriter, ErrorKind, Write},
    path::Path,
    time::{SystemTime, UNIX_EPOCH},
};

use na_seq::{SEQ_DESCRIPTION_KEY, SeqType, Sequence, SequenceData};

// Metadata keys for the header fields we read and write. The `DEFINITION` line is stored under
// `na_seq::SEQ_DESCRIPTION_KEY`.
pub const KEY_ACCESSION: &str = "Accession";
pub const KEY_VERSION: &str = "Version";
pub const KEY_KEYWORDS: &str = "Keywords";
pub const KEY_SOURCE: &str = "Source";
pub const KEY_ORGANISM: &str = "Organism";
/// The lineage lines following `ORGANISM`.
pub const KEY_TAXONOMY: &str = "Taxonomy";
pub const KEY_DBLINK: &str = "DB link";
pub const KEY_COMMENT: &str = "Comment";
// These are from the `LOCUS` line.
/// E.g. "mRNA", or "ss-DNA".
pub const KEY_MOL_TYPE: &str = "Molecule type";
/// "linear" or "circular".
pub const KEY_TOPOLOGY: &str = "Topology";
/// A three-letter code, e.g. "PLN".
pub const KEY_DIVISION: &str = "Division";
/// E.g. "21-JUN-1999"
pub const KEY_DATE: &str = "Date";

/// Header values start at this column.
const INDENT: usize = 12;
/// Header lines are wrapped to this width when writing.
const LINE_WIDTH: usize = 79;
const RESIDUES_PER_LINE: usize = 60;
const RESIDUES_PER_GROUP: usize = 10;

const MONTHS: [&str; 12] = [
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
];

#[derive(Clone, Debug, Default)]
pub struct GenBank {
    pub records: Vec<Sequence>,
    /// When loading: the number of residue letters left out across all records, as they can't
    /// be represented. E.g. ambiguity codes like N or X.
    pub skipped_residues: usize,
}

/// The header field currently being read; continuation lines are appended to it.
#[derive(Clone, Copy, PartialEq)]
enum Field {
    /// A header value stored under a metadata key. `true` to join its lines with newlines instead
    /// of spaces.
    Meta(&'static str, bool),
    Taxonomy,
    /// Anything we don't store, e.g. references and the feature table.
    Skip,
    Origin,
}

/// One record as it's being read.
#[derive(Default)]
struct RecordRaw {
    name: String,
    /// From the `LOCUS` line.
    seq_type: Option<SeqType>,
    metadata: HashMap<String, String>,
    residues: String,
}

impl RecordRaw {
    fn append(&mut self, key: &str, text: &str, newline: bool) {
        let text = text.trim();
        if text.is_empty() {
            return;
        }

        let entry = self.metadata.entry(key.to_owned()).or_default();
        if !entry.is_empty() {
            entry.push(if newline { '\n' } else { ' ' });
        }
        entry.push_str(text);
    }

    fn parse_locus(&mut self, value: &str) {
        let mut tokens = value.split_whitespace();
        self.name = tokens.next().unwrap_or_default().to_owned();

        let mut is_protein = false;
        for token in tokens {
            let lower = token.to_ascii_lowercase();

            if lower == "aa" {
                is_protein = true;
            } else if lower == "linear" || lower == "circular" {
                self.metadata.insert(KEY_TOPOLOGY.to_owned(), lower);
            } else if lower.contains("dna") || lower.contains("rna") {
                self.metadata
                    .insert(KEY_MOL_TYPE.to_owned(), token.to_owned());
            } else if token.len() == 11 && token.matches('-').count() == 2 {
                self.metadata.insert(KEY_DATE.to_owned(), token.to_owned());
            } else if token.len() == 3 && token.chars().all(|c| c.is_ascii_uppercase()) {
                self.metadata
                    .insert(KEY_DIVISION.to_owned(), token.to_owned());
            }
            // Otherwise, e.g. the length, which we get from the sequence itself.
        }

        self.seq_type = if is_protein {
            Some(SeqType::AminoAcid)
        } else {
            match self.metadata.get(KEY_MOL_TYPE) {
                Some(t) if t.to_ascii_lowercase().contains("rna") => Some(SeqType::Rna),
                Some(_) => Some(SeqType::Dna),
                None => None,
            }
        };
    }

    fn finish(mut self, skipped_residues: &mut usize) -> Sequence {
        // A value of "." means there is none, e.g. for keywords.
        self.metadata.retain(|_, v| v != ".");

        let seq_type = self
            .seq_type
            .unwrap_or_else(|| SeqType::infer(&self.residues));
        let (data, skipped) = SequenceData::from_letters(&self.residues, seq_type);
        *skipped_residues += skipped;

        let mut result = Sequence::new(data, self.name);
        result.metadata = self.metadata;
        result
    }
}

impl GenBank {
    pub fn new(text: &str) -> io::Result<Self> {
        let mut records = Vec::new();
        let mut skipped_residues = 0;

        let mut current: Option<RecordRaw> = None;
        let mut field = Field::Skip;

        for line in text.lines() {
            let line = line.trim_end();

            if line.starts_with("//") {
                if let Some(record) = current.take() {
                    records.push(record.finish(&mut skipped_residues));
                }
                field = Field::Skip;
                continue;
            }

            if line.trim().is_empty() {
                continue;
            }

            if field == Field::Origin && line.starts_with(' ') {
                if let Some(record) = &mut current {
                    // Includes the position numbers and spaces; these are filtered when parsing.
                    record.residues.push_str(line);
                }
                continue;
            }

            let (keyword, value) = split_keyword(line);

            if keyword == "LOCUS" {
                // A new record, even if the previous one was missing its terminator.
                if let Some(record) = current.take() {
                    records.push(record.finish(&mut skipped_residues));
                }

                let mut record = RecordRaw::default();
                record.parse_locus(value);
                current = Some(record);
                field = Field::Skip;
                continue;
            }

            let Some(record) = &mut current else {
                continue;
            };

            if keyword.is_empty() {
                // A continuation of the current field.
                match field {
                    Field::Meta(key, newline) => record.append(key, value, newline),
                    Field::Taxonomy => record.append(KEY_TAXONOMY, value, false),
                    Field::Skip | Field::Origin => (),
                }
                continue;
            }

            if keyword == "ORGANISM" {
                record.append(KEY_ORGANISM, value, false);
                // The lines that follow are the lineage.
                field = Field::Taxonomy;
                continue;
            }

            field = match keyword {
                "DEFINITION" => Field::Meta(SEQ_DESCRIPTION_KEY, false),
                "ACCESSION" => Field::Meta(KEY_ACCESSION, false),
                "VERSION" => Field::Meta(KEY_VERSION, false),
                "KEYWORDS" => Field::Meta(KEY_KEYWORDS, false),
                "SOURCE" => Field::Meta(KEY_SOURCE, false),
                "DBLINK" => Field::Meta(KEY_DBLINK, true),
                "COMMENT" => Field::Meta(KEY_COMMENT, true),
                "ORIGIN" => Field::Origin,
                _ => Field::Skip,
            };

            if let Field::Meta(key, newline) = field {
                record.append(key, value, newline);
            }
        }

        if let Some(record) = current.take() {
            records.push(record.finish(&mut skipped_residues));
        }

        if records.is_empty() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "No LOCUS records found in the GenBank data",
            ));
        }

        Ok(Self {
            records,
            skipped_residues,
        })
    }

    pub fn load(path: &Path) -> io::Result<Self> {
        let text = fs::read_to_string(path)?;

        let mut result = Self::new(&text)?;
        for record in &mut result.records {
            record.path = Some(path.to_owned());
        }

        Ok(result)
    }

    pub fn to_text(&self) -> String {
        let mut result = String::new();

        for record in &self.records {
            write_record(record, &mut result);
        }

        result
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = BufWriter::new(File::create(path)?);
        file.write_all(self.to_text().as_bytes())?;
        file.flush()
    }
}

/// Split a header line into its keyword and value. The keyword is empty for continuation lines.
/// Sub-keywords, e.g. `  ORGANISM`, are returned without their indent.
fn split_keyword(line: &str) -> (&str, &str) {
    let (head, value) = if line.len() > INDENT && line.is_char_boundary(INDENT) {
        line.split_at(INDENT)
    } else {
        (line, "")
    };

    let keyword = head.trim();
    // A keyword longer than the indent, e.g. in a malformed file, or a line with no value.
    if keyword.contains(char::is_whitespace) {
        let trimmed = line.trim_start();
        return match trimmed.split_once(char::is_whitespace) {
            Some((k, v)) => (k, v.trim_start()),
            None => (trimmed, ""),
        };
    }

    (keyword, value)
}

/// Write a header field, wrapping its value, with continuation lines indented.
fn write_field(keyword: &str, value: &str, out: &mut String) {
    let mut first = true;

    // Kept verbatim if it fits, e.g. to preserve the spacing in `VERSION` values.
    if !value.trim().is_empty() && !value.contains('\n') && INDENT + value.len() <= LINE_WIDTH {
        push_header_line(keyword, value.trim(), &mut first, out);
        return;
    }

    for paragraph in value.lines() {
        let mut line = String::new();

        for word in paragraph.split_whitespace() {
            if !line.is_empty() && INDENT + line.len() + 1 + word.len() > LINE_WIDTH {
                push_header_line(keyword, &line, &mut first, out);
                line.clear();
            }

            if !line.is_empty() {
                line.push(' ');
            }
            line.push_str(word);
        }

        push_header_line(keyword, &line, &mut first, out);
    }

    if first {
        // No value.
        push_header_line(keyword, ".", &mut first, out);
    }
}

fn push_header_line(keyword: &str, text: &str, first: &mut bool, out: &mut String) {
    let lead = if *first { keyword } else { "" };
    *first = false;

    out.push_str(&format!("{lead:<width$}{text}\n", width = INDENT));
}

fn write_record(record: &Sequence, out: &mut String) {
    let meta = |key: &str| {
        record
            .metadata
            .get(key)
            .map(String::as_str)
            .filter(|v| !v.trim().is_empty())
    };

    let seq_type = record.seq_type();
    let name = if record.name.trim().is_empty() {
        "unnamed".to_owned()
    } else {
        record.name.trim().replace(char::is_whitespace, "_")
    };

    let (unit, mol_type) = match seq_type {
        SeqType::AminoAcid => ("aa", ""),
        SeqType::Dna => ("bp", meta(KEY_MOL_TYPE).unwrap_or("DNA")),
        SeqType::Rna => ("bp", meta(KEY_MOL_TYPE).unwrap_or("RNA")),
    };

    // E.g. "ss-DNA": The strandedness prefix has its own columns.
    let (strand, mol_type) = mol_type.split_at(mol_type.find('-').map_or(0, |i| i + 1));

    let topology = meta(KEY_TOPOLOGY).unwrap_or("linear");
    let division = meta(KEY_DIVISION).unwrap_or("UNA");
    let date = meta(KEY_DATE).map(str::to_owned).unwrap_or_else(date_today);

    // Column positions follow the GenBank release notes, section 3.4.4.
    out.push_str(&format!(
        "LOCUS       {name:<16} {:>11} {unit} {strand:>3}{mol_type:<6}  {topology:<8} {division} {date}\n",
        record.data.len(),
    ));

    write_field("DEFINITION", record.description().unwrap_or("."), out);
    write_field("ACCESSION", meta(KEY_ACCESSION).unwrap_or(&name), out);

    if let Some(v) = meta(KEY_VERSION) {
        write_field("VERSION", v, out);
    }
    if let Some(v) = meta(KEY_DBLINK) {
        write_field("DBLINK", v, out);
    }

    write_field("KEYWORDS", meta(KEY_KEYWORDS).unwrap_or("."), out);
    write_field("SOURCE", meta(KEY_SOURCE).unwrap_or("."), out);
    write_field("  ORGANISM", meta(KEY_ORGANISM).unwrap_or("."), out);

    if let Some(v) = meta(KEY_TAXONOMY) {
        write_field("", v, out);
    }
    if let Some(v) = meta(KEY_COMMENT) {
        write_field("COMMENT", v, out);
    }

    out.push_str("FEATURES             Location/Qualifiers\n");
    out.push_str("ORIGIN\n");

    let letters = record.data.to_letters().to_ascii_lowercase();
    for (line_i, line) in letters.as_bytes().chunks(RESIDUES_PER_LINE).enumerate() {
        out.push_str(&format!("{:>9}", line_i * RESIDUES_PER_LINE + 1));

        for group in line.chunks(RESIDUES_PER_GROUP) {
            out.push(' ');
            // Residue letters are ASCII.
            out.push_str(std::str::from_utf8(group).unwrap_or_default());
        }
        out.push('\n');
    }

    out.push_str("//\n");
}

/// Today's date (UTC) in GenBank's format, e.g. "26-SEP-2026".
fn date_today() -> String {
    let days = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() / 86_400)
        .unwrap_or_default() as i64;

    // Civil date from days since the Unix epoch; Howard Hinnant's `civil_from_days`.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + if month <= 2 { 1 } else { 0 };

    format!("{day:02}-{}-{year}", MONTHS[(month - 1) as usize])
}
