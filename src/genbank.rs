//! For reading and writing GenBank (and GenPept) flat files: one or more records, each with a
//! header, a feature table, and the sequence, ending with `//`.
//!
//! Parsing and writing is done by the [gb_io](https://docs.rs/gb-io) crate; we convert between
//! its records and [`Sequence`]. Header fields become metadata, the feature table becomes
//! features, and references and comments are kept alongside, in [`GenBankRecord`].
//!
//! `gb_io` only handles nucleotide records; we adapt protein (GenPept) records' `LOCUS` lines, which
//! use "aa" in place of "bp", on the way in and out.
//!
//! [Format reference](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)

use std::{
    borrow::Cow,
    fs,
    io::{self, ErrorKind},
    path::Path,
    time::{SystemTime, UNIX_EPOCH},
};

use gb_io::{
    reader::SeqReader,
    seq::{After, Before, Date, Feature, Location, Reference, Source, Topology},
    writer::SeqWriter,
};
use na_seq::{SEQ_DESCRIPTION_KEY, SeqFeature, SeqRange, SeqTopology, SeqType, Sequence, Strand};

// Metadata keys for the header fields we read and write. The `DEFINITION` line is stored under
// `na_seq::SEQ_DESCRIPTION_KEY`. Multi-line values keep their line breaks.
pub const KEY_ACCESSION: &str = "Accession";
pub const KEY_VERSION: &str = "Version";
pub const KEY_KEYWORDS: &str = "Keywords";
pub const KEY_SOURCE: &str = "Source";
/// The organism's name, followed by its lineage on the lines after.
pub const KEY_ORGANISM: &str = "Organism";
pub const KEY_DBLINK: &str = "DB link";
// These are from the `LOCUS` line.
/// E.g. "mRNA", or "ss-DNA".
pub const KEY_MOL_TYPE: &str = "Molecule type";
/// A three-letter code, e.g. "PLN".
pub const KEY_DIVISION: &str = "Division";
/// E.g. "21-JUN-1999".
pub const KEY_DATE: &str = "Date";

/// The division code we write when none is set: "unannotated".
const DIVISION_DEFAULT: &str = "UNA";
/// `gb_io`'s placeholder for a missing division. Not stored as metadata.
const DIVISION_GB_IO_DEFAULT: &str = "UNK";

/// Put in place of a protein record's molecule type, so `gb_io` can parse its `LOCUS` line.
const PROTEIN_MOL_TYPE: &str = "PROTEIN";

/// These qualifiers aren't standard GenBank, but are written by e.g. SnapGene and ApE. We read
/// them into the feature's label and strand, and write them back.
const QUAL_LABEL: &str = "label";
const QUAL_DIRECTION: &str = "direction";

const MONTHS: [&str; 12] = [
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
];

/// A publication cited by a record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenBankReference {
    /// E.g. "1  (bases 1 to 5028)".
    pub description: String,
    pub authors: Option<String>,
    pub consortium: Option<String>,
    pub title: String,
    pub journal: Option<String>,
    pub pubmed: Option<String>,
    pub remark: Option<String>,
}

impl From<&Reference> for GenBankReference {
    fn from(r: &Reference) -> Self {
        Self {
            description: r.description.clone(),
            authors: r.authors.clone(),
            consortium: r.consortium.clone(),
            title: r.title.clone(),
            journal: r.journal.clone(),
            pubmed: r.pubmed.clone(),
            remark: r.remark.clone(),
        }
    }
}

impl From<&GenBankReference> for Reference {
    fn from(r: &GenBankReference) -> Self {
        Self {
            description: r.description.clone(),
            authors: r.authors.clone(),
            consortium: r.consortium.clone(),
            title: r.title.clone(),
            journal: r.journal.clone(),
            pubmed: r.pubmed.clone(),
            remark: r.remark.clone(),
        }
    }
}

/// One record: the sequence with its header fields and features, and the parts of the header
/// that don't fit in its metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct GenBankRecord {
    pub seq: Sequence,
    pub references: Vec<GenBankReference>,
    /// Each `COMMENT` field. Line breaks are kept.
    pub comments: Vec<String>,
}

impl From<Sequence> for GenBankRecord {
    fn from(seq: Sequence) -> Self {
        Self {
            seq,
            references: Vec::new(),
            comments: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct GenBank {
    pub records: Vec<GenBankRecord>,
    /// When loading: the number of residue letters left out across all records, as they can't
    /// be represented. E.g. ambiguity codes like N or X. Feature ranges are adjusted to match.
    pub skipped_residues: usize,
}

impl GenBank {
    pub fn new(text: &str) -> io::Result<Self> {
        let (text, is_protein) = adapt_protein_locus_lines(text);

        let gb_records = SeqReader::new(text.as_bytes())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, format!("GenBank: {e}")))?;

        if gb_records.is_empty() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "No records found in the GenBank data",
            ));
        }

        let mut records = Vec::with_capacity(gb_records.len());
        let mut skipped_residues = 0;

        for (i, gb) in gb_records.iter().enumerate() {
            let is_protein = is_protein.get(i).copied().unwrap_or_default();
            let (record, skipped) = record_from_gb(gb, is_protein);

            records.push(record);
            skipped_residues += skipped;
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
            record.seq.path = Some(path.to_owned());
        }

        Ok(result)
    }

    pub fn to_text(&self) -> io::Result<String> {
        let mut result = String::new();

        for record in &self.records {
            let gb = record_to_gb(record);

            let mut buf = Vec::new();
            SeqWriter::new(&mut buf).write(&gb)?;

            let mut text = String::from_utf8(buf)
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e.to_string()))?;

            if record.seq.seq_type() == SeqType::AminoAcid {
                text = protein_locus_line(&text);
            }

            result.push_str(&text);
        }

        Ok(result)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        fs::write(path, self.to_text()?)
    }
}

/// `gb_io` only parses `LOCUS` lines with "bp" units. For each record's `LOCUS` line with "aa"
/// units, i.e. a protein, substitute one it can parse. Returns the adapted text, and whether each
/// record, in order, is a protein.
fn adapt_protein_locus_lines(text: &str) -> (Cow<'_, str>, Vec<bool>) {
    let mut is_protein = Vec::new();
    let mut adapted = String::new();
    let mut any_protein = false;

    for line in text.split_inclusive('\n') {
        let Some(value) = line.strip_prefix("LOCUS") else {
            adapted.push_str(line);
            continue;
        };

        let tokens: Vec<&str> = value.split_whitespace().collect();
        let protein = tokens.get(2).is_some_and(|t| t.eq_ignore_ascii_case("aa"));
        is_protein.push(protein);

        if !protein {
            adapted.push_str(line);
            continue;
        }
        any_protein = true;

        // Name, length, "aa", then e.g. "linear PLN 12-APR-1996". NCBI's protein records all
        // state their topology, but we don't depend on it.
        let rest = &tokens[3..];
        let topology = match rest.first() {
            Some(&"linear") | Some(&"circular") => "",
            _ => "linear ",
        };

        adapted.push_str(&format!(
            "LOCUS       {} {} bp {PROTEIN_MOL_TYPE} {topology}{}\n",
            tokens[0],
            tokens[1],
            rest.join(" ")
        ));
    }

    if any_protein {
        (Cow::Owned(adapted), is_protein)
    } else {
        (Cow::Borrowed(text), is_protein)
    }
}

/// `gb_io` writes "bp" units on the `LOCUS` line. Use "aa", for a protein record.
fn protein_locus_line(record_text: &str) -> String {
    let (locus, rest) = record_text.split_once('\n').unwrap_or((record_text, ""));

    // The name can't contain spaces, so this is the units.
    let locus = locus.replacen(" bp ", " aa ", 1);

    format!("{locus}\n{rest}")
}

fn record_from_gb(gb: &gb_io::seq::Seq, is_protein: bool) -> (GenBankRecord, usize) {
    let letters = String::from_utf8_lossy(&gb.seq);

    let mol_type = gb.molecule_type.as_deref().filter(|_| !is_protein);

    let seq_type = if is_protein {
        SeqType::AminoAcid
    } else {
        match mol_type {
            Some(t) if t.to_ascii_uppercase().contains("RNA") => SeqType::Rna,
            Some(_) => SeqType::Dna,
            None => SeqType::infer(&letters),
        }
    };

    let circular = gb.topology == Topology::Circular && !is_protein;

    let features = gb
        .features
        .iter()
        .map(|f| feature_from_gb(f, circular, gb.seq.len()))
        .collect();
    let name = gb.name.clone().unwrap_or_default();

    let (mut seq, skipped) =
        Sequence::from_letters_with_features(&letters, seq_type, name, features);

    if !is_protein {
        seq.topology = Some(match gb.topology {
            Topology::Linear => SeqTopology::Linear,
            Topology::Circular => SeqTopology::Circular,
        });
    }

    let mut meta = |key: &str, val: Option<&str>| {
        if let Some(v) = val.map(str::trim).filter(|v| !v.is_empty() && *v != ".") {
            seq.metadata.insert(key.to_owned(), v.to_owned());
        }
    };

    meta(SEQ_DESCRIPTION_KEY, gb.definition.as_deref());
    meta(KEY_ACCESSION, gb.accession.as_deref());
    meta(KEY_VERSION, gb.version.as_deref());
    meta(KEY_DBLINK, gb.dblink.as_deref());
    meta(KEY_KEYWORDS, gb.keywords.as_deref());
    meta(KEY_MOL_TYPE, mol_type);

    if let Some(source) = &gb.source {
        meta(KEY_SOURCE, Some(&source.source));
        meta(KEY_ORGANISM, source.organism.as_deref());
    }

    if gb.division != DIVISION_GB_IO_DEFAULT {
        meta(KEY_DIVISION, Some(&gb.division));
    }

    let date = gb.date.as_ref().map(ToString::to_string);
    meta(KEY_DATE, date.as_deref());

    let record = GenBankRecord {
        seq,
        references: gb.references.iter().map(Into::into).collect(),
        comments: gb.comments.clone(),
    };

    (record, skipped)
}

fn feature_from_gb(feature: &Feature, circular: bool, seq_len: usize) -> SeqFeature {
    let mut ranges = Vec::new();
    let mut complement = false;
    add_location_ranges(&feature.location, &mut ranges, &mut complement);

    // A feature spanning the origin of a circular sequence is written as e.g.
    // `join(120..130,1..10)`; represent it as one range that wraps.
    if circular {
        let mut merged: Vec<SeqRange> = Vec::with_capacity(ranges.len());
        for r in ranges {
            match merged.last_mut() {
                Some(prev) if prev.end == seq_len && r.start == 1 && prev.start > r.end => {
                    prev.end = r.end;
                }
                _ => merged.push(r),
            }
        }
        ranges = merged;
    }

    let mut strand = if complement {
        Strand::Reverse
    } else {
        Strand::None
    };

    let mut label = String::new();
    let mut qualifiers = Vec::new();

    for (key, val) in &feature.qualifiers {
        let val = val.clone().unwrap_or_default();

        match key.as_ref() {
            QUAL_LABEL if label.is_empty() => label = val,
            QUAL_DIRECTION => match val.to_ascii_lowercase().as_str() {
                "right" => strand = Strand::Forward,
                "left" => strand = Strand::Reverse,
                _ => qualifiers.push((key.to_string(), val)),
            },
            _ => qualifiers.push((key.to_string(), val)),
        }
    }

    SeqFeature {
        kind: feature.kind.to_string(),
        label,
        ranges,
        strand,
        color: None,
        qualifiers,
    }
}

/// Flatten a location into ranges. `gb_io` positions are 0-based, with exclusive ends; ours are
/// 1-based and inclusive.
fn add_location_ranges(loc: &Location, ranges: &mut Vec<SeqRange>, complement: &mut bool) {
    match loc {
        // `end` before `start` wraps past the origin of a circular sequence. Not standard GenBank,
        // but some programs write it.
        Location::Range((start, _), (end, _)) => {
            if *start >= 0 && *end >= 0 && end != start {
                ranges.push(SeqRange::new(*start as usize + 1, *end as usize));
            }
        }
        // A site between two adjacent positions, e.g. a cut site. `gb_io` stores the two
        // positions, 0-based.
        Location::Between(a, b) => {
            if *a >= 0 && *b >= 0 {
                ranges.push(SeqRange::new(*a as usize + 1, *b as usize + 1));
            }
        }
        Location::Complement(inner) => {
            *complement = true;
            add_location_ranges(inner, ranges, complement);
        }
        Location::Join(locs)
        | Location::Order(locs)
        | Location::Bond(locs)
        | Location::OneOf(locs) => {
            for l in locs {
                add_location_ranges(l, ranges, complement);
            }
        }
        // In other records, or of unknown length.
        Location::External(..) | Location::Gap(_) => (),
    }
}

fn record_to_gb(record: &GenBankRecord) -> gb_io::seq::Seq {
    let seq = &record.seq;

    let meta = |key: &str| {
        seq.metadata
            .get(key)
            .map(String::as_str)
            .filter(|v| !v.trim().is_empty())
            .map(str::to_owned)
    };

    let mut gb = gb_io::seq::Seq::empty();

    gb.name = Some(if seq.name.trim().is_empty() {
        "unnamed".to_owned()
    } else {
        seq.name.trim().to_owned()
    });

    gb.topology = match seq.topology {
        Some(SeqTopology::Circular) => Topology::Circular,
        _ => Topology::Linear,
    };

    gb.molecule_type = match seq.seq_type() {
        SeqType::AminoAcid => None,
        SeqType::Dna => Some(meta(KEY_MOL_TYPE).unwrap_or_else(|| "DNA".to_owned())),
        SeqType::Rna => Some(meta(KEY_MOL_TYPE).unwrap_or_else(|| "RNA".to_owned())),
    };

    gb.division = meta(KEY_DIVISION).unwrap_or_else(|| DIVISION_DEFAULT.to_owned());
    gb.date = meta(KEY_DATE)
        .and_then(|d| parse_date(&d))
        .and_then(|(y, m, d)| Date::from_ymd(y, m as u32, d as u32).ok())
        .or_else(|| Some(date_today()));

    gb.definition = Some(meta(SEQ_DESCRIPTION_KEY).unwrap_or_else(|| ".".to_owned()));
    gb.accession = meta(KEY_ACCESSION);
    gb.version = meta(KEY_VERSION);
    gb.dblink = meta(KEY_DBLINK);
    gb.keywords = Some(meta(KEY_KEYWORDS).unwrap_or_else(|| ".".to_owned()));

    gb.source = Some(Source {
        source: meta(KEY_SOURCE).unwrap_or_else(|| ".".to_owned()),
        organism: Some(meta(KEY_ORGANISM).unwrap_or_else(|| ".".to_owned())),
    });

    gb.references = record.references.iter().map(Into::into).collect();
    gb.comments = record.comments.clone();

    gb.seq = seq.data.to_letters().to_ascii_lowercase().into_bytes();
    gb.len = Some(gb.seq.len());

    gb.features = seq
        .features
        .iter()
        .filter_map(|f| feature_to_gb(f, gb.seq.len()))
        .collect();

    gb
}

fn feature_to_gb(feature: &SeqFeature, seq_len: usize) -> Option<Feature> {
    let range_loc = |start: usize, end: usize| {
        Location::Range(
            (start.saturating_sub(1) as i64, Before(false)),
            (end as i64, After(false)),
        )
    };

    let mut locs = Vec::new();
    for r in &feature.ranges {
        if r.end >= r.start {
            locs.push(range_loc(r.start, r.end));
        } else {
            // Wraps past the origin of a circular sequence.
            locs.push(range_loc(r.start, seq_len));
            locs.push(range_loc(1, r.end));
        }
    }

    let mut location = match locs.len() {
        0 => return None,
        1 => locs.remove(0),
        _ => Location::Join(locs),
    };

    if feature.strand == Strand::Reverse {
        location = Location::Complement(Box::new(location));
    }

    let mut qualifiers = Vec::new();

    if !feature.label.is_empty() {
        qualifiers.push((Cow::Borrowed(QUAL_LABEL), Some(feature.label.clone())));
    }

    for (key, val) in &feature.qualifiers {
        // A qualifier with no value, e.g. `/pseudo`, is stored with an empty one.
        let val = (!val.is_empty()).then(|| val.clone());
        qualifiers.push((Cow::Owned(key.clone()), val));
    }

    // A location that isn't `complement(...)` doesn't say whether it's on the forward strand, or
    // not stranded, so this does. For the reverse strand, the location says so already.
    if feature.strand == Strand::Forward {
        qualifiers.push((Cow::Borrowed(QUAL_DIRECTION), Some("RIGHT".to_owned())));
    }

    Some(Feature {
        kind: Cow::Owned(feature.kind.clone()),
        location,
        qualifiers,
    })
}

/// Parse a date in GenBank's format, e.g. "21-JUN-1999", as stored under [`KEY_DATE`]. Returns
/// (year, month, day).
pub fn parse_date(text: &str) -> Option<(i32, u8, u8)> {
    let mut parts = text.trim().split('-');

    let day: u8 = parts.next()?.parse().ok()?;
    let month = parts.next()?.to_ascii_uppercase();
    let month = MONTHS.iter().position(|m| *m == month)? as u8 + 1;
    let year = parts.next()?.parse().ok()?;

    (1..=31).contains(&day).then_some((year, month, day))
}

/// Format a date in GenBank's format, e.g. "21-JUN-1999", as stored under [`KEY_DATE`].
pub fn format_date(year: i32, month: u8, day: u8) -> String {
    let month = MONTHS
        .get((month as usize).saturating_sub(1))
        .copied()
        .unwrap_or("JAN");

    format!("{day:02}-{month}-{year:04}")
}

/// Today's date (UTC).
fn date_today() -> Date {
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

    Date::from_ymd(year as i32, month as u32, day as u32)
        .unwrap_or_else(|_| Date::from_ymd(1970, 1, 1).unwrap())
}
