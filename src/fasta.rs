//! For reading and writing FASTA files: one or more sequence records, each a `>` header line
//! followed by residue letters.
//!
//! The header is `>ID description`; we store the ID as the sequence's name, and the description
//! under [`SEQ_DESCRIPTION_KEY`]. FASTA can't hold other metadata, so it's not written.
//!
//! FASTA doesn't state whether a record is protein, DNA, or RNA. We infer it from the letters,
//! unless given explicitly. Text with no header line is read as one unnamed record, so plain
//! sequence files load as well.

use std::{
    fs,
    fs::File,
    io::{self, BufWriter, ErrorKind, Write},
    path::Path,
};

use na_seq::{SEQ_DESCRIPTION_KEY, SeqType, Sequence, SequenceData};

/// Residues per line when writing.
const LINE_LEN: usize = 60;

#[derive(Clone, Debug, Default)]
pub struct Fasta {
    pub records: Vec<Sequence>,
    /// When loading: the number of residue letters left out across all records, as they can't
    /// be represented. E.g. ambiguity codes like N or X.
    pub skipped_residues: usize,
}

impl Fasta {
    /// `seq_type` sets the type of every record; if `None`, it's inferred per record.
    pub fn new(text: &str, seq_type: Option<SeqType>) -> io::Result<Self> {
        // (header, residue letters)
        let mut raw: Vec<(Option<&str>, String)> = Vec::new();

        for line in text.lines() {
            let line = line.trim_end();

            if let Some(header) = line.strip_prefix('>') {
                raw.push((Some(header.trim()), String::new()));
                continue;
            }

            // Comment lines, from the original FASTA format.
            if line.starts_with(';') || line.trim().is_empty() {
                continue;
            }

            match raw.last_mut() {
                Some((_, residues)) => residues.push_str(line),
                // Residues before any header: a bare sequence.
                None => raw.push((None, line.to_owned())),
            }
        }

        if raw.is_empty() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "No sequences found in the FASTA data",
            ));
        }

        let mut records = Vec::with_capacity(raw.len());
        let mut skipped_residues = 0;

        for (header, residues) in raw {
            let type_ = seq_type.unwrap_or_else(|| SeqType::infer(&residues));
            let (data, skipped) = SequenceData::from_letters(&residues, type_);
            skipped_residues += skipped;

            let (name, descrip) = match header {
                Some(h) => match h.split_once(char::is_whitespace) {
                    Some((id, descrip)) => (id.to_owned(), descrip.trim()),
                    None => (h.to_owned(), ""),
                },
                None => (String::new(), ""),
            };

            let mut record = Sequence::new(data, name);
            if !descrip.is_empty() {
                record
                    .metadata
                    .insert(SEQ_DESCRIPTION_KEY.to_owned(), descrip.to_owned());
            }

            records.push(record);
        }

        Ok(Self {
            records,
            skipped_residues,
        })
    }

    /// Infers each record's type, except for extensions that are protein-specific, e.g. `faa`.
    pub fn load(path: &Path) -> io::Result<Self> {
        let text = fs::read_to_string(path)?;

        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();

        let seq_type = match ext.as_str() {
            "faa" | "mpfa" => Some(SeqType::AminoAcid),
            _ => None,
        };

        let mut result = Self::new(&text, seq_type)?;
        for record in &mut result.records {
            record.path = Some(path.to_owned());
        }

        Ok(result)
    }

    pub fn to_text(&self) -> String {
        let mut result = String::new();

        for record in &self.records {
            // The ID ends at the first whitespace, so keep it to one token for a clean round trip.
            let id: String = record
                .name
                .trim()
                .chars()
                .map(|c| if c.is_whitespace() { '_' } else { c })
                .collect();

            result.push('>');
            result.push_str(&id);

            if let Some(descrip) = record.description() {
                result.push(' ');
                // A description can't span lines in FASTA.
                result.push_str(&descrip.replace(['\r', '\n'], " "));
            }
            result.push('\n');

            let letters = record.data.to_letters();
            // Residue letters are ASCII, so these byte chunks are valid UTF-8.
            for chunk in letters.as_bytes().chunks(LINE_LEN) {
                result.push_str(std::str::from_utf8(chunk).unwrap_or_default());
                result.push('\n');
            }
        }

        result
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let mut file = BufWriter::new(File::create(path)?);
        file.write_all(self.to_text().as_bytes())?;
        file.flush()
    }
}
