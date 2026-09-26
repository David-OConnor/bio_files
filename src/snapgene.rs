//! For reading and writing SnapGene DNA files (`.dna`): sequence, topology, features, primers,
//! and notes.
//!
//! [Unofficial file format description](https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf)
//!
//! Files are divided into packets. Each has a byte indicating its type, a big-endian `u32` payload
//! length, then the payload. The DNA packet is binary; features, primers, and notes are XML.

use std::{
    collections::HashMap,
    fs,
    io::{self, ErrorKind},
    path::Path,
    str,
};

use na_seq::{
    SEQ_DESCRIPTION_KEY, Seq, SeqFeature, SeqRange, SeqTopology, SeqType, Sequence, SequenceData,
    Strand,
};
use quick_xml::{
    de::from_str, escape::resolve_predefined_entity, events::Event, reader::Reader, se::to_string,
};

use crate::snapgene::feature_xml::{
    FeatureSnapGene, Features, PrimerSnapGene, Primers, Qualifier, QualifierValue, Segment,
};

/// "SnapGene", then the sequence type, export version, and import version, as `u16`s.
const COOKIE_PACKET_LEN: usize = 14;
/// In the cookie. We only handle DNA.
const SEQ_TYPE_DNA: u16 = 1;

/// Bit 0 of the DNA packet's flags byte.
const FLAG_CIRCULAR: u8 = 0x01;

/// SnapGene inserts these into feature qualifier and note values; we remove them.
const HTML_TAGS: [&str; 8] = [
    "<html>", "</html>", "<body>", "</body>", "<i>", "</i>", "<b>", "</b>",
];

/// The packets we read or write. There are others, e.g. 0x3, 0x8, 0xd, 0xe, 0x11, and 0x1c.
#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(u8)]
enum PacketType {
    /// The cookie is always the first packet.
    Cookie = 0x09,
    Dna = 0x00,
    Primers = 0x05,
    Notes = 0x06,
    Features = 0x0a,
}

impl PacketType {
    fn from_u8(v: u8) -> Option<Self> {
        Some(match v {
            0x09 => Self::Cookie,
            0x00 => Self::Dna,
            0x05 => Self::Primers,
            0x06 => Self::Notes,
            0x0a => Self::Features,
            _ => return None,
        })
    }
}

/// A primer stored with the sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct SnapGenePrimer {
    pub sequence: Seq,
    pub name: String,
    pub description: Option<String>,
}

#[derive(Clone, Debug)]
pub struct SnapGene {
    /// DNA, with its topology and features. Notes are in its metadata; the description, if any,
    /// under [`SEQ_DESCRIPTION_KEY`].
    pub seq: Sequence,
    pub primers: Vec<SnapGenePrimer>,
    /// When loading: the number of residue letters left out, as they can't be represented. E.g.
    /// ambiguity codes like N. Feature ranges are adjusted to match.
    pub skipped_residues: usize,
}

impl SnapGene {
    pub fn new(buf: &[u8]) -> io::Result<Self> {
        let mut dna = None;
        let mut features = Vec::new();
        let mut primers = Vec::new();
        let mut notes = HashMap::new();

        let mut i = 0;
        let mut first = true;

        while i + 5 <= buf.len() {
            let packet_type = buf[i];
            let payload_len = u32::from_be_bytes(buf[i + 1..i + 5].try_into().unwrap()) as usize;
            i += 5;

            if i + payload_len > buf.len() {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "SnapGene packet length exceeds the file length",
                ));
            }

            let payload = &buf[i..i + payload_len];
            i += payload_len;

            let packet_type = PacketType::from_u8(packet_type);

            if first {
                first = false;
                check_cookie(packet_type, payload)?;
                continue;
            }

            // Packets other than the DNA are optional; if one can't be read, the sequence is still
            // usable.
            match packet_type {
                Some(PacketType::Dna) => {
                    if payload.is_empty() {
                        return Err(io::Error::new(ErrorKind::InvalidData, "Empty DNA packet"));
                    }
                    // Flags, then the sequence.
                    dna = Some((payload[0], &payload[1..]));
                }
                Some(PacketType::Features) => match parse_features(payload) {
                    Ok(v) => features = v,
                    Err(e) => eprintln!("Error parsing SnapGene features: {e}"),
                },
                Some(PacketType::Primers) => match parse_primers(payload) {
                    Ok(v) => primers = v,
                    Err(e) => eprintln!("Error parsing SnapGene primers: {e}"),
                },
                Some(PacketType::Notes) => match parse_notes(payload) {
                    Ok(v) => notes = v,
                    Err(e) => eprintln!("Error parsing SnapGene notes: {e}"),
                },
                _ => (),
            }
        }

        if first {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Empty SnapGene file",
            ));
        }

        let Some((flags, letters)) = dna else {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "No DNA packet in the SnapGene file",
            ));
        };

        let letters = str::from_utf8(letters).map_err(|e| {
            io::Error::new(ErrorKind::InvalidData, format!("Invalid DNA packet: {e}"))
        })?;

        let (mut seq, skipped_residues) =
            Sequence::from_letters_with_features(letters, SeqType::Dna, String::new(), features);

        seq.topology = Some(if flags & FLAG_CIRCULAR != 0 {
            SeqTopology::Circular
        } else {
            SeqTopology::Linear
        });
        seq.metadata = notes;

        Ok(Self {
            seq,
            primers,
            skipped_residues,
        })
    }

    /// The sequence is named after the file.
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut result = Self::new(&fs::read(path)?)?;

        result.seq.name = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        result.seq.path = Some(path.to_owned());

        Ok(result)
    }

    /// Writes the sequence, topology, features, and primers. Notes aren't written.
    pub fn to_bytes(&self) -> io::Result<Vec<u8>> {
        if self.seq.seq_type() != SeqType::Dna {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "Only DNA can be saved as a SnapGene DNA file",
            ));
        }

        let mut buf = Vec::new();

        let mut cookie = Vec::with_capacity(COOKIE_PACKET_LEN);
        cookie.extend(b"SnapGene");
        cookie.extend(SEQ_TYPE_DNA.to_be_bytes());
        // Export and import versions.
        cookie.extend([0; 4]);
        push_packet(&mut buf, PacketType::Cookie, &cookie);

        let flags = match self.seq.topology {
            Some(SeqTopology::Circular) => FLAG_CIRCULAR,
            _ => 0,
        };

        let mut dna = vec![flags];
        dna.extend(self.seq.data.to_letters().to_ascii_lowercase().bytes());
        push_packet(&mut buf, PacketType::Dna, &dna);

        push_packet(
            &mut buf,
            PacketType::Features,
            features_xml(&self.seq.features)?.as_bytes(),
        );
        push_packet(
            &mut buf,
            PacketType::Primers,
            primers_xml(&self.primers)?.as_bytes(),
        );

        Ok(buf)
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        fs::write(path, self.to_bytes()?)
    }
}

fn check_cookie(packet_type: Option<PacketType>, payload: &[u8]) -> io::Result<()> {
    if packet_type != Some(PacketType::Cookie)
        || payload.len() < COOKIE_PACKET_LEN
        || &payload[..8] != b"SnapGene"
    {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Not a SnapGene file: missing its header",
        ));
    }

    let seq_type = u16::from_be_bytes([payload[8], payload[9]]);
    if seq_type != SEQ_TYPE_DNA {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            format!(
                "Unsupported SnapGene sequence type: {seq_type}. Only DNA files are supported."
            ),
        ));
    }

    Ok(())
}

fn push_packet(buf: &mut Vec<u8>, packet_type: PacketType, payload: &[u8]) {
    buf.push(packet_type as u8);
    buf.extend((payload.len() as u32).to_be_bytes());
    buf.extend(payload);
}

fn payload_str(payload: &[u8]) -> io::Result<&str> {
    str::from_utf8(payload).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert payload to string: {e}"),
        )
    })
}

fn strip_html(text: &str) -> String {
    let mut result = text.to_owned();
    for tag in HTML_TAGS {
        result = result.replace(tag, "");
    }
    result
}

mod feature_xml {
    use std::str::FromStr;

    use serde::{Deserialize, Deserializer, Serialize};

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Features {
        #[serde(rename = "Feature", default)]
        pub inner: Vec<FeatureSnapGene>,
    }

    // Workaround for parsing "" into None not supported natively by quick-xml/serde.
    fn deserialize_directionality<'de, D>(deserializer: D) -> Result<Option<u8>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: Option<String> = Option::deserialize(deserializer)?;
        match s {
            Some(s) if !s.is_empty() => {
                u8::from_str(&s).map(Some).map_err(serde::de::Error::custom)
            }
            _ => Ok(None),
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct FeatureSnapGene {
        #[serde(rename = "@type")]
        pub feature_type: Option<String>,
        /// 1: forward. 2: reverse.
        #[serde(
            rename = "@directionality",
            deserialize_with = "deserialize_directionality",
            default,
            skip_serializing_if = "Option::is_none"
        )]
        pub directionality: Option<u8>,
        #[serde(rename = "@name", default)]
        pub name: Option<String>,
        // Other Feature attributes: allowSegmentOverlaps (0/1), consecutiveTranslationNumbering (0/1)
        #[serde(rename = "Segment", default)]
        pub segments: Vec<Segment>,
        #[serde(rename = "Q", default)]
        pub qualifiers: Vec<Qualifier>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Segment {
        #[serde(rename = "@type", default, skip_serializing_if = "Option::is_none")]
        pub segment_type: Option<String>,
        /// E.g. "1-120", 1-based and inclusive.
        #[serde(rename = "@range", default)]
        pub range: Option<String>,
        #[serde(rename = "@name", default, skip_serializing_if = "Option::is_none")]
        pub name: Option<String>,
        /// Hex.
        #[serde(rename = "@color", default, skip_serializing_if = "Option::is_none")]
        pub color: Option<String>,
        // Other fields: "translated": 0/1
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Qualifier {
        #[serde(rename = "@name")]
        pub name: String,
        #[serde(rename = "V", default)]
        pub values: Vec<QualifierValue>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct QualifierValue {
        #[serde(rename = "@text", default, skip_serializing_if = "Option::is_none")]
        pub text: Option<String>,
        #[serde(rename = "@predef", default, skip_serializing_if = "Option::is_none")]
        pub predef: Option<String>,
        #[serde(rename = "@int", default, skip_serializing_if = "Option::is_none")]
        pub int: Option<i32>,
    }

    #[derive(Debug, Serialize, Deserialize)]
    pub struct Primers {
        #[serde(rename = "Primer", default)]
        pub inner: Vec<PrimerSnapGene>,
    }

    // Note; We have left out the binding site and a number of other fields. This includes melting
    // temperature, which can be calculated.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct PrimerSnapGene {
        #[serde(rename = "@sequence")]
        pub sequence: String,
        #[serde(rename = "@name")]
        pub name: String,
        #[serde(rename = "@description", default)]
        pub description: String,
    }
}

fn parse_features(payload: &[u8]) -> io::Result<Vec<SeqFeature>> {
    let features: Features = from_str(payload_str(payload)?).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse features: {e}"),
        )
    })?;

    let mut result = Vec::with_capacity(features.inner.len());

    for feature_sg in features.inner {
        let strand = match feature_sg.directionality {
            Some(1) => Strand::Forward,
            Some(2) => Strand::Reverse,
            _ => Strand::None,
        };

        let mut qualifiers = Vec::new();
        for qual in &feature_sg.qualifiers {
            // Generally one value per qualifier, which has one of text, int, and predef. If there
            // are several, they're stored as separate qualifiers.
            for val in &qual.values {
                let v = val
                    .text
                    .clone()
                    .or_else(|| val.predef.clone())
                    .or_else(|| val.int.map(|v| v.to_string()))
                    .unwrap_or_default();

                qualifiers.push((qual.name.clone(), strip_html(&v)));
            }
        }

        let ranges = feature_sg
            .segments
            .iter()
            .filter_map(|s| s.range.as_deref().and_then(range_from_str))
            .collect();

        // One color per feature; SnapGene has one per segment.
        let color = feature_sg.segments.iter().find_map(|s| s.color.clone());

        result.push(SeqFeature {
            kind: feature_sg.feature_type.unwrap_or_default(),
            label: feature_sg.name.unwrap_or_default(),
            ranges,
            strand,
            color,
            qualifiers,
        });
    }

    Ok(result)
}

fn parse_primers(payload: &[u8]) -> io::Result<Vec<SnapGenePrimer>> {
    let primers: Primers = from_str(payload_str(payload)?).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to parse primers: {e}"),
        )
    })?;

    let result = primers
        .inner
        .into_iter()
        .map(|p| {
            let sequence = match SequenceData::from_letters(&p.sequence, SeqType::Dna).0 {
                SequenceData::Dna(s) => s,
                _ => unreachable!(),
            };

            SnapGenePrimer {
                sequence,
                name: p.name,
                description: (!p.description.is_empty()).then_some(p.description),
            }
        })
        .collect();

    Ok(result)
}

/// The notes packet is a `<Notes>` element with a child element per field, e.g. `<Type>`,
/// `<CreatedBy>`, or `<Description>`. We read each field with text content into metadata. The
/// description goes under [`SEQ_DESCRIPTION_KEY`].
fn parse_notes(payload: &[u8]) -> io::Result<HashMap<String, String>> {
    // Not trimming text: entities split it into several events, and trimming each would remove
    // the spaces around them. Values are trimmed once complete.
    let mut reader = Reader::from_str(payload_str(payload)?);

    let mut result = HashMap::new();
    let mut depth = 0;
    // The field currently being read, and its text.
    let mut field: Option<(String, String)> = None;

    loop {
        let event = reader.read_event().map_err(notes_err)?;

        match event {
            Event::Start(e) => {
                depth += 1;
                if depth == 2 {
                    let name = String::from_utf8_lossy(e.name().as_ref()).into_owned();
                    field = Some((name, String::new()));
                }
            }
            Event::Text(t) if depth == 2 => {
                if let Some((_, text)) = &mut field {
                    text.push_str(&t.decode().map_err(notes_err)?);
                }
            }
            // Entities, e.g. `&lt;`, and character references are separate from the text around
            // them.
            Event::GeneralRef(r) if depth == 2 => {
                if let Some((_, text)) = &mut field {
                    if let Some(c) = r.resolve_char_ref().map_err(notes_err)? {
                        text.push(c);
                    } else if let Some(v) =
                        resolve_predefined_entity(&r.decode().map_err(notes_err)?)
                    {
                        text.push_str(v);
                    }
                }
            }
            Event::End(_) => {
                if depth == 2
                    && let Some((name, text)) = field.take()
                {
                    let text = strip_html(&text).trim().to_owned();

                    if !text.is_empty() {
                        let key = if name == "Description" {
                            SEQ_DESCRIPTION_KEY.to_owned()
                        } else {
                            name
                        };
                        result.insert(key, text);
                    }
                }
                depth -= 1;
            }
            Event::Eof => break,
            _ => (),
        }
    }

    Ok(result)
}

fn notes_err(e: impl std::fmt::Display) -> io::Error {
    io::Error::new(
        ErrorKind::InvalidData,
        format!("Unable to parse notes: {e}"),
    )
}

/// E.g. "1-120".
fn range_from_str(range: &str) -> Option<SeqRange> {
    let (start, end) = range.split_once('-')?;
    Some(SeqRange::new(
        start.trim().parse().ok()?,
        end.trim().parse().ok()?,
    ))
}

fn features_xml(features: &[SeqFeature]) -> io::Result<String> {
    let mut features_sg = Features { inner: Vec::new() };

    for feature in features {
        let directionality = match feature.strand {
            Strand::Forward => Some(1),
            Strand::Reverse => Some(2),
            Strand::None => None,
        };

        let segments = feature
            .ranges
            .iter()
            .map(|r| Segment {
                segment_type: None,
                range: Some(format!("{}-{}", r.start, r.end)),
                name: None,
                color: feature.color.clone(),
            })
            .collect();

        let qualifiers = feature
            .qualifiers
            .iter()
            .map(|(key, value)| Qualifier {
                name: key.clone(),
                values: vec![QualifierValue {
                    text: Some(value.clone()),
                    predef: None,
                    int: None,
                }],
            })
            .collect();

        features_sg.inner.push(FeatureSnapGene {
            feature_type: Some(feature.kind.clone()),
            directionality,
            name: Some(feature.label.clone()),
            segments,
            qualifiers,
        });
    }

    to_string(&features_sg).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert features to XML: {e}"),
        )
    })
}

fn primers_xml(primers: &[SnapGenePrimer]) -> io::Result<String> {
    let primers_sg = Primers {
        inner: primers
            .iter()
            .map(|p| PrimerSnapGene {
                sequence: SequenceData::Dna(p.sequence.clone())
                    .to_letters()
                    .to_ascii_lowercase(),
                name: p.name.clone(),
                description: p.description.clone().unwrap_or_default(),
            })
            .collect(),
    };

    to_string(&primers_sg).map_err(|e| {
        io::Error::new(
            ErrorKind::InvalidData,
            format!("Unable to convert primers to XML: {e}"),
        )
    })
}
