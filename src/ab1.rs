//! For reading AB1 trace files. (Applied Biosystem's sequencing)
//! [BioPython docs](https://biopython.org/wiki/ABI_traces)
//!
//! Adapted directly from this [BioPython code](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py)
//!
//! We are unable to find the official format spec for AB1 files.

use std::{
    collections::HashMap,
    fs::File,
    io::{self, ErrorKind, Read, Seek, SeekFrom},
    path::Path,
};

#[cfg(feature = "encode")]
use bincode::{Decode, Encode};
use bio::io::fastq;
use na_seq::{Seq, seq_from_str};

const HEADER_SIZE: usize = 26;
const DIR_SIZE: usize = 28;

/// The data structure representing AB1 data.
#[cfg_attr(feature = "encode", derive(Encode, Decode))]
#[derive(Clone, Debug, Default)]
pub struct SeqRecordAb1 {
    pub id: String,
    pub name: String,
    pub description: String,
    pub sequence: Seq,
    pub sequence_user: Option<Seq>,
    pub annotations: HashMap<String, String>,
    pub quality: Option<Vec<u8>>,
    pub quality_user: Option<Vec<u8>>,
    pub peak_heights: Vec<u16>,
    /// Analyzed data, for each channel.
    /// G
    pub data_ch1: Vec<u16>,
    /// A
    pub data_ch2: Vec<u16>,
    /// T
    pub data_ch3: Vec<u16>,
    /// C
    pub data_ch4: Vec<u16>,
    /// Peak locations.
    pub peak_locations: Vec<u16>,
    /// Peak locations edited by user.
    pub peak_locations_user: Option<Vec<u16>>,
}

#[derive(Debug)]
struct Header {
    pub file_version: u16,
    pub tag_name: Seq, // 4 bytes always.
    pub tag_number: u32,
    pub element_type_code: u16,
    pub element_size: u16,
    pub num_elements: usize,
    pub data_size: u32,
    pub data_offset: u32,
}

impl Header {
    pub fn from_bytes(bytes: [u8; HEADER_SIZE]) -> io::Result<Self> {
        let seq_str = std::str::from_utf8(&bytes[2..6]).unwrap().to_owned(); // todo: Handle
        Ok(Self {
            file_version: u16::from_be_bytes(bytes[0..2].try_into().unwrap()),
            tag_name: seq_from_str(&seq_str),
            tag_number: u32::from_be_bytes(bytes[6..10].try_into().unwrap()),
            element_type_code: u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
            element_size: u16::from_be_bytes(bytes[12..14].try_into().unwrap()),
            num_elements: u32::from_be_bytes(bytes[14..18].try_into().unwrap()) as usize,
            data_size: u32::from_be_bytes(bytes[18..22].try_into().unwrap()),
            data_offset: u32::from_be_bytes(bytes[22..26].try_into().unwrap()),
        })
    }
}

#[derive(Debug)]
struct Dir {
    // todo: This breakdown into fields is wrong, but I'm not sure how.
    pub tag_name: String, // 4 bytes
    pub tag_number: u32,
    pub elem_code: u16,
    // pub a: u16, // placeholder
    pub num_elements: usize,
    pub data_size: usize,
    pub data_offset: usize,
    // pub b: u32, // placeholder
    pub tag_offset: usize,
    // todo: Tag offset??
}

impl Dir {
    pub fn from_bytes(bytes: [u8; DIR_SIZE], tag_offset: usize) -> io::Result<Self> {
        Ok(Self {
            tag_name: std::str::from_utf8(&bytes[..4]).unwrap().to_owned(), // todo: Handle
            tag_number: u32::from_be_bytes(bytes[4..8].try_into().unwrap()),
            elem_code: u16::from_be_bytes(bytes[8..10].try_into().unwrap()),
            // a: u16::from_be_bytes(bytes[10..12].try_into().unwrap()),
            num_elements: u32::from_be_bytes(bytes[12..16].try_into().unwrap()) as usize,
            data_size: u32::from_be_bytes(bytes[16..20].try_into().unwrap()) as usize,
            data_offset: u32::from_be_bytes(bytes[20..24].try_into().unwrap()) as usize,
            // b: u32::from_be_bytes(bytes[24..28].try_into().unwrap()),
            tag_offset,
        })
    }
}

#[derive(Debug)]
struct AbiIterator<R: Read + Seek> {
    stream: R,
}

impl<R: Read + Seek> AbiIterator<R> {
    pub fn new(mut stream: R) -> io::Result<Self> {
        let mut marker = [0; 4];
        stream.read_exact(&mut marker)?;
        if &marker != b"ABIF" {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Invalid AB1 file start marker",
            ));
        }
        Ok(Self { stream })
    }

    pub fn next(&mut self) -> io::Result<Option<SeqRecordAb1>> {
        let mut result = SeqRecordAb1::default();
        let mut header_data = [0; HEADER_SIZE];

        if self.stream.read(&mut header_data)? == 0 {
            return Ok(None); // EOF
        }

        let header = Header::from_bytes(header_data)?;

        for i in 0..header.num_elements {
            // todo: QC data_offset; coming out much too high.
            // Note: Element size should always be DIR_SIZE.
            let start = header.data_offset as usize + i * header.element_size as usize;

            self.stream.seek(SeekFrom::Start(start as u64))?;
            let mut dir_buf = [0; DIR_SIZE];
            if self.stream.read(&mut dir_buf)? == 0 {
                return Ok(None); // EOF
            };

            let mut dir = Dir::from_bytes(dir_buf, start)?;

            let key = format!("{}{}", dir.tag_name, dir.tag_number);
            // println!("DIR: {:?}, KEY: {:?}", dir, key);

            if dir.data_size <= 4 {
                dir.data_offset = dir.tag_offset + 20;
            }

            self.stream.seek(SeekFrom::Start(dir.data_offset as u64))?;
            let mut tag_buf = vec![0; dir.data_size];
            if self.stream.read(&mut tag_buf)? == 0 {
                return Ok(None); // EOF
            };

            let tag_data = parse_tag_data(dir.elem_code, dir.num_elements, &tag_buf)?;

            // todo: This section is repetative.
            match key.as_str() {
                "PBAS1" => match tag_data {
                    TagData::Str(s) => {
                        result.sequence_user = Some(seq_from_str(&s));
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid PBAS sequence",
                        ));
                    }
                },
                "PBAS2" => match tag_data {
                    TagData::Str(s) => {
                        result.sequence = seq_from_str(&s);
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid PBAS sequence",
                        ));
                    }
                },
                "PCON1" => {
                    match tag_data {
                        TagData::Str(s) => {
                            // Note: We have reversed the above conversion from bytes; we get
                            // the "char" type for both.
                            result.quality_user = Some(s.as_bytes().to_vec());
                        }
                        _ => {
                            return Err(io::Error::new(
                                ErrorKind::InvalidData,
                                "Invalid quality data",
                            ));
                        }
                    }
                }
                // Quality values
                "PCON2" => {
                    match tag_data {
                        TagData::Str(s) => {
                            // Note: We have reversed the above conversion from bytes; we get
                            // the "char" type for both.
                            result.quality = Some(s.as_bytes().to_vec());
                        }
                        _ => {
                            return Err(io::Error::new(
                                ErrorKind::InvalidData,
                                "Invalid quality data",
                            ));
                        }
                    }
                }
                // Sample ID
                "SMPL1" => match tag_data {
                    TagData::Str(s) => result.id = s,
                    _ => return Err(io::Error::new(ErrorKind::InvalidData, "Invalid sample ID")),
                },
                "PLOC1" => match tag_data {
                    TagData::U16(d) => {
                        result.peak_locations_user = Some(d);
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid peak location data",
                        ));
                    }
                },
                "PLOC2" => match tag_data {
                    TagData::U16(d) => {
                        result.peak_locations = d;
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid peak location data",
                        ));
                    }
                },
                "DATA9" => match tag_data {
                    TagData::U16(d) => {
                        result.data_ch1 = d;
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid height data",
                        ));
                    }
                },
                "DATA10" => match tag_data {
                    TagData::U16(d) => {
                        result.data_ch2 = d;
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid height data",
                        ));
                    }
                },
                "DATA11" => match tag_data {
                    TagData::U16(d) => {
                        result.data_ch3 = d;
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid height data",
                        ));
                    }
                },
                "DATA12" => match tag_data {
                    TagData::U16(d) => {
                        result.data_ch4 = d;
                    }
                    _ => {
                        return Err(io::Error::new(
                            ErrorKind::InvalidData,
                            "Invalid height data",
                        ));
                    }
                },
                _ => {
                    // todo: Implement others A/R.
                    eprintln!("Invalid key in AB1 file: {key:?}");
                    eprintln!("Tag data for this key: {tag_data:?}");
                }
            }
        }

        Ok(Some(result))
    }
}

// Helper function to parse ABI tags
// fn parse_abi_tag(data: &[u8]) -> Result<(String, String), Box<dyn Error>> {
fn parse_abi_tag(data: &[u8]) -> io::Result<(String, String)> {
    let tag_name = String::from_utf8_lossy(&data[0..4]).to_string();
    let tag_number = u32::from_be_bytes(
        data[4..8]
            .try_into()
            .map_err(|err| io::Error::new(ErrorKind::InvalidData, err))?,
    );
    Ok((tag_name, tag_number.to_string()))
}

fn abi_trim(seq_record: &fastq::Record) -> fastq::Record {
    // Richard Mott's modified trimming algorithm.

    let segment = 20; // Minimum sequence length
    let cutoff = 0.05; // Default cutoff value for calculating base score

    // If the length of the sequence is less than or equal to the segment size, return as is.
    if seq_record.seq().len() <= segment {
        return seq_record.clone();
    }

    // Calculate base scores from quality values.
    let score_list: Vec<f64> = seq_record
        .qual()
        .iter()
        .map(|&qual| cutoff - 10f64.powf((qual as f64) / -10.0))
        .collect();

    // Calculate cumulative score, initialize with zero.
    let mut cumulative_scores: Vec<f64> = vec![0.0];
    let mut trim_start = 0;
    let mut start_flag = false;

    for i in 1..score_list.len() {
        let score = cumulative_scores[i - 1] + score_list[i];
        if score < 0.0 {
            cumulative_scores.push(0.0);
        } else {
            cumulative_scores.push(score);
            if !start_flag {
                // Set trim_start when cumulative score is first greater than zero
                trim_start = i;
                start_flag = true;
            }
        }
    }

    // Find the index of the highest cumulative score to mark the end of the trimming segment
    let trim_finish = cumulative_scores
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);

    // Extract the trimmed sequence
    let trimmed_seq = &seq_record.seq()[trim_start..trim_finish];
    let trimmed_qual = &seq_record.qual()[trim_start..trim_finish];

    // Create a new trimmed record and return it
    fastq::Record::with_attrs(seq_record.id(), None, trimmed_seq, trimmed_qual)
}

#[derive(Debug)]
enum TagData {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    Str(String),
}

fn parse_tag_data(elem_code: u16, _elem_num: usize, data: &[u8]) -> io::Result<TagData> {
    //     1: "b",  # byte
    //     2: "s",  # char
    //     3: "H",  # word
    //     4: "h",  # short
    //     5: "i",  # long
    //     6: "2i",  # rational, legacy unsupported
    //     7: "f",  # float
    //     8: "d",  # double
    //     10: "h2B",  # date
    //     11: "4B",  # time
    //     12: "2i2b",  # thumb
    //     13: "B",  # bool
    //     14: "2h",  # point, legacy unsupported
    //     15: "4h",  # rect, legacy unsupported
    //     16: "2i",  # vPoint, legacy unsupported
    //     17: "4i",  # vRect, legacy unsupported
    //     18: "s",  # pString
    //     19: "s",  # cString
    //     20: "2i",  # tag, legacy unsupported

    match elem_code {
        // 2 => Some(TagData::U8(data.to_vec())),
        2 => Ok(TagData::Str(
            std::str::from_utf8(data).unwrap_or("").to_string(),
        )),
        4 => {
            let as_u16 = data
                .chunks_exact(2)
                .map(|chunk| u16::from_be_bytes([chunk[0], chunk[1]]))
                .collect();
            Ok(TagData::U16(as_u16))
        }
        5 => {
            let as_u32 = data
                .chunks_exact(4)
                .map(|chunk| u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                .collect();
            Ok(TagData::U32(as_u32))
        }

        _ => {
            // todo: Handle appropriately.
            Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid data type in AB1 file: {elem_code}"),
            ))
        }
    }
}

fn read_string<R: Read>(reader: &mut R, length: usize) -> io::Result<String> {
    let mut buffer = vec![0; length];
    reader.read_exact(&mut buffer)?;
    Ok(String::from_utf8_lossy(&buffer)
        .trim_end_matches(char::from(0))
        .to_string())
}

/// Read a file in the GenBank format.
/// [Rust docs ref of fields](https://docs.rs/gb-io/latest/gb_io/seq/struct.Seq.html)
pub fn import_ab1(path: &Path) -> io::Result<Vec<SeqRecordAb1>> {
    let file = File::open(path)?;
    let mut iterator = AbiIterator::new(file)?;

    let mut results = Vec::new();

    while let Some(record) = iterator.next()? {
        // println!("{:?}", record);
        results.push(record);
    }
    Ok(results)
}
