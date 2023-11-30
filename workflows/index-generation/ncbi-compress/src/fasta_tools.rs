pub mod fasta_tools {
    use std::path::Path;
    use std::fs::OpenOptions;
    use std::{borrow::BorrowMut, fs};
    use std::collections::HashMap;
    use std::sync::Mutex;


    use bio::io::fasta;
    use rayon::prelude::*;
    use lazy_static::lazy_static;

    // Mutex-protected HashMap for file handles
    lazy_static! {
        static ref FILE_HANDLES: Mutex<HashMap<String, std::fs::File>> = Mutex::new(HashMap::new());
    }

    pub fn bin_number_to_floor(sequence_length: &usize, bin_size: &i64) -> i64 {
        let floored = (*sequence_length as f64 / *bin_size as f64).floor();
        (floored * *bin_size as f64) as i64
    }

    pub fn get_filename(sequence_length: usize, bin_size: &i64) -> String {
        let min_length = bin_number_to_floor(&sequence_length, bin_size);
        let max_length = min_length + bin_size - 1;


        // pad with leading zeros so that we can concatenate the files together in order with cat later
        let min_length_padded = format!("{:0>12}", min_length);
        let max_length_padded = format!("{:0>12}", max_length);
        return format!("sequences_{}-{}.fa", min_length_padded, max_length_padded)
    }

    pub fn process_seq_chunk<P: std::fmt::Display>(record: &fasta::Record, output_directory: &P, bin_size: &i64) {
        let sequence_length = record.seq().len();
        let filename = get_filename(sequence_length, bin_size);
        let output_path = format!("{}/{}", output_directory, filename);

        let mut handles = FILE_HANDLES.lock().unwrap(); // Lock the mutex here
        let file = handles.entry(output_path.clone()).or_insert_with(|| {
            OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(&output_path)
                .expect("Error opening fasta file")
        });

        let mut writer = fasta::Writer::new(file);
        writer.write_record(&record).expect("Error writing record");
    }
    
    pub fn break_up_fasta_by_sequence_length<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + std::marker::Send + std::marker::Sync>(
        input_fasta_path: P,
        output_directory: P,
        total_sequence_count: usize,
        chunk_size: usize,
        bin_size: i64,
     ) {
        let mut current_count = 0;
        fs::create_dir_all(&output_directory).expect("Error creating output directory");
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap(); 
        let mut records_iter = reader.records();

        // create initial chunk of records
        let mut chunk = records_iter
        .borrow_mut()
        .take(chunk_size)
        .collect::<Vec<_>>();

        while chunk.len() > 0 {
            chunk.par_iter().filter_map(|r| {
                let rec = r.as_ref().expect("Error reading record");
                process_seq_chunk(rec, &output_directory, &bin_size);
                Some(())
            }).collect::<Vec<_>>();
                
            // update current count and log progress
            current_count += chunk.len();
            let processed_percentage = (current_count as f64 / total_sequence_count as f64) * 100.0;
            log::info!("{} of sequences processed", processed_percentage);

            // refill chunk with new records from iterator
            chunk = records_iter
            .borrow_mut()
            .take(chunk_size)
            .collect::<Vec<_>>();
        }
        log::info!("all sequences processed");
    }
}


// testing starts here
use std::{fs, path::Path};
use std::cmp::Ordering;

use bio::io::fasta;
use tempfile::tempdir;
// use crate::fasta_tools::fasta_tools;

// Define a struct to hold both ID and sequence
#[derive(Eq, PartialEq)]
struct FastaRecord {
    id: String,
    seq: String,
}

// Implement Ord and PartialOrd for FastaRecord to enable sorting by ID
impl Ord for FastaRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        // order based on id
        self.id.cmp(&other.id)
    }
}
impl PartialOrd for FastaRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn create_fasta_records_from_file<P: AsRef<Path> + std::fmt::Debug>(input_fasta_path: P) -> Vec<FastaRecord> {
    let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
    let mut records_iter = reader.records();
    let mut records: Vec<FastaRecord> = Vec::new();
    while let Some(record) = records_iter.next() {
        let record = record.unwrap();
        records.push(FastaRecord {
            id: record.id().to_string(),
            seq: String::from_utf8(record.seq().to_vec()).expect("Invalid UTF-8 sequence"),
        });
    }
    return records;
}

#[test]
fn test_break_up_fasta_by_sequence_length() {
    // Setup
    let bin_size = 1000;
    let temp_dir = tempdir().unwrap();
    let truth_dir_outputs = "test_data/fasta_tools/truth_outputs/break_up_fasta_by_sequence_length";
    let input_fasta_file = "test_data/fasta_tools/inputs/nt";

    let temp_dir_path_str = temp_dir.path().to_str().unwrap();
    fasta_tools::break_up_fasta_by_sequence_length(input_fasta_file, temp_dir_path_str, 1000, 1000, bin_size);

    //Compare files in from test with truth
    for entry in fs::read_dir(temp_dir_path_str).unwrap() {
        let entry = entry.unwrap();
        let test_file_path = entry.path();

        // get file name so we can compare to the file in the test directory
        let test_file_name = test_file_path.file_name().unwrap().to_str().unwrap();
        let truth_file_path = format!("{}/{}", truth_dir_outputs, test_file_name);

        let mut test_data_records = create_fasta_records_from_file(&test_file_path);
        let mut truth_data_records = create_fasta_records_from_file(&truth_file_path);

        assert_eq!(truth_data_records.sort(), test_data_records.sort(), "Files do not match: {:?}", test_file_name);
    }
}

#[test]
fn test_get_filename() {
    let bin_size = 25;
    let sequence_length = 1763;
    let filename = fasta_tools::get_filename(sequence_length, &bin_size);
    assert_eq!(filename, "sequences_000000001750-000000001774.fa");
}

#[test]
fn test_bin_number_to_floor() {
    let bin_size = 25;

    let sequence_length_short = 1763;
    let floored_short = fasta_tools::bin_number_to_floor(&sequence_length_short, &bin_size);
    assert_eq!(floored_short, 1750);

    let sequence_length_long = 1774;
    let floored_long = fasta_tools::bin_number_to_floor(&sequence_length_long, &bin_size);
    assert_eq!(floored_long, 1750);
}

#[test]
fn test_process_seq_chunk() {
    let bin_size = 25;
    let truth_directory = "test_data/fasta_tools/truth_outputs/process_seq_chunk/";
    let test_directory = tempdir().unwrap();
    let temp_dir_path_str = test_directory.path().to_str().unwrap();


    let record = fasta::Record::with_attrs("id", None, b"ACGT");
    fasta_tools::process_seq_chunk(&record, &temp_dir_path_str, &bin_size);

    let test_file_name = "sequences_000000000000-000000000024.fa";
    let test_file_path = format!("{}/{}", temp_dir_path_str, test_file_name);
    let truth_file_path = format!("{}/{}", truth_directory, test_file_name);
    let mut test_data_records = create_fasta_records_from_file(&test_file_path);
    let mut truth_data_records = create_fasta_records_from_file(&truth_file_path);
    assert_eq!(truth_data_records.sort(), test_data_records.sort(), "Files do not match: {:?}", test_file_name);
}

