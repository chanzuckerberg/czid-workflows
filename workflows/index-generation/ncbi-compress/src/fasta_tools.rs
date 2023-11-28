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


    pub fn get_filename(sequence_length: usize) -> String {
        let min_length = (sequence_length/1000)*1000;
        let max_length = min_length + 1000;

        // pad with leading zeros so that we can concatenate the files together in order with cat later
        let min_length_padded = format!("{:0>12}", min_length);
        let max_length_padded = format!("{:0>12}", max_length);
        return format!("sequences_{}-{}.fa", min_length_padded, max_length_padded)
    }

    fn process_seq_chunk<P: std::fmt::Display>(record: &fasta::Record, output_directory: &P) {
        let sequence_length = record.seq().len();
        let filename = get_filename(sequence_length);
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
                process_seq_chunk(rec, &output_directory);
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


// use std::{fs, path::Path};
// use tempfile::tempdir;

// #[test]
// fn test_break_up_fasta_by_sequence_lengthy() {
//     // Setup
//     let temp_dir = tempdir().unwrap();
//     let test_dir_outputs = Path::new("test_data/fasta_tools/outputs");
//     let fasta_file = Path::new("test_data/fasta_tools/inputs/nt");

//     break_up_fasta_by_sequence_length(fasta_file, temp_dir.path(), 1000, 1000);

//     // Compare files
//     for entry in fs::read_dir(test_dir).unwrap() {
//         let entry = entry.unwrap();
//         let test_file_path = entry.path();
//         let test_file_name = test_file_path.file_name().unwrap();
//         let output_file_path = temp_dir.path().join(test_file_name);

//         let test_data = fs::read(test_file_path).unwrap();
//         let output_data = fs::read(output_file_path).unwrap();

//         assert_eq!(test_data, output_data, "Files do not match: {:?}", test_file_name);
//     }
// }
