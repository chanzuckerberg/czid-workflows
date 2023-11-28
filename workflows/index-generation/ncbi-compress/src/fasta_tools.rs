pub mod fasta_tools {
    use std::path::Path;
    use std::fs::OpenOptions;
    use std::{borrow::BorrowMut, fs};

    use bio::io::fasta;
    use rayon::prelude::*;
    // use rayon::iter::IntoParallelRefIterator;

    pub fn get_filename(sequence_length: usize) -> String {
        let min_length = (sequence_length/1000)*1000;
        let max_length = min_length + 1000;
    
        // pad with leading zeros so that we can concatenate the files together in order with cat
        let min_length_padded = format!("{:0>12}", min_length);
        let max_length_padded = format!("{:0>12}", max_length);
        return format!("sequences_{}-{}.fa", min_length_padded, max_length_padded)
    }

    fn process_seq_chunk<P: std::fmt::Display>(record: &fasta::Record, output_directory: &P) -> () {
        //let record = c.unwrap();
        let sequence_length = record.seq().len();
        let filename = get_filename(sequence_length);
        let output_path = format!("{}/{}", output_directory, filename);
        let file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open(output_path)
            .expect("Error opening fasta file");
        let mut writer = fasta::Writer::new(file);
        writer.write_record(&record).expect("Error writing record");
        // current_count += 1;

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
            current_count += chunk.len();
            if current_count % 100000 == 0 {
                let processed_percentage = (current_count as f64 / total_sequence_count as f64) * 100.0;
                log::info!("{} of sequences processed", processed_percentage);
            }
            // refill chunk with new records from iterator
            chunk = records_iter
            .borrow_mut()
            .take(chunk_size)
            .collect::<Vec<_>>();
        }
        log::info!("all sequences processed");
    }
}