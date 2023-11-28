pub mod fasta_tools {
    use std::fs;
    use std::path::Path;
    use std::fs::OpenOptions;

    use bio::io::fasta;

    use crate::logging::logging;

    pub fn get_filename(sequence_length: usize) -> String {
        let min_length = (sequence_length/1000)*1000;
        let max_length = min_length + 1000;
    
        // pad with leading zeros so that we can concatenate the files together in order with cat
        let min_length_padded = format!("{:0>12}", min_length);
        let max_length_padded = format!("{:0>12}", max_length);
        return format!("sequences_{}-{}.fa", min_length_padded, max_length_padded)
    }

    pub fn break_up_fasta_by_sequence_length<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display>(
        input_fasta_path: P,
        output_directory: P,
        total_sequence_count: usize,
     ) {
        let mut current_count = 0;
        fs::create_dir_all(&output_directory).expect("Error creating output directory");
        let reader = fasta::Reader::from_file(input_fasta_path).expect("Error reading fasta file"); 
        reader.records().enumerate().for_each(|(_i, record)| {
            let record = record.expect("Error reading record");
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
            current_count += 1;
            if current_count % 100000 == 0 {
                let processed_percentage = (current_count as f64 / total_sequence_count as f64) * 100.0;
                log::info!("{} of sequences processed", processed_percentage);
            }
        });

    }




}