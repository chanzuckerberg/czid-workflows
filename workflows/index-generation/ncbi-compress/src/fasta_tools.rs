pub mod fasta_tools {
    use std::collections::HashMap;
    use std::fs;
    use std::io::Write;
    use std::os::unix::prelude::FileExt;

    use rand::seq::SliceRandom;
    use rand::thread_rng;

    use bio::io::fasta;
    use rayon::prelude::*;
    use rayon::slice::ParallelSliceMut;

    struct OffsetWriter<'a> {
        file: &'a mut fs::File,
        global_offset: u64,
        local_offset: u64,
    }

    impl<'a> OffsetWriter<'a> {
        pub fn new(file: &'a mut fs::File, offset: u64) -> Self {
            Self {
                file,
                global_offset: offset,
                local_offset: 0,
            }
        }
    }

    impl Write for OffsetWriter<'_> {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let n = self
                .file
                .write_at(buf, self.global_offset + self.local_offset)?;
            self.local_offset += n as u64;
            Ok(n)
        }

        fn flush(&mut self) -> std::io::Result<()> {
            self.file.flush()
        }
    }

    pub fn count_accessions_by_taxid(
        input_taxid_dir: &str,
        output_tsv_path: &str,
    ) -> std::io::Result<()> {
        // loop through all files in the input_taxid_dir
        for entry in fs::read_dir(input_taxid_dir)? {
            let entry = entry?;
            let path = entry.path();
            // get the taxid from the filename:
            // the filename is the taxid with a .fasta extension
            let taxid = path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .parse::<u32>()
                .unwrap();
            // for each file, get the number of sequences in it
            let records = fasta::Reader::from_file(&path)
                .unwrap() // unwrap left here because this is an anyhow error
                .records();

            // write taxid and number of sequences to output_tsv_path
            let mut writer = fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(&output_tsv_path)?;

            writer.write_all(format!("{}\t{}\n", taxid, records.count()).as_bytes())?;
        }
        Ok(())
    }

    pub fn shuffle_fasta_by_sequence_index(
        input_fasta_path: &str,
        output_fasta_path: &str,
    ) -> std::io::Result<()> {
        let mut n_bytes_by_sequence_index = HashMap::new();
        // pre-allocate scratch space to write to so we know how many bytes a sequence is when
        // written. pre-allocate it so we can clear it and we don't need a fresh allocation for
        // every record.
        let mut scratch = Vec::new();

        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();

        for (index, record) in reader.records().enumerate() {
            let record = record.unwrap();
            fasta::Writer::new(&mut scratch).write_record(&record)?;
            n_bytes_by_sequence_index.insert(index, scratch.len() as u64);
            // Clear the scratch space so we can re-use it for the next record
            scratch.clear();
        }

        let mut shuffled_indexes: Vec<usize> = n_bytes_by_sequence_index.keys().cloned().collect();
        // set seed for testing purposes. This should be a parameter
        let mut rng = thread_rng();
        shuffled_indexes.shuffle(&mut rng);

        // Save memory by creating our offset map in place on top of our n_bytes map
        let mut offset_by_sequence_index = n_bytes_by_sequence_index;

        let mut offset = 0;
        for index in shuffled_indexes {
            // Save how many bytes we will need for it
            // it is safe to unwrap here because sorted_lengths are the keys of this map
            let n_bytes = *offset_by_sequence_index.get(&index).unwrap();
            // Overwrite the length with the offset
            offset_by_sequence_index.insert(index, offset);
            // Add the saved number of bytes to the offset, so the next offset starts at the end of
            // the space that this sequence will occupy
            offset += n_bytes;
        }

        let mut writer = fs::File::create(&output_fasta_path)?;

        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        for (index, record) in reader.records().enumerate() {
            let record = record.unwrap();
            // It is safe to unwrap here because all indexes in the file are in the map from the
            // first pass
            let offset = offset_by_sequence_index.get(&index).unwrap();
            let mut offset_writer = OffsetWriter::new(&mut writer, *offset);

            // This scope lets us create a fasta_writer from borrowing offset_writer then drops
            // fasta_writer so we can get offset_writer back and read it's local offset
            {
                let mut fasta_writer = fasta::Writer::new(&mut offset_writer);
                fasta_writer.write_record(&record)?;
            }
        }
        Ok(())
    }

    pub fn sort_fasta_by_sequence_length(
        input_fasta_path: &str,
        output_fasta_path: &str,
    ) -> std::io::Result<()> {
        let mut n_bytes_by_sequence_length = HashMap::new();
        // pre-allocate scratch space to write to so we know how many bytes a sequence is when
        // written. pre-allocate it so we can clear it and we don't need a fresh allocation for
        // every record.
        let mut scratch = Vec::new();

        let mut records = fasta::Reader::from_file(&input_fasta_path)
            .unwrap() // unwrap left here because this is an anyhow error
            .records();

        while let Some(Ok(record)) = records.next() {
            fasta::Writer::new(&mut scratch).write_record(&record)?;
            let sequence_length = record.seq().len();
            let existing = n_bytes_by_sequence_length
                .get(&sequence_length)
                .unwrap_or(&0);
            n_bytes_by_sequence_length.insert(sequence_length, existing + scratch.len() as u64);
            // Clear the scratch space so we can re-use it for the next record
            scratch.clear()
        }
        let mut sorted_lengths: Vec<usize> = n_bytes_by_sequence_length.keys().cloned().collect();

        // sort lengths in descending order
        // subtracting a usize from the largest possible usize inverts the order
        sorted_lengths.par_sort_unstable_by_key(|n| usize::MAX - n);

        // Save memory by creating our offset map in place on top of our n_bytes map
        let mut offset_by_sequence_length = n_bytes_by_sequence_length;
        let mut offset = 0;
        // For each length in descending order
        for length in sorted_lengths {
            // Save how many bytes we will need for it
            // it is safe to unwrap here because sorted_lengths are the keys of this map
            let n_bytes = *offset_by_sequence_length.get(&length).unwrap();
            // Overwrite the length with the offset
            offset_by_sequence_length.insert(length, offset);
            // Add the saved number of bytes to the offset, so the next offset starts at the end of
            // the space that this sequence will occupy
            offset += n_bytes;

            // for example if our longest sequences needs 100 bytes it will be at offset 0 since
            // they will come first, but our second longest sequences will start at offset 100
            // because the longest sequences will occupy the first 100 bytes
        }

        let records = fasta::Reader::from_file(&input_fasta_path)
            .unwrap() // unwrap left here because this is an anyhow error
            .records();
        let mut writer = fs::File::create(&output_fasta_path)?;
        for maybe_record in records {
            let record = maybe_record?;
            let sequence_length = record.seq().len();
            // It is safe to unwrap here because all lengths in the file are in the map from the
            // first pass
            let offset = offset_by_sequence_length.get(&sequence_length).unwrap();
            let mut offset_writer = OffsetWriter::new(&mut writer, *offset);

            // This scope lets us create a fasta_writer from borrowing offset_writer then drops
            // fasta_writer so we can get offset_writer back and read it's local offset
            {
                let mut fasta_writer = fasta::Writer::new(&mut offset_writer);
                fasta_writer.write_record(&record)?;
            }
            // Add the local offset from the writer to the sequence lengths offset so the next
            // sequence of the same length will be written after this sequence
            offset_by_sequence_length.insert(sequence_length, offset + offset_writer.local_offset);
        }
        Ok(())
    }

    pub fn sort_taxid_dir_by_sequence_length(input_taxid_dir: &str, output_taxid_dir: &str) {
        fs::create_dir_all(&output_taxid_dir).expect("Error creating output directory");

        // Read the directory and collect entries
        let entries: Vec<_> = fs::read_dir(input_taxid_dir)
            .expect("Failed to read directory")
            .filter_map(Result::ok)
            .collect();

        entries.par_iter().for_each(|entry| {
            let path = entry.path();
            let input_fasta_path = path.to_str().unwrap();
            let input_fasta_basename = path.file_name().unwrap().to_str().unwrap();
            let output_fasta_path = format!("{}/{}", output_taxid_dir, input_fasta_basename);
            let _ = sort_fasta_by_sequence_length(input_fasta_path, &output_fasta_path);
        });
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use csv::{ReaderBuilder, Trim};
    use tempfile::tempdir;

    use crate::fasta_tools::fasta_tools;
    use crate::util::util;

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

    #[test]
    fn test_sort_fasta_by_sequence_length() {
        let temp_dir = tempdir().unwrap();
        let temp_file = temp_dir.path().join("sorted.fa");
        let output_truth_fasta_file =
            "test_data/fasta_tools/truth_outputs/sorted_taxids/9771.fasta";
        let input_fasta_file = "test_data/fasta_tools/inputs/unordered_taxids/9771.fasta";

        let temp_file_path_str = temp_file.to_str().unwrap();
        fasta_tools::sort_fasta_by_sequence_length(input_fasta_file, temp_file_path_str).unwrap();

        util::compare_fasta_records_from_files(output_truth_fasta_file, temp_file_path_str);
    }

    #[test]
    fn test_count_accessions_by_taxid() {
        let output_truth_tsv_file =
            "test_data/fasta_tools/truth_outputs/count_accessions_by_taxid/output_counts.tsv";

        let test_truth_tsv_file = tempfile::NamedTempFile::new().unwrap();
        let test_truth_tsv_file_path_str = test_truth_tsv_file.path().to_str().unwrap();

        let _ = fasta_tools::count_accessions_by_taxid(
            "test_data/commands/fasta_compress_from_taxid_dir/inputs",
            test_truth_tsv_file_path_str,
        );

        // read truth tsv file and create vec of records
        let mut truth_rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .trim(Trim::All)
            .from_path(output_truth_tsv_file)
            .unwrap();
        let mut truth_records = Vec::new();
        for result in truth_rdr.records() {
            let record = result.unwrap();
            if let (Some(col1), Some(col2)) = (record.get(0), record.get(1)) {
                truth_records.push((col1.to_string(), col2.to_string()));
            } else {
                println!("Error parsing record: {:?}", record);
                continue;
            }
        }

        // read test tsv file and create vec of records
        let mut test_rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .trim(Trim::All)
            .from_path(test_truth_tsv_file_path_str)
            .unwrap();
        let mut test_records = Vec::new();
        for result in test_rdr.records() {
            let record = result.unwrap();
            if let (Some(col1), Some(col2)) = (record.get(0), record.get(1)) {
                test_records.push((col1.to_string(), col2.to_string()));
            } else {
                println!("Error parsing record: {:?}", record);
                continue;
            }
        }

        // verify that test_records and truth_records contain the same records
        for truth_record in &truth_records {
            assert!(test_records.contains(&truth_record));
        }
        for test_record in &test_records {
            assert!(truth_records.contains(&test_record));
        }
    }

    #[test]
    fn test_sort_taxid_dir_by_sequence_length() {
        let temp_dir = tempdir().unwrap();
        let temp_dir_path_str = temp_dir.path().to_str().unwrap();
        let output_truth_taxid_dir = "test_data/fasta_tools/truth_outputs/sorted_taxids";
        let input_taxid_dir = "test_data/fasta_tools/inputs/unordered_taxids";

        fasta_tools::sort_taxid_dir_by_sequence_length(input_taxid_dir, temp_dir_path_str);

        for entry in std::fs::read_dir(temp_dir_path_str).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            let output_fasta_file = path.to_str().unwrap();
            let truth_fasta_file = format!(
                "{}/{}",
                output_truth_taxid_dir,
                path.file_name().unwrap().to_str().unwrap()
            );
            assert!(util::are_files_equal(output_fasta_file, &truth_fasta_file));
        }
    }

    #[test]
    fn test_shuffle_fasta() {
        let input_fasta_file = "test_data/fasta_tools/inputs/nt";

        let temp_dir = tempdir().unwrap();
        let temp_file = temp_dir.path().join("shuffled.fa");
        let temp_file_path_str = temp_file.to_str().unwrap();

        let _ = fasta_tools::shuffle_fasta_by_sequence_index(input_fasta_file, temp_file_path_str);

        // make sure the shuffled file contains the same records as the original
        util::compare_fasta_records_from_files(input_fasta_file, temp_file_path_str);

        // test that the shuffled file has accessions out of order compared to the input
        let input_records = util::create_fasta_records_from_file(input_fasta_file);
        let shuffled_records = util::create_fasta_records_from_file(temp_file_path_str);
        // the test file has 201 records so the chance is very small that the ordering will be the same
        assert_ne!(input_records, shuffled_records);
    }
}
