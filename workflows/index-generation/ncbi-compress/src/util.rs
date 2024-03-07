// sequence and file utility functions, mainly for testing purposes

pub mod util {
    
    use std::fs;
    use std::cmp::Ordering;
    use std::io::Read;
    use std::io::Write;

    use bio::io::fasta;
    use rand::Rng;

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

    pub fn create_fasta_records_from_file(input_fasta_path: &str) -> Vec<fasta::Record> {
        let reader = fasta::Reader::from_file(&input_fasta_path).unwrap();
        let mut records_iter = reader.records();
        let mut records: Vec<fasta::Record> = Vec::new();
        while let Some(record) = records_iter.next() {
            let record = record.unwrap();
            records.push(record);
        }
        return records;
    }

    pub fn are_files_equal(file_path1: &str, file_path2: &str) -> bool {
        if let Ok(contents1) = fs::read(file_path1) {
            if let Ok(contents2) = fs::read(file_path2) {
                return contents1 == contents2;
            }
        }
        false
    }

    pub fn are_file_records_similar(file1: &str, file2: &str) -> f64 {

        let mut records1 = create_fasta_records_from_file(file1);
        let mut records2 = create_fasta_records_from_file(file2);

        records1.sort();
        records2.sort();
    
        let total_elements = records1.len();
        let mut similar_elements = 0;
    
        for (elem1, elem2) in records1.iter().zip(records2.iter()) {
            if elem1 == elem2 {
                similar_elements += 1;
            }
        }
    
        let similarity = similar_elements as f64 / total_elements as f64;
        similarity
    }

    pub fn compare_fasta_records_from_files(file1: &str, file2: &str) {

        let mut records1 = create_fasta_records_from_file(file1);
        let mut records2 = create_fasta_records_from_file(file2);

        records1.sort();
        records2.sort();

        assert_eq!(
            records1,
            records2,
            "Records do not match",
        );
    }

    pub fn write_to_file(filename: &str, content: &str) -> std::io::Result<()> {
        let mut file = fs::File::create(filename)?;
        println!("{}", format!("Wrote to file: {}", filename));
        file.write_all(content.as_bytes())?;
        Ok(())
    }

    pub fn read_contents(input_fasta_path: &str) -> String {
        let mut file_content = String::new();
        let mut file = std::fs::File::open(&input_fasta_path).expect("Failed to open file");
        file.read_to_string(&mut file_content)
            .expect("Failed to read from file");
        file_content
    }

    pub fn read_and_write_to_file(input_fasta_path: &str, output_fasta_path: &str) {
            let file_content = read_contents(input_fasta_path);
            println!("{}", format!("writing to file: {}", output_fasta_path));
            let _ = write_to_file(&output_fasta_path, &file_content);
    }

    fn split_string_into_chunks(s: &str, chunk_size: usize) -> Vec<String> {
        let chars: Vec<char> = s.chars().collect();
        let mut chunks = Vec::new();
    
        for i in (0..s.len()).step_by(chunk_size) {
            let chunk: String = chars[i..i + chunk_size].iter().collect();
            chunks.push(chunk);
        }
        chunks
    }

    pub fn create_sequence_w_containment(ksize: usize, haystack: &str, containment: f64) -> String {
        // check that haystack can be broken down into ksize chunks evenly
        if haystack.len() % ksize != 0 {
            panic!("Haystack length must be evenly divisible by ksize")
        }
        // break apart haystack into ksize chunks
        let haystack_kmers = split_string_into_chunks(haystack, ksize);

        // calculate the number of kmers to replace from haystack to get requested containment
        let num_kmers_to_replace = (haystack_kmers.len() as f64 * (1.0 - containment)) as usize;
        println!("num_kmers_to_replace: {}", num_kmers_to_replace);
        // create a vector of random indices to replace the length of num_kmers_to_replace
        let mut rng = rand::thread_rng();
        let mut indices = Vec::new();
        for _ in 0..num_kmers_to_replace {
            indices.push(rng.gen_range(0..haystack_kmers.len()));
        }

        // replace the kmers at the random indices with random kmers
        let mut sequence = String::new();
        for (i, kmer) in haystack_kmers.into_iter().enumerate() {
            if indices.contains(&i) {
                let random_kmer = create_random_kmer(ksize);
                sequence.push_str(&random_kmer);
            } else {
                sequence.push_str(&kmer);
            }
        }
        sequence
    }

    pub fn create_random_kmer(ksize: usize) -> String {
        let letters = ['A', 'C', 'T', 'G'];
        let mut rng = rand::thread_rng();
        let mut sequence = String::new();
        for _ in 0..ksize {
            let random_letter = letters[rng.gen_range(0..letters.len())];
            sequence.push(random_letter);
        }
        sequence
    }

    pub fn create_random_sequence(ksize: usize, total_kmers: u32) -> String {
        let mut sequence = String::new();
        for _ in 0..total_kmers {
            let random_kmer = create_random_kmer(ksize);
            sequence.push_str(&random_kmer);
        }
        sequence
    }
}

