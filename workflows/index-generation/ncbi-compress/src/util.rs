pub mod util {
    use std::fs;
    use std::io::Read;
    use std::io::Write;

    use rand::Rng;

    pub fn are_files_equal(file_path1: &str, file_path2: &str) -> bool {
        if let Ok(contents1) = fs::read(file_path1) {
            if let Ok(contents2) = fs::read(file_path2) {
                return contents1 == contents2;
            }
        }
        false
    }

    pub fn write_to_file(filename: &str, content: &str) -> std::io::Result<()> {
        println!("starting write to file in write_to_file");
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
        return file_content
    }

    pub fn read_and_write_to_file(input_fasta_path: &str, output_fasta_path: &str) -> () {
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
        return sequence
    }

    pub fn create_random_kmer(ksize: usize) -> String {
        let letters = ['A', 'C', 'T', 'G'];
        let mut rng = rand::thread_rng();
        let mut sequence = String::new();
        for _ in 0..ksize {
            let random_letter = letters[rng.gen_range(0..letters.len())];
            sequence.push(random_letter);
        }
        return sequence
    }

    pub fn create_random_sequence(ksize: usize, total_kmers: u32) -> String {
        let mut sequence = String::new();
        for _ in 0..total_kmers {
            let random_kmer = create_random_kmer(ksize);
            sequence.push_str(&random_kmer);
        }
        return sequence
    }
}

