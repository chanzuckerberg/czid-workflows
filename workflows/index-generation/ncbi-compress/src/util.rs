pub mod util {
    use std::fs;
    use std::io::Read;
    use std::io::Write;

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
            let mut file_content = read_contents(input_fasta_path);
            println!("{}", format!("writing to file: {}", output_fasta_path));
            let _ = write_to_file(&output_fasta_path, &file_content);
    }
}

