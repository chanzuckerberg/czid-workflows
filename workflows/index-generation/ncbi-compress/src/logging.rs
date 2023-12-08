
pub mod logging {
    use std::fs::OpenOptions;
    use std::io::Write;
    use std::fs;

    use csv::WriterBuilder;
    use log::LevelFilter;
    use chrono::Local;
    use env_logger;


    pub fn init_stdout_logging() {
        env_logger::Builder::new()
        .format(|buf, record| {
            writeln!(
                buf,
                "{} {}: {}",
                record.level(),
                Local::now().format("%Y-%m-%d %H:%M:%S%.3f"),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
    }

    pub fn write_to_file(data: Vec<&str>, file_name: &str) {
        // Open the file in append mode
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(file_name)
            .expect("Failed to open the file");

        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(file);

        // TODO: handle errors here
        let _ = writer.write_record(&data);

        let _ = writer.flush();
    }

    pub fn initialize_tsv(file_name: &str, headers: Vec<&str>) {

        println!("Removing file {} if it exists to start fresh", file_name);
        match fs::remove_file(file_name) {
            Ok(()) => {
                println!("File removed successfully.");
            }
            Err(_e) => {} // logging file doesn't exist so we don't need to remove it
        };

        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .open(file_name)
            .expect("Failed to open the file");

        let mut writer = WriterBuilder::new().delimiter(b'\t').from_writer(file);
        let _ = writer.write_record(&headers).expect("Failed to write headers");


        let _ = writer.flush().expect("Failed to flush writer");
    }


}
