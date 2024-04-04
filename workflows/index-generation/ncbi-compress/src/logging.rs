pub mod logging {
    use std::io::Write;

    use chrono::Local;
    use env_logger;
    use log::LevelFilter;

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
}
