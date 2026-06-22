use env_logger::fmt::Color;
use env_logger::Builder;
use log::{Level, LevelFilter};
use std::ffi::OsStr;
use std::io::Write;
use std::path::Path;

pub const _CRATE_VERSION: &'static str = env!(
    "CARGO_PKG_VERSION",
    "Failed to get the crate version at compile time - this is not good!"
);

pub fn init_logger() {
    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style
                .set_color(Color::Rgb(255, 102, 0))
                .set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style
                .set_color(Color::Rgb(255, 195, 0))
                .set_bold(true);

            let msg = match record.level() {
                Level::Warn => (
                    orange_style.value(record.level()),
                    orange_style.value(record.args()),
                ),
                Level::Info => (
                    green_style.value(record.level()),
                    white_style.value(record.args()),
                ),
                Level::Debug => (
                    apricot_style.value(record.level()),
                    apricot_style.value(record.args()),
                ),
                Level::Error => (
                    red_style.value(record.level()),
                    red_style.value(record.args()),
                ),
                _ => (
                    white_style.value(record.level()),
                    white_style.value(record.args()),
                ),
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                white_style.value(timestamp),
                msg.0,
                msg.1
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
}
