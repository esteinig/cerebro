use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use csv::{Reader, ReaderBuilder, Writer, WriterBuilder};
use env_logger::Builder;
use log::{LevelFilter, Level};
use niffler::{get_reader, get_writer};
use plotters::coord::Shift;
use plotters::prelude::SVGBackend;
use plotters_bitmap::BitMapBackend;
use serde::{Deserialize, Serialize};
use plotters::drawing::IntoDrawingArea;
use glob::glob;
use itertools::Itertools;
use plotters::prelude::*;
use plotters::style::Color;
use std::hash::Hash;
use std::hash::Hasher;
use plotters::style::text_anchor::*;
use crate::error::CiqaError;
use crate::plate::DiagnosticData;
use crate::terminal::PlotDiagStripArgs;


#[derive(Clone, Copy, Eq, PartialEq, Hash, Debug)]
pub enum Metric { Sens, Spec }


#[derive(Clone, Debug)]
pub struct Row {
    size: String,     // "4b".."32b"
    quant: String,    // "q2".."q8"
    clinical: String, // "clinical"|"noclinical"| "na"
    scheme: String,   // "tiered"|"simple"| "na"
    metric: Metric,
    value: f64,       // 0..1
    is_consensus: bool,
    sample_key: String,
}

pub fn run_plot_diag_strip(args: PlotDiagStripArgs) -> anyhow::Result<()> {
    // collect input paths
    let mut paths = vec![];
    for pat in &args.inputs {
        for ent in glob(pat)? {
            paths.push(ent?);
        }
    }
    if paths.is_empty() {
        anyhow::bail!("no inputs matched");
    }

    // build rows
    let mut rows: Vec<Row> = vec![];
    for p in paths {
        let fname = p.file_name().and_then(|s| s.to_str()).unwrap_or_default();
        let (size, quant, clinical, scheme) = parse_diag_fname(fname)?;
        let dd = DiagnosticData::from_json(&p)?;

        // consensus rows if present
        if dd.consensus.sensitivity > 0.0 || dd.consensus.specificity > 0.0 {
            rows.push(Row { size: size.clone(), quant: quant.clone(), clinical: clinical.clone(), scheme: scheme.clone(),
                metric: Metric::Sens, value: dd.consensus.sensitivity, is_consensus: true, sample_key: "consensus".into() });
            rows.push(Row { size: size.clone(), quant: quant.clone(), clinical: clinical.clone(), scheme: scheme.clone(),
                metric: Metric::Spec, value: dd.consensus.specificity, is_consensus: true, sample_key: "consensus".into() });
        }
        // per-replicate rows
        for s in dd.stats.iter().filter(|s| s.name != "consensus") {
            rows.push(Row { size: size.clone(), quant: quant.clone(), clinical: clinical.clone(), scheme: scheme.clone(),
                metric: Metric::Sens, value: s.sensitivity, is_consensus: false, sample_key: s.name.clone() });
            rows.push(Row { size: size.clone(), quant: quant.clone(), clinical: clinical.clone(), scheme: scheme.clone(),
                metric: Metric::Spec, value: s.specificity, is_consensus: false, sample_key: s.name.clone() });
        }
    }

    // facets: clinical × scheme if available, else single panel
    let has_clin = rows.iter().any(|r| r.clinical == "clinical" || r.clinical == "noclinical");
    let has_scheme = rows.iter().any(|r| r.scheme == "tiered" || r.scheme == "simple");

    let facet_groups = rows.iter().group_by(|r| (
        if has_clin { r.clinical.clone() } else { "all".into() },
        if has_scheme { r.scheme.clone() } else { "all".into() },
    ));
    let facets: Vec<((String,String), Vec<Row>)> =
        facet_groups.into_iter().map(|(k, it)| (k, it.cloned().collect())).collect();

    // backend
    let ext = args.output.extension().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
    let is_svg = ext == "svg";
    if is_svg {
        let root = SVGBackend::new(&args.output, (args.width, args.height)).into_drawing_area();
        root.fill(&WHITE)?;
        draw_facets(root, facets, &args)?;
    } else {
        let root = BitMapBackend::new(&args.output, (args.width, args.height)).into_drawing_area();
        root.fill(&WHITE)?;
        draw_facets(root, facets, &args)?;
    }
    Ok(())
}



fn draw_facets<DB: DrawingBackend>(
    root: DrawingArea<DB, Shift>,
    facets: Vec<((String,String), Vec<Row>)>,
    args: &PlotDiagStripArgs,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {

    // layout grid
    let n = facets.len().max(1);
    let cols = n.min(3);
    let rows = (n + cols - 1) / cols;
    let mut areas = root.split_evenly((rows, cols)).into_iter();

    for ((clin, scheme), mut data) in facets {
        // enforce size ordering
        let order = ["4b","8b","14b","32b"];
        data.sort_by_key(|r| order.iter().position(|s| s==&r.size).unwrap_or(0));

        let area = areas.next().unwrap();
        draw_panel(area, &data, &clin, &scheme, args)?;
    }

    root.present()?;
    Ok(())
}

fn draw_panel<DB: DrawingBackend>(
    area: DrawingArea<DB, Shift>,
    data: &[Row],
    clin: &str,
    scheme: &str,
    args: &PlotDiagStripArgs,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {

    let sizes = ["4b","8b","14b","32b"];
    let size_idx = |s:&str| sizes.iter().position(|x| x==&s).unwrap_or(0) as i32;
    let y_max = 1.0;

    let sizes = ["4b","8b","14b","32b"];
    let idx_of = |s:&str| sizes.iter().position(|x| x==&s).unwrap_or(0) as f64;

    let x_min = 0.0;
    let x_max = (sizes.len() - 1) as f64;
    let y_min = 0.0;
    let y_max = 1.0;
    
    let mut chart =
    ChartBuilder::on(&area)
        .margin(16)
        .caption(panel_title(args.title.as_deref(), clin, scheme), ("DejaVu Sans", 16))
        .x_label_area_size(32)
        .y_label_area_size(44)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // axis and grid
    chart.configure_mesh()
        .disable_x_mesh()
        .y_desc(if args.percent { "Value (%)" } else { "Value" })
        .x_labels(sizes.len())
        .x_label_formatter(&|i| sizes[*i as usize].to_string())
        .y_label_formatter(&|v| if args.percent { format!("{:.0}", v*100.0) } else { format!("{:.2}", v) })
        .light_line_style(&RGBColor(230,230,230))
        .label_style(("DejaVu Sans", 12))
        .draw()?;

    // N per size (under ticks)
    
    let n_by_size = data.iter()
    .filter(|r| !r.is_consensus)
    .group_by(|r| r.size.clone())
    .into_iter()
    .map(|(s, g)| (s, g.count()))
    .collect::<std::collections::HashMap<_, _>>();

    for (i, s) in sizes.iter().enumerate() {
    if let Some(n) = n_by_size.get(&s.to_string()) {
        let x = i as f64;
        let y = 0.02_f64; // inside range [0,1]
        chart.draw_series(std::iter::once(
            Text::new(format!("n={}", n / 2), (x, y), ("DejaVu Sans", 10).into_font()),
        ))?;
    }
    }

    // colors
    let teal = RGBColor(91,142,125);   // sensitivity
    let brown = RGBColor(156,109,91);  // specificity
    let charcoal = RGBColor(40,50,60); // consensus + lines

    // offsets
    let off_sens = -0.12;
    let off_spec =  0.12;

    // strip points
    for (metric, col, off) in [
        (Metric::Sens, teal, off_sens),
        (Metric::Spec, brown, off_spec),
    ] {
        for r in data.iter().filter(|r| r.metric==metric && !r.is_consensus) {
            let x = size_idx(&r.size) as f64 + off + jitter_for(&r.sample_key, args.jitter);
            chart.draw_series(std::iter::once(Circle::new((x, r.value), 3, ShapeStyle::from(&col).filled())))?;
        }
    }

    // means and CI
    if args.mean || args.ci {
        for (metric, col, off) in [
            (Metric::Sens, teal, off_sens),
            (Metric::Spec, brown, off_spec),
        ] {
            for s in sizes {
                let vals: Vec<f64> = data.iter().filter(|r| r.metric==metric && r.size==s && !r.is_consensus).map(|r| r.value).collect();
                if vals.is_empty() { continue; }
                let (mean, lo, hi) = mean_ci95(&vals);
                let x = size_idx(s) as f64 + off;

                if args.mean {
                    // narrow bar from 0 to mean (muted fill)
                    let w = 0.10;
                    let rect = Rectangle::new(
                        [(x-w/2.0, 0.0), (x+w/2.0, mean)],
                        ShapeStyle::from(&col.mix(0.35)).filled(),
                    );
                    chart.draw_series(std::iter::once(rect))?;
                }
                if args.ci {
                    // whisker
                    let style = ShapeStyle::from(&col).stroke_width(1);
                    chart.draw_series(std::iter::once(PathElement::new(vec![(x, lo), (x, hi)], style)))?;
                    chart.draw_series(std::iter::once(PathElement::new(vec![(x-0.04, lo), (x+0.04, lo)], style)))?;
                    chart.draw_series(std::iter::once(PathElement::new(vec![(x-0.04, hi), (x+0.04, hi)], style)))?;
                }
            }
        }
    }

    // consensus horizontal lines
    if args.consensus {
        for metric in [Metric::Sens, Metric::Spec] {
            if let Some(v) = consensus_value(data, metric) {
                let style = ShapeStyle::from(&charcoal).stroke_width(1);
                let xs: Vec<(f64, f64)> = (0..sizes.len()).map(|i| (i as f64, v)).collect();
                chart.draw_series(std::iter::once(PathElement::new(xs, style)))?;
            }
        }
    }

    // quantization overlay lines
    if args.lines_by_quant {
        let quants = ["q2", "q4", "q8"];
        let lw = args.line_width;              // Copy for legend closure
        let col = charcoal;                    // Copy color

        for (metric, off) in [(Metric::Sens, off_sens), (Metric::Spec, off_spec)] {
            for q in quants {
                // mean per size for this quant and metric
                let series: Vec<(f64, f64)> = data.iter()
                    .filter(|r| r.metric == metric && !r.is_consensus && r.quant == q)
                    .group_by(|r| r.size.clone())
                    .into_iter()
                    .map(|(size, grp)| {
                        let vals: Vec<f64> = grp.map(|r| r.value).collect();
                        let mean = if vals.is_empty() { f64::NAN } else { vals.iter().sum::<f64>() / (vals.len() as f64) };
                        (idx_of(&size), mean)
                    })
                    .sorted_by(|a, b| a.0.total_cmp(&b.0))
                    .collect();

                if series.is_empty() { continue; }

                let path_pts: Vec<(f64, f64)> = series.into_iter()
                    .filter(|&(_, y)| y.is_finite())
                    .map(|(x, y)| (x + off, y))
                    .collect();

                let style = ShapeStyle::from(&col).stroke_width(lw);
                chart
                    .draw_series(std::iter::once(PathElement::new(path_pts, style)))?
                    .label(format!("{} {}", q, match metric { Metric::Sens => "sens", Metric::Spec => "spec" }))
                    .legend(move |(x, y)| {
                        // Rebuild style here to avoid borrowing from outer scope
                        PathElement::new(
                            vec![(x, y), (x + 20, y)],
                            ShapeStyle::from(&col).stroke_width(lw),
                        )
                    });
            }
        }

        chart.configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.9))
            .border_style(&RGBColor(200, 200, 200))
            .label_font(("DejaVu Sans", 10))
            .draw()?;
    }

    Ok(())
}

// === helpers ===

fn panel_title(main: Option<&str>, clin:&str, scheme:&str) -> String {
    match (main, clin, scheme) {
        (Some(t), "all", "all") => t.to_string(),
        (Some(t), c, "all") => format!("{t} — {c}"),
        (Some(t), "all", s) => format!("{t} — {s}"),
        (Some(t), c, s) => format!("{t} — {c}, {s}"),
        (None, c, s) => {
            if c=="all" && s=="all" { "Diagnostic performance".into() }
            else if s=="all" { format!("Diagnostic performance — {}", c) }
            else if c=="all" { format!("Diagnostic performance — {}", s) }
            else { format!("Diagnostic performance — {}, {}", c, s) }
        }
    }
}

fn jitter_for(key:&str, half_width:f64)->f64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    key.hash(&mut h);
    let u = (h.finish() as f64 % 10_000.0) / 10_000.0; // 0..1
    (u*2.0 - 1.0) * half_width
}

fn consensus_value(data:&[Row], metric:Metric)->Option<f64>{
    data.iter().find(|r| r.is_consensus && r.metric==metric).map(|r| r.value)
}

fn mean_ci95(xs:&[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let mean = xs.iter().copied().sum::<f64>() / n.max(1.0);
    if xs.len() < 2 { return (mean, mean, mean); }
    let var = xs.iter().map(|v| (v-mean)*(v-mean)).sum::<f64>() / (n-1.0);
    let se = (var / n).sqrt();
    // use normal approx 1.96; swap to t if you like
    let z = 1.96;
    (mean, (mean - z*se).max(0.0), (mean + z*se).min(1.0))
}

// qwen3-14b-q2-kl_clinical_tiered_tiered-threshold_default.ct.json
fn parse_diag_fname(fname:&str)->anyhow::Result<(String,String,String,String)>{
    let stem = fname.strip_suffix(".ct.json").unwrap_or(fname);
    let (head, rest) = stem.split_once('_').ok_or_else(|| anyhow::anyhow!("bad name: {fname}"))?;
    let mut head_parts = head.split('-');
    let _model = head_parts.next().unwrap_or("qwen3");
    let size = head_parts.next().unwrap_or("14b").to_string();
    let quant = head_parts.next().unwrap_or("q2").to_string();
    let _ctx = head_parts.next().unwrap_or("kl");

    let mut rparts = rest.split('_');
    let clinical = rparts.next().unwrap_or("na").to_string();
    let scheme = rparts.next().unwrap_or("na").to_string();
    Ok((size, quant, clinical, scheme))
}

pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}


/// Enum to specify the type of file component to retrieve
pub enum FileComponent {
    /// The full file name including the extension
    FileName,
    /// The file name without the extension
    FileStem,
}

/// Extracts the specified file component from a `PathBuf` and returns it as a `String`.
///
/// # Arguments
///
/// * `path` - A `PathBuf` representing the file path.
/// * `component` - A `FileComponent` specifying whether to get the file name or the file stem.
///
/// # Returns
///
/// * `Result<String, DatabaseError>` - A `Result` containing the specified file component as a `String`
///   if successful, or a `DatabaseError` if an error occurs.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use cerebro_ciqa::utils::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(&path, FileComponent::FileName) {
///     Ok(file_name) => println!("File name: {}", file_name),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
///
/// ```
/// use std::path::PathBuf;
/// use cerebro_ciqa::utils::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(&path, FileComponent::FileStem) {
///     Ok(file_stem) => println!("File stem: {}", file_stem),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn get_file_component(path: &PathBuf, component: FileComponent) -> Result<String, CiqaError> {
    match component {
        FileComponent::FileName => {
            path.file_name()
                .ok_or(CiqaError::FileNameConversionError(path.to_path_buf()))
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(CiqaError::FileNameConversionError(path.to_path_buf())))
        }
        FileComponent::FileStem => {
            path.file_stem()
                .ok_or(CiqaError::FileNameConversionError(path.to_path_buf()))
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(CiqaError::FileNameConversionError(path.to_path_buf())))
        }
    }
}



pub trait StringUtils {
    fn substring(&self, start: usize, len: usize) -> Self;
}

impl StringUtils for String {
    fn substring(&self, start: usize, len: usize) -> Self {
        self.chars().skip(start).take(len).collect()
    }
}


pub trait UuidUtils {
    fn shorten(&self, len: usize) -> String;
}

impl UuidUtils for uuid::Uuid {
    fn shorten(&self, len: usize) -> String {
        self.to_string().substring(0, len)
    }
}

pub fn init_logger() {

    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(env_logger::fmt::Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(env_logger::fmt::Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(env_logger::fmt::Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style.set_color(env_logger::fmt::Color::Rgb(255, 102, 0)).set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style.set_color(env_logger::fmt::Color::Rgb(255, 195, 0)).set_bold(true);

            let msg = match record.level(){
                Level::Warn => (orange_style.value(record.level()), orange_style.value(record.args())),
                Level::Info => (green_style.value(record.level()), white_style.value(record.args())),
                Level::Debug => (apricot_style.value(record.level()), apricot_style.value(record.args())),
                Level::Error => (red_style.value(record.level()), red_style.value(record.args())),
                _ => (white_style.value(record.level()), white_style.value(record.args()))
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

pub fn get_tsv_reader(file: &Path, flexible: bool, header: bool) -> Result<Reader<Box<dyn Read>>, CiqaError> {

    let buf_reader = BufReader::new(File::open(&file)?);
    let (reader, _format) = get_reader(Box::new(buf_reader))?;

    let tsv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .flexible(flexible) // Allows records with a different number of fields
        .from_reader(reader);

    Ok(tsv_reader)
}

pub fn get_tsv_writer(file: &Path, header: bool) -> Result<Writer<Box<dyn Write>>, CiqaError> {
    
    let buf_writer = BufWriter::new(File::create(&file)?);
    let writer = get_writer(Box::new(buf_writer), niffler::Format::from_path(file), niffler::compression::Level::Six)?;

    let csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .from_writer(writer);

    Ok(csv_writer)
}

pub fn write_tsv<T: Serialize>(data: &Vec<T>, file: &Path, header: bool) -> Result<(), CiqaError> {

    let mut writer = get_tsv_writer(file, header)?;

    for value in data {
        // Serialize each value in the vector into the writer
        writer.serialize(&value)?;
    }

    // Flush and complete writing
    writer.flush()?;
    
    Ok(())
}

pub fn read_tsv<T: for<'de>Deserialize<'de>>(file: &Path, flexible: bool, header: bool) -> Result<Vec<T>, CiqaError> {

    let mut reader = get_tsv_reader(file, flexible, header)?;

    let mut records = Vec::new();
    for record in reader.deserialize() {
        records.push(record?)
    }

    Ok(records)
}

