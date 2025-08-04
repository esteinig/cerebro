use std::{
    io::{self, Stdout},
    time::{Duration, Instant},
};

use clap::{Parser, Subcommand};
use ratatui::crossterm::{
    event::{self, Event, KeyCode, KeyEvent, KeyModifiers},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style, Stylize},
    symbols,
    text::{Line, Span, Text},
    widgets::{
        Axis, Bar, BarChart, BarGroup, Block, Borders, Cell, Chart, Dataset, GraphType,
        LegendPosition, Paragraph, Row, Table, TableState,
    },
    Terminal,
};
use rand::{SeedableRng, rngs::StdRng, Rng};

#[derive(Copy, Clone)]
enum CellClass { TP, FP, TN, FN, Excluded, Control }

impl CellClass {
    fn glyph(self) -> &'static str {
        match self {
            CellClass::TP => "●",
            CellClass::FP => "●",
            CellClass::TN => "●",
            CellClass::FN => "●",
            CellClass::Excluded => "·",
            CellClass::Control => "◆",
        }
    }
}

#[derive(Clone)]
struct Model {
    name: String,
    gpu_times: Vec<(String, u64)>,
    run_times: Vec<(String, u64)>,
    parameters: String,
    size_gb: f32,
    tp: u64,
    fp: u64,
    tn: u64,
    fn_: u64,
    // classification matrix data
    replicates: usize,          // number of replicate columns displayed
    sample_classes: Vec<CellClass>,            // length 96 consensus (final classification)
    replicate_classes: Vec<CellClass>,         // length 96 * replicates
}

impl Model {
    fn gpu_stats(&self) -> (u64, f64, u64) {
        stats(&self.gpu_times.iter().map(|(_, v)| *v).collect::<Vec<_>>())
    }
    fn run_stats(&self) -> (u64, f64, u64) {
        stats(&self.run_times.iter().map(|(_, v)| *v).collect::<Vec<_>>())
    }
    fn sensitivity(&self) -> f64 {
        let denom = self.tp + self.fn_;
        if denom == 0 { 0.0 } else { self.tp as f64 / denom as f64 }
    }
    fn specificity(&self) -> f64 {
        let denom = self.tn + self.fp;
        if denom == 0 { 0.0 } else { self.tn as f64 / denom as f64 }
    }
    fn replicate_cell(&self, r: usize, sample_idx: usize) -> CellClass {
        // replicate index r in [0, replicates)
        let base = r * 96;
        self.replicate_classes[base + sample_idx]
    }
}

#[derive(Clone, Debug)]
struct Perf {
    name: String,
    sensitivity: f64,
    specificity: f64,
}

fn perf_for_models(models: &[Model]) -> Vec<Perf> {
    models.iter().map(|m| Perf {
        name: m.name.clone(),
        sensitivity: m.sensitivity(),
        specificity: m.specificity(),
    }).collect()
}

fn stats(vals: &[u64]) -> (u64, f64, u64) {
    if vals.is_empty() { return (0, 0.0, 0); }
    let mut v = vals.to_vec();
    v.sort_unstable();
    let min = v[0];
    let max = v[v.len()-1];
    let median = if v.len() % 2 == 0 {
        let a = v[v.len()/2 - 1];
        let b = v[v.len()/2];
        (a as f64 + b as f64)/2.0
    } else { v[v.len()/2] as f64 };
    (min, median, max)
}

// ---------- Animated Sine ----------
#[derive(Clone)]
struct SinSignal { x: f64, interval: f64, period: f64, amplitude: f64, phase: f64 }
impl SinSignal { const fn new(interval: f64, period: f64, amplitude: f64, phase: f64) -> Self { Self { x:0.0, interval, period, amplitude, phase } } }
impl Iterator for SinSignal {
    type Item=(f64,f64);
    fn next(&mut self)->Option<Self::Item>{
        let y=((self.x/self.period)*std::f64::consts::TAU + self.phase).sin()*self.amplitude;
        let p=(self.x,y); self.x+=self.interval; Some(p)
    }
}

// ---------- Themes ----------
#[derive(Clone)]
struct PastelTheme {
    header_models: Style,
    header_metrics: Style,
    header_stats: Style,
    header_detail: Style,
    footer: Style,
    summary: Style,
    perf_sens: Style,
    perf_spec: Style,
    dna_a: Style,
    dna_b: Style,
    matrix_tp: Style,
    matrix_tn: Style,
    matrix_fp: Style,
    matrix_fn: Style,
    matrix_ex: Style,
    matrix_ctrl: Style,
    matrix_consensus: Style,
    matrix_header: Style,
}

fn hex_rgb(hex: u32) -> Color {
    Color::Rgb(((hex>>16)&0xFF) as u8, ((hex>>8)&0xFF) as u8, (hex&0xFF) as u8)
}

fn make_themes()->Vec<PastelTheme>{
    let base = PastelTheme{
        header_models: Style::default().fg(hex_rgb(0xFF66CC)).add_modifier(Modifier::BOLD),
        header_metrics: Style::default().fg(hex_rgb(0x66FFCC)).add_modifier(Modifier::BOLD),
        header_stats: Style::default().fg(hex_rgb(0xC4A1FF)).add_modifier(Modifier::BOLD),
        header_detail: Style::default().fg(hex_rgb(0xFFBF80)).add_modifier(Modifier::BOLD),
        footer: Style::default().fg(hex_rgb(0x80859A)).add_modifier(Modifier::DIM),
        summary: Style::default().fg(hex_rgb(0xFFF8A6)),
        perf_sens: Style::default().fg(hex_rgb(0x5BE3FF)),
        perf_spec: Style::default().fg(hex_rgb(0xF26BFF)),
        dna_a: Style::default().fg(hex_rgb(0x5BE3FF)),
        dna_b: Style::default().fg(hex_rgb(0xF26BFF)),
        matrix_tp: Style::default().fg(hex_rgb(0x66FFC9)).add_modifier(Modifier::BOLD),
        matrix_tn: Style::default().fg(hex_rgb(0x66FFC9)),
        matrix_fp: Style::default().fg(hex_rgb(0xFF5A9E)),
        matrix_fn: Style::default().fg(hex_rgb(0xFFCF5A)),
        matrix_ex: Style::default().fg(Color::DarkGray),
        matrix_ctrl: Style::default().fg(hex_rgb(0x53E0FF)),
        matrix_consensus: Style::default().fg(hex_rgb(0xFFFFFF)).add_modifier(Modifier::BOLD),
        matrix_header: Style::default().fg(hex_rgb(0x8890A0)).add_modifier(Modifier::BOLD),
    };
    let alt = PastelTheme { header_models: Style::default().fg(hex_rgb(0xFF75D1)).add_modifier(Modifier::BOLD), ..base.clone() };
    vec![base, alt]
}
struct App {
   pub models: Vec<Model>,
   pub selected: TableState,
   pub themes: Vec<PastelTheme>,
   pub theme_index: usize,
   pub strand_a: SinSignal,
   pub strand_b: SinSignal,
   pub dna_a: Vec<(f64,f64)>,
   pub dna_b: Vec<(f64,f64)>,
   pub window: [f64;2],
   pub perf_mode: PerfMode,
}

impl App {
    fn new()->Self{
        let replicate_count = 6;
        let mut models = vec![
            Model{
                name:"qwen3-32b-q8-0".into(),
                gpu_times: vec![("rep1".into(),210),("rep2".into(),205),("rep3".into(),220),("rep4".into(),215)],
                run_times: vec![("rep1".into(),62),("rep2".into(),64),("rep3".into(),60),("rep4".into(),63)],
                parameters:"32B".into(), size_gb:23.5,
                tp:120, fp:0, tn:130, fn_:10,
                replicates: replicate_count,
                sample_classes: vec![],
                replicate_classes: vec![],
            },
            Model{
                name:"qwen3-8b-q8-0".into(),
                gpu_times: vec![("rep1".into(),95),("rep2".into(),98),("rep3".into(),92),("rep4".into(),97)],
                run_times: vec![("rep1".into(),25),("rep2".into(),26),("rep3".into(),24),("rep4".into(),25)],
                parameters:"8B".into(), size_gb:6.1,
                tp:110, fp:1, tn:129, fn_:12,
                replicates: replicate_count,
                sample_classes: vec![],
                replicate_classes: vec![],
            },
            Model{
                name:"deepseekr1-qwen32b-q8-0".into(),
                gpu_times: vec![("rep1".into(),260),("rep2".into(),255),("rep3".into(),268),("rep4".into(),259)],
                run_times: vec![("rep1".into(),72),("rep2".into(),74),("rep3".into(),70),("rep4".into(),73)],
                parameters:"32B hybrid".into(), size_gb:24.7,
                tp:125, fp:1, tn:129, fn_:8,
                replicates: replicate_count,
                sample_classes: vec![],
                replicate_classes: vec![],
            },
        ];

        // Populate synthetic classification matrices
        for (mi, m) in models.iter_mut().enumerate() {
            synthesize_matrix(m, 96, replicate_count, mi as u64 + 42);
        }

        let mut state=TableState::default(); state.select(Some(0));
        let mut strand_a=SinSignal::new(0.15,4.0,10.0,0.0);
        let mut strand_b=SinSignal::new(0.15,4.0,10.0,std::f64::consts::PI);
        let dna_a:Vec<_>=strand_a.by_ref().take(200).collect();
        let dna_b:Vec<_>=strand_b.by_ref().take(200).collect();
        Self{
            models,
            selected:state,
            themes:make_themes(),
            theme_index:0,
            strand_a,
            strand_b,
            dna_a,
            dna_b,
            window:[0.0,30.0],
            perf_mode: PerfMode::Both
        }
    }
    fn selected_model(&self)->&Model{
        let idx=self.selected.selected().unwrap_or(0); &self.models[idx]
    }
    fn next_model(&mut self){
        let i=match self.selected.selected(){
            Some(i) if i+1<self.models.len()=>i+1,_=>0
        }; self.selected.select(Some(i));
    }
    fn prev_model(&mut self){
        let i=match self.selected.selected(){
            Some(0)|None=>self.models.len()-1,
            Some(i)=>i-1
        }; self.selected.select(Some(i));
    }
    fn cycle_theme(&mut self){ self.theme_index=(self.theme_index+1)%self.themes.len(); }
    fn theme(&self)->&PastelTheme{ &self.themes[self.theme_index] }
    fn on_tick(&mut self){
        for _ in 0..5 {
            if let Some(p)=self.strand_a.next(){ self.dna_a.push(p); }
            if let Some(p)=self.strand_b.next(){ self.dna_b.push(p); }
        }
        let start=self.window[0];
        self.dna_a.retain(|(x,_)|*x>=start-1.0);
        self.dna_b.retain(|(x,_)|*x>=start-1.0);
        self.window[0]+=0.5; self.window[1]+=0.5;
    }
}
fn synthesize_matrix(m: &mut Model, samples: usize, replicates: usize, seed: u64) {
    use rand::{SeedableRng, rngs::StdRng, Rng};
    let mut rng = StdRng::seed_from_u64(seed);

    let sens = m.sensitivity();   // copy scalar values
    let spec = m.specificity();

    // Build consensus classifications first without referencing m after mutation
    let mut sample_classes: Vec<CellClass> = Vec::with_capacity(samples);

    for i in 0..samples {
        let base = rng.random::<f64>();
        let p_pos = sens * 0.6 + 0.2;
        let class = if base < 0.05 {
            CellClass::Excluded
        } else if base < 0.08 {
            CellClass::Control
        } else if base < p_pos {
            if rng.random::<f64>() < 0.9 { CellClass::TP } else { CellClass::FN }
        } else {
            if rng.random::<f64>() < 0.95 { CellClass::TN } else { CellClass::FP }
        };
        // deterministic tweak
        let class = if i % 37 == 0 { CellClass::Control } else { class };
        sample_classes.push(class);
    }

    // Assign consensus to model now
    m.sample_classes = sample_classes;

    // Now build replicate classes referencing the (already assigned) consensus vector
    let consensus_ref = m.sample_classes.clone(); // clone so closure doesn’t borrow m mutably
    let mut replicate_classes: Vec<CellClass> = Vec::with_capacity(samples * replicates);

    for _r in 0..replicates {
        for s in 0..samples {
            let consensus = consensus_ref[s];
            let disagree = match consensus {
                CellClass::Excluded | CellClass::Control => false,
                _ => rng.random_bool(0.05),
            };
            let cell = if !disagree {
                consensus
            } else {
                match consensus {
                    CellClass::TP => CellClass::FN,
                    CellClass::FN => CellClass::TP,
                    CellClass::TN => CellClass::FP,
                    CellClass::FP => CellClass::TN,
                    other => other,
                }
            };
            replicate_classes.push(cell);
        }
    }

    m.replicate_classes = replicate_classes;
}


pub fn start_tui()->io::Result<()>{
    let mut app=App::new();
    enable_raw_mode()?;
    let mut stdout=io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend=CrosstermBackend::new(stdout);
    let mut terminal:Terminal<CrosstermBackend<Stdout>>=Terminal::new(backend)?;
    let tick=Duration::from_millis(150);
    let mut last=Instant::now();
    let res=loop{
        terminal.draw(|f| ui(f,&app))?;
        let timeout=tick.checked_sub(last.elapsed()).unwrap_or(Duration::from_secs(0));
        if event::poll(timeout)? {
            if let Event::Key(k)=event::read()? {
                if handle_key(&mut app,k)? { break Ok(()); }
            }
        }
        if last.elapsed()>=tick { app.on_tick(); last=Instant::now(); }
    };
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;
    res
}

fn handle_key(app:&mut App, key:KeyEvent)->io::Result<bool>{
    if key.modifiers.contains(KeyModifiers::CONTROL) { return Ok(false); }
    match key.code {
        KeyCode::Char('q')|KeyCode::Esc => return Ok(true),
        KeyCode::Down|KeyCode::Char('j') => app.next_model(),
        KeyCode::Up|KeyCode::Char('k') => app.prev_model(),
        KeyCode::Char('m') => app.cycle_theme(),
        _=>{}
    }
    Ok(false)
}

fn ui(f:&mut ratatui::Frame, app:&App){
    let area=f.area();
    let theme=app.theme();
    
    let main=Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(70), Constraint::Percentage(30)])
        .split(area);

    let left=Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(30), Constraint::Percentage(70)])
        .split(main[0]);

    let upper_split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(left[0]);

    let lower_split=Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(100)])
        .split(left[1]);

    let right=Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(100)])
        .split(main[1]);

    draw_models_table(f,upper_split[0],app,theme);
    draw_statistics(f,lower_split[0],app,theme);
    draw_detail(f,upper_split[1],app,theme);
    draw_performance(f,right[0],app,theme);
}

fn draw_models_table(f:&mut ratatui::Frame, area:Rect, app:&App, theme:&PastelTheme){
    let header=Row::new(vec![
        Cell::from(Span::styled("Model", theme.header_models)),
        Cell::from(Span::styled("Params", theme.header_models)),
        Cell::from(Span::styled("Size (GB)", theme.header_models)),
    ])
    .style(Style::default().add_modifier(Modifier::BOLD))
    .height(1);

    let mut state=app.selected.clone();
    let rows=app.models.iter().enumerate().map(|(i,m)|{
        let style= if Some(i)==app.selected.selected(){
            Style::default().fg(Color::Black).bg(Color::Rgb(255,228,196)).add_modifier(Modifier::BOLD)
        } else { Style::default() };
        Row::new(vec![
            Cell::from(m.name.clone()),
            Cell::from(m.parameters.clone()),
            Cell::from(format!("{:.1}", m.size_gb)),
        ]).style(style).height(1)
    });

    let table=Table::new(rows, [Constraint::Percentage(50), Constraint::Length(10), Constraint::Length(10)])
        .header(header)
        .block(Block::default().borders(Borders::ALL)
            .title(Span::styled("Models", theme.header_models)))
        .highlight_symbol("▶ ")
        .row_highlight_style(Style::default().fg(Color::Black).bg(Color::Rgb(255,235,205)).add_modifier(Modifier::BOLD));
    f.render_stateful_widget(table, area, &mut state);
}
fn draw_statistics(f:&mut ratatui::Frame, area:Rect, app:&App, theme:&PastelTheme){
    
    let chunks=Layout::default().direction(Direction::Vertical)
        .constraints([Constraint::Percentage(80), Constraint::Percentage(20)]).split(area);
    
    // Classification matrix
    let class_inner=Rect{ x:chunks[0].x+1, y:chunks[0].y+1, width:chunks[0].width.saturating_sub(2), height:chunks[0].height.saturating_sub(2) };
    draw_class_matrix(f, class_inner, app.selected_model(), theme);

    // Dna helix animation
    draw_dna_chart(f,chunks[1], app, theme);

    f.render_widget(Block::default().borders(Borders::ALL)
        .title(Span::styled("Reference Plate", theme.header_stats)), area);
}

fn draw_dna_chart(f:&mut ratatui::Frame, area:Rect, app:&App, theme:&PastelTheme){
    let datasets=vec![
        Dataset::default().name("Strand A").marker(symbols::Marker::Braille)
            .style(theme.dna_a).graph_type(GraphType::Line).data(&app.dna_a),
        Dataset::default().name("Strand B").marker(symbols::Marker::Braille)
            .style(theme.dna_b).graph_type(GraphType::Line).data(&app.dna_b),
    ];
    let chart=Chart::new(datasets)
        .block(Block::default().borders(Borders::ALL))
        .x_axis(Axis::default().title("t").style(Style::default().fg(Color::Gray)).bounds(app.window).labels([
            Span::raw(format!("{:.0}", app.window[0])),
            Span::raw(format!("{:.0}", (app.window[0]+app.window[1])/2.0)),
            Span::raw(format!("{:.0}", app.window[1])),
        ]))
        .y_axis(Axis::default().title("Amplitude").style(Style::default().fg(Color::Gray))
            .bounds([-12.0,12.0]).labels(["-10","0","10"]))
        .legend_position(Some(LegendPosition::TopRight));
    f.render_widget(chart, area);
}

fn draw_line_demo(f:&mut ratatui::Frame, area:Rect){
    let ds=Dataset::default().name("Line 2pts").marker(symbols::Marker::Braille)
        .style(Style::default().fg(Color::Yellow)).graph_type(GraphType::Line)
        .data(&[(1.,1.),(4.,4.)]);
    let chart=Chart::new(vec![ds])
        .block(Block::default().borders(Borders::ALL).title(Line::from("Line").cyan().bold().centered()))
        .x_axis(Axis::default().bounds([0.0,5.0]).labels(["0","2.5","5"]).style(Style::default().fg(Color::Gray)))
        .y_axis(Axis::default().bounds([0.0,5.0]).labels(["0","2.5","5"]).style(Style::default().fg(Color::Gray)));
    f.render_widget(chart, area);
}

fn draw_bar_demo(f:&mut ratatui::Frame, area:Rect){
    let ds=Dataset::default().graph_type(GraphType::Bar).marker(symbols::Marker::HalfBlock)
        .style(Style::default().fg(Color::Blue)).data(&[
            (0.,0.4),(10.,2.9),(20.,13.5),(30.,41.1),(40.,80.1),(50.,100.0),
            (60.,80.1),(70.,41.1),(80.,13.5),(90.,2.9),(100.,0.4),
        ]);
    let chart=Chart::new(vec![ds])
        .block(Block::default().borders(Borders::ALL).title(Line::from("Bar").cyan().bold().centered()))
        .x_axis(Axis::default().bounds([0.0,100.0]).labels(["0","50","100"]).style(Style::default().fg(Color::Gray)))
        .y_axis(Axis::default().bounds([0.0,100.0]).labels(["0","50","100"]).style(Style::default().fg(Color::Gray)));
    f.render_widget(chart, area);
}

fn draw_scatter_demo(f:&mut ratatui::Frame, area:Rect){
    let ds1=Dataset::default().name("S").graph_type(GraphType::Scatter).marker(symbols::Marker::Dot)
        .style(Style::default().fg(Color::Cyan)).data(&[(1.,3.),(2.,4.),(3.,5.),(4.,4.)]);
    let ds2=Dataset::default().name("M").graph_type(GraphType::Scatter).marker(symbols::Marker::Braille)
        .style(Style::default().fg(Color::Magenta)).data(&[(1.,1.),(2.,1.5),(3.,2.),(4.,1.8)]);
    let ds3=Dataset::default().name("H").graph_type(GraphType::Scatter).marker(symbols::Marker::Block)
        .style(Style::default().fg(Color::Yellow)).data(&[(1.,2.),(2.,2.2),(3.,2.6),(4.,3.0)]);
    let chart=Chart::new(vec![ds1,ds2,ds3])
        .block(Block::default().borders(Borders::ALL).title(Line::from("Scatter").on_cyan().bold().centered()))
        .x_axis(Axis::default().bounds([0.0,5.0]).labels(["0","2.5","5"]).style(Style::default().fg(Color::Gray)))
        .y_axis(Axis::default().bounds([0.0,6.0]).labels(["0","3","6"]).style(Style::default().fg(Color::Gray)))
        .legend_position(Some(LegendPosition::BottomLeft));
    f.render_widget(chart, area);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerfMode {
    Both,
    SensOnly,
    SpecOnly,
}
impl PerfMode {
    pub fn cycle(self) -> Self {
        match self {
            PerfMode::Both => PerfMode::SensOnly,
            PerfMode::SensOnly => PerfMode::SpecOnly,
            PerfMode::SpecOnly => PerfMode::Both,
        }
    }
}

fn build_performance_barchart(app: &App, theme: &PastelTheme) -> BarChart<'static> {
    // Get & sort model performances
    let mut perfs = perf_for_models(&app.models);
    perfs.sort_by(|a, b| b.sensitivity.partial_cmp(&a.sensitivity).unwrap());

    let selected = app.selected_model().name.as_str();
    let maxl = 9usize;

    // Build groups: one group per model
    let groups: Vec<BarGroup<'static>> = perfs
        .into_iter()
        .map(|p| {
            let short = if p.name.len() > maxl {
                format!("{}…", &p.name[..maxl.saturating_sub(1)])
            } else {
                p.name.clone()
            };

            let sel = p.name == selected;
            let sens = (p.sensitivity * 100.0).round() as u64;
            let spec = (p.specificity * 100.0).round() as u64;

            let sens_style = theme
                .perf_sens
                .add_modifier(if sel { Modifier::BOLD } else { Modifier::DIM });
            let spec_style = theme
                .perf_spec
                .add_modifier(if sel { Modifier::BOLD } else { Modifier::DIM });

            // Bars inside the group, depending on current mode
            let mut bars = Vec::with_capacity(2);
            if matches!(app.perf_mode, PerfMode::Both | PerfMode::SensOnly) {
                bars.push(
                    Bar::default()
                        .label(Line::from("S"))
                        .value(sens)
                        .text_value(format!("{sens}%"))
                        .style(sens_style)
                        .value_style(sens_style.reversed()),
                );
            }
            if matches!(app.perf_mode, PerfMode::Both | PerfMode::SpecOnly) {
                bars.push(
                    Bar::default()
                        .label(Line::from("P"))
                        .value(spec)
                        .text_value(format!("{spec}%"))
                        .style(spec_style)
                        .value_style(spec_style.reversed()),
                );
            }

            BarGroup::default()
                .label(Line::from(short))
                .bars(&bars)
        })
        .collect();

    let title = match app.perf_mode {
        PerfMode::Both => "Model Comparison (Sens & Spec)",
        PerfMode::SensOnly => "Model Comparison (Sensitivity)",
        PerfMode::SpecOnly => "Model Comparison (Specificity)",
    };

    let mut chart = BarChart::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(Line::from(title).centered().add_modifier(Modifier::BOLD)),
        )
        .direction(Direction::Horizontal)
        .bar_width(1)   // width of each bar
        .bar_gap(0)     // gap between bars inside a group
        .group_gap(2);   // gap between groups (requires ratatui >= 0.26)

    for group in groups {
        chart = chart.data(group);
    }

    chart
        
}

// fn build_performance_barchart(app:&App, theme:&PastelTheme)->BarChart<'static>{
//     let mut perfs=perf_for_models(&app.models);
//     perfs.sort_by(|a,b| b.sensitivity.partial_cmp(&a.sensitivity).unwrap());
//     let selected=app.selected_model().name.as_str();
//     let mut bars=Vec::with_capacity(perfs.len()*2);
//     for p in perfs {
//         let sens=(p.sensitivity*100.0).round() as u64;
//         let spec=(p.specificity*100.0).round() as u64;
//         let sel=p.name==selected;
//         let sens_style=theme.perf_sens.add_modifier(if sel {Modifier::BOLD}else{Modifier::DIM});
//         let spec_style=theme.perf_spec.add_modifier(if sel {Modifier::BOLD}else{Modifier::DIM});
//         let maxl=9usize;
//         let short=if p.name.len()>maxl { format!("{}…",&p.name[..maxl.saturating_sub(1)]) } else { p.name.clone() };
//         bars.push(Bar::default().value(sens).label(Line::from(format!("{short} S"))).text_value(format!("{sens}%")).style(sens_style).value_style(sens_style.reversed()));
//         bars.push(Bar::default().value(spec).label(Line::from(format!("{short} P"))).text_value(format!("{spec}%")).style(spec_style).value_style(spec_style.reversed()));
//     }
//     BarChart::default()
//         .block(Block::default().borders(Borders::ALL)
//             .title(Line::from("Model Comparison").centered().add_modifier(Modifier::BOLD)))
//         .data(BarGroup::default().bars(&bars))
//         .direction(Direction::Horizontal)
//         .bar_width(1)
//         .bar_gap(1)
// }

fn draw_performance(f:&mut ratatui::Frame, area:Rect, app:&App, theme:&PastelTheme){
    f.render_widget(build_performance_barchart(app, theme), area);
}

fn draw_confusion_matrix(f:&mut ratatui::Frame, area:Rect, model:&Model, theme:&PastelTheme){
    
    let inner= Rect{ x:area.x+1, y:area.y+1, width:area.width.saturating_sub(2), height:area.height.saturating_sub(2) };

    let rows= Layout::default().direction(Direction::Vertical)
        .constraints([Constraint::Length(1),Constraint::Length(1),Constraint::Length(1),Constraint::Min(0)])
        .split(inner);

    if rows.len()<3 { return; }

    let header_cols=Layout::default().direction(Direction::Horizontal)
        .constraints([Constraint::Length(7),Constraint::Percentage(50),Constraint::Percentage(50)]).split(rows[0]);

    let mut render=|r:Rect,s:&str,st:Style| { 
        f.render_widget(Paragraph::new(s).style(st), r); 
    };

    if header_cols.len()>=3 {
        render(header_cols[0],"", theme.matrix_header);
        render(header_cols[1],"True", theme.matrix_header);
        render(header_cols[2],"False", theme.matrix_header);
    }

    let tcols=Layout::default().direction(Direction::Horizontal)
        .constraints([Constraint::Length(7),Constraint::Percentage(50),Constraint::Percentage(50)]).split(rows[1]);

    let fcols=Layout::default().direction(Direction::Horizontal)
        .constraints([Constraint::Length(7),Constraint::Percentage(50),Constraint::Percentage(50)]).split(rows[2]);

    let tp=model.tp; let fn_=model.fn_; let fp=model.fp; let tn=model.tn;
    let sens=model.sensitivity(); let spec=model.specificity();

    if tcols.len()>=3 {
        render(tcols[0],"True", Style::default().fg(Color::White).add_modifier(Modifier::BOLD));
        render(tcols[1], &format!("TP {tp}"), theme.matrix_tp);
        render(tcols[2], &format!("FN {fn_}"), theme.matrix_fn);
    }
    if fcols.len()>=3 {
        render(fcols[0],"False", Style::default().fg(Color::White).add_modifier(Modifier::BOLD));
        render(fcols[1], &format!("FP {fp}"), theme.matrix_fp);
        render(fcols[2], &format!("TN {tn}"), theme.matrix_tn);
    }
}

fn draw_class_matrix(f:&mut ratatui::Frame, area:Rect, model:&Model, theme:&PastelTheme){
    
    // Panels layout 2 rows x 4 cols
    let panel_rows=Layout::default().direction(Direction::Vertical)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(Rect{ x:area.x+1, y:area.y+1, width:area.width.saturating_sub(2), height:area.height.saturating_sub(2) });

    if panel_rows.len()<2 { return; }

    let top_panels = Layout::default().direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(25),Constraint::Percentage(25),
                      Constraint::Percentage(25),Constraint::Percentage(25)])
        .split(panel_rows[0]);

    let bottom_panels = Layout::default().direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(25),Constraint::Percentage(25),
                      Constraint::Percentage(25),Constraint::Percentage(25)])
        .split(panel_rows[1]);

    let mut panel_rects = Vec::new();
    panel_rects.extend_from_slice(&top_panels);
    panel_rects.extend_from_slice(&bottom_panels);

    // Each panel gets 12 samples
    for (pi, prect) in panel_rects.into_iter().enumerate() {
        draw_class_panel(f, prect, model, theme, pi);
    }
}

fn draw_class_panel(f:&mut ratatui::Frame, area:Rect, model:&Model, theme:&PastelTheme, panel_idx:usize){
    
    // Panel border (thin)
    f.render_widget(Block::default().borders(Borders::NONE), area);

    if area.width < 6 || area.height < 4 { return; }

    let start = panel_idx * 12;
    let end = ((panel_idx + 1)*12).min(96);

    if start >= end { return; }

    let inner = Rect{ x:area.x+1, y:area.y+1, width:area.width.saturating_sub(2), height:area.height.saturating_sub(2) };

    let replicates = model.replicates.min(10); // limit width

    let mut lines: Vec<Line> = Vec::new();
    for sample in start..end {
        let consensus = model.sample_classes[sample];
        let consensus_style = style_for_cell(consensus, theme).patch(theme.matrix_consensus);
        let mut spans = vec![
            Span::styled(con_char(consensus), consensus_style),
            Span::raw(" "),
        ];
        for r in 0..replicates {
            let cell = model.replicate_cell(r, sample);
            spans.push(Span::styled(cell.glyph(), style_for_cell(cell, theme)));
        }
        lines.push(Line::from(spans));
    }

    f.render_widget(Paragraph::new(Text::from(lines)), inner);
}

fn con_char(c:CellClass)->&'static str {
    match c {
        CellClass::TP => "■",
        CellClass::TN => "■",
        CellClass::FP => "■",
        CellClass::FN => "■",
        CellClass::Excluded => "□",
        CellClass::Control => "◆",
    }
}

fn style_for_cell(c:CellClass, theme:&PastelTheme)->Style {
    match c {
        CellClass::TP => theme.matrix_tp,
        CellClass::TN => theme.matrix_tn,
        CellClass::FP => theme.matrix_fp,
        CellClass::FN => theme.matrix_fn,
        CellClass::Excluded => theme.matrix_ex,
        CellClass::Control => theme.matrix_ctrl,
    }
}

fn draw_detail(f:&mut ratatui::Frame, area:Rect, app:&App, theme:&PastelTheme){
    
    let chunks=Layout::default().direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(50),
            Constraint::Percentage(50),
        ]).split(area);

    f.render_widget(Block::default().borders(Borders::ALL)
        .title(Span::styled("Detail", theme.header_detail)), area);

    // Confusion

    let conf_inner=Rect{ x:chunks[1].x+1, y:chunks[1].y+1, width:chunks[1].width.saturating_sub(2), height:chunks[1].height.saturating_sub(2) };
    draw_confusion_matrix(f, conf_inner, app.selected_model(), theme);

    // Text summary
    let m=app.selected_model();
    let (gmin,gmed,gmax)=m.gpu_stats();
    let (rmin,rmed,rmax)=m.run_stats();
    let sens=m.sensitivity()*100.0;
    let spec=m.specificity()*100.0;
    let text=Text::from(vec![
        Line::from(format!("Name: {}", m.name)),
        Line::from(""),
        Line::from(format!("Sensitivity: {:.2}%", sens,)),
        Line::from(format!("Specificity: {:.2}%", spec,)),
        Line::from(""),
        Line::from(format!("GPU {gmin}/{gmed:.1}/{gmax}  RUN {rmin}/{rmed:.1}/{rmax}")),
        Line::from(""),
        Line::from("Keys: ↑/↓ select | m theme | q quit"),
    ]);

    let text_inner=Rect{ x:chunks[0].x+1, y:chunks[0].y+1, width:chunks[0].width.saturating_sub(2), height:chunks[0].height.saturating_sub(2) };
    f.render_widget(Paragraph::new(text), text_inner);

}

fn draw_footer(_f:&mut ratatui::Frame, _size:Rect, _theme:&PastelTheme){
}
