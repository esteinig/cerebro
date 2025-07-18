
use std::{
    io::{self, Stdout},
    time::{Duration, Instant},
};
use clap::{Parser, Subcommand};
use ratatui::crossterm::{
    event::{self, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::Span,
    widgets::{BarChart, Block, Borders, Paragraph},
    Terminal,
};

/// Application state holding our two series of timings.
struct App {
    gpu_times: Vec<(&'static str, u64)>,
    run_times: Vec<(&'static str, u64)>,
}

impl App {
    fn new() -> Self {
        // Dummy data: replace with real measurements later
        let gpu_times = vec![("rep1", 120), ("rep2", 150), ("rep3", 90), ("rep4", 200)];
        let run_times = vec![("rep1", 30), ("rep2", 45), ("rep3", 25), ("rep4", 60)];
        App { gpu_times, run_times }
    }
}

pub fn start_tui() -> io::Result<()> {
    // Prepare app state
    let app = App::new();

    // Terminal init
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal: Terminal<CrosstermBackend<Stdout>> = Terminal::new(backend)?;

    // Event loop
    let tick_rate = Duration::from_millis(250);
    let mut last_tick = Instant::now();
    loop {
        terminal.draw(|f| ui(f, &app))?;

        let timeout = tick_rate.checked_sub(last_tick.elapsed()).unwrap_or_default();
        if event::poll(timeout)? {
            if let Event::Key(key) = event::read()? {
                if key.code == KeyCode::Char('q') {
                    break;
                }
            }
        }
        if last_tick.elapsed() >= tick_rate {
            // could update `app` here if it were mutable
            last_tick = Instant::now();
        }
    }

    // Restore terminal
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;
    Ok(())
}

fn ui(f: &mut ratatui::Frame, app: &App) {
    let size = f.size();

    // Predefine pastel styles for headers
    let pastel_pink   = Style::default().fg(Color::Rgb(255,182,193)).add_modifier(Modifier::BOLD);
    let pastel_mint   = Style::default().fg(Color::Rgb(152,251,152)).add_modifier(Modifier::BOLD);
    let pastel_lav    = Style::default().fg(Color::Rgb(230,230,250)).add_modifier(Modifier::BOLD);
    let pastel_peach  = Style::default().fg(Color::Rgb(255,218,185)).add_modifier(Modifier::BOLD);
    let pastel_yellow = Style::default().fg(Color::Rgb(255,255,224)).add_modifier(Modifier::BOLD);

    // Split: left 70%, right 30%
    let main_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(70), Constraint::Percentage(30)])
        .split(size);

    // On the left: Controls, Metrics, Stats
    let left_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(7),  // Controls
            Constraint::Length(16), // Metrics (we’ll subdivide further)
            Constraint::Length(7),  // Stats
            Constraint::Min(0),
        ])
        .split(main_chunks[0]);

    // Controls
    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .title(Span::styled("Controls", pastel_pink)),
        left_chunks[0],
    );

    // Metrics (split into 3: two charts + summary paragraph)
    let metrics_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage(40),
            Constraint::Percentage(40),
            Constraint::Percentage(20),
        ])
        .split(left_chunks[1]);

    // 1) GPU time bar chart
    let gpu_chart = BarChart::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(Span::styled("GPU Time (ms)", pastel_mint)),
        )
        .data(&app.gpu_times)
        .bar_width(5)
        .bar_gap(2)
        .bar_style(Style::default().add_modifier(Modifier::ITALIC))
        .value_style(Style::default());
    f.render_widget(gpu_chart, metrics_chunks[0]);

    // 2) Run time bar chart
    let run_chart = BarChart::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(Span::styled("Run Time (ms)", pastel_lav)),
        )
        .data(&app.run_times)
        .bar_width(5)
        .bar_gap(2)
        .bar_style(Style::default().add_modifier(Modifier::ITALIC))
        .value_style(Style::default());
    f.render_widget(run_chart, metrics_chunks[1]);

    // 3) Summary paragraph with min/median/max
    // Compute stats
    let mut gvs: Vec<u64> = app.gpu_times.iter().map(|(_,v)| *v).collect();
    let mut rvs: Vec<u64> = app.run_times.iter().map(|(_,v)| *v).collect();
    gvs.sort_unstable();
    rvs.sort_unstable();
    let (gmin, gmax) = (*gvs.first().unwrap(), *gvs.last().unwrap());
    let (rmin, rmax) = (*rvs.first().unwrap(), *rvs.last().unwrap());
    let gmed = (gvs[gvs.len()/2 - 1] + gvs[gvs.len()/2]) as f64 / 2.0;
    let rmed = (rvs[rvs.len()/2 - 1] + rvs[rvs.len()/2]) as f64 / 2.0;
    let summary = format!(
        "GPU: min {gmin} ms │ median {gmed:.1} ms │ max {gmax} ms\n\
         Run: min {rmin} ms │ median {rmed:.1} ms │ max {rmax} ms"
    );

    let para = Paragraph::new(summary)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(Span::styled("Summary", pastel_peach)),
        );
    f.render_widget(para, metrics_chunks[2]);

    // Stats
    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .title(Span::styled("Stats", pastel_yellow)),
        left_chunks[2],
    );

    // Right detail pane
    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .title("Detail"),
        main_chunks[1],
    );

    // Footer across bottom
    let footer = Rect {
        x: size.x,
        y: size.height.saturating_sub(3),
        width: size.width,
        height: 3,
    };
    f.render_widget(
        Block::default()
            .borders(Borders::ALL)
            .title("Footer"),
        footer,
    );
}