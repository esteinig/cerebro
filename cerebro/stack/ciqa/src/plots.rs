
use plotters::{prelude::*, style::text_anchor::{HPos, Pos, VPos}};
use std::{f64::consts::PI, num::ParseIntError};

const CATS: [&str; 5] = ["RPC", "NPV", "Specificity", "Sensitivity", "PPV"];

fn clamp01(x: f64) -> f64 {
    if x < 0.0 { 0.0 } else if x > 1.0 { 1.0 } else { x }
}

fn pt(cx: i32, cy: i32, r: f64, theta: f64) -> (i32, i32) {
    let x = cx as f64 + r * theta.cos();
    let y = cy as f64 - r * theta.sin(); // screen y grows down
    (x.round() as i32, y.round() as i32)
}
   
   
pub fn draw_radar_chart<B: DrawingBackend>(
    backend: B,
    series: &[(String, [f64; 5])],
    colors: Option<Vec<String>>,
    size: (u32, u32),
    extra_level: Option<f64>, 
    gridweight: GridWeight,
    labelsize: LabelSize,
    lineweight: Option<u32>
) -> Result<(), DrawingAreaErrorKind<B::ErrorType>> {
    let (w, h) = size;
    let root = backend.into_drawing_area();
    root.fill(&WHITE)?;

    // Layout
    let margin = 90i32;
    let cx = (w as i32) / 2;
    let cy = (h as i32) / 2 + 10;
    let radius = ((w.min(h) as i32) / 2 - margin) as f64;

    // Angles
    let k = CATS.len();
    let angles: Vec<f64> = (0..k).map(|i| 2.0 * PI * (i as f64) / (k as f64) - PI / 2.0).collect();

    let (grid_rgba, grid_stroke, label_rgb, grid_label_rgba) = grid_look(gridweight);
    let s = label_scale(labelsize);

    // fonts; if s==0.0 labels are off
    let cat_px   = (16.0 * s).round() as i32;
    let tick_px  = (12.0 * s).round() as i32;
    let legend_px= (14.0 * s).round() as i32;

    // Grid
    if !matches!(gridweight, GridWeight::Off) {
        let levels = 5;
        for lvl in 1..=levels {
            let r = radius * (lvl as f64) / (levels as f64);
            let mut poly: Vec<(i32, i32)> = angles.iter().map(|&ang| pt(cx, cy, r, ang)).collect();
            poly.push(poly[0]);
            root.draw(&PathElement::new(
                poly,
                ShapeStyle { color: grid_rgba, filled: false, stroke_width: grid_stroke },
            ))?;
        }
    }

    // Spokes + category labels
    if !matches!(gridweight, GridWeight::Off) {
        let axis_style = ShapeStyle { color: grid_rgba, filled: false, stroke_width: grid_stroke.max(1) };
        for (i, &ang) in angles.iter().enumerate() {
            let end = pt(cx, cy, radius, ang);
            root.draw(&PathElement::new(vec![(cx, cy), end], axis_style.clone()))?;
        }
    }

    if s > 0.0 {
        let label_font = ("monospace", cat_px).into_font().style(FontStyle::Bold).color(&label_rgb);
        for (i, &ang) in angles.iter().enumerate() {
            let (lx, ly) = pt(cx, cy, radius + 22.0, ang);
            let hpos = if ang.cos() > 0.2 { HPos::Left } else if ang.cos() < -0.2 { HPos::Right } else { HPos::Center };
            let vpos = if ang.sin() > 0.2 { VPos::Bottom } else if ang.sin() < -0.2 { VPos::Top } else { VPos::Center };
            let style = TextStyle::from(label_font.clone()).pos(Pos::new(hpos, vpos));
            root.draw(&Text::new(CATS[i].to_string(), (lx, ly), style))?;
        }
    }

    // Radial ticks
    if s > 0.0 && !matches!(gridweight, GridWeight::Off) {
        let tick_font = ("monospace", tick_px).into_font().color(&grid_label_rgba);
        let levels = 5;
        for lvl in 1..=levels {
            let r = radius * (lvl as f64) / (levels as f64);
            let p = pt(cx, cy, r, -PI / 2.0);
            root.draw(&Text::new(
                format!("{:.0}%", (lvl as f64 / levels as f64)*100.0),
                (p.0 + 6, p.1),
                tick_font.clone().pos(Pos::new(HPos::Left, VPos::Center)),
            ))?;
        }
    }

    // Pre-pick colors so we can draw legend after the polygons
    let mut picked: Vec<(RGBAColor, RGBAColor)> = Vec::with_capacity(series.len());
    for i in 0..series.len() {
        let (stroke_rgba, fill_rgba) = if let Some(ref cs) = colors {
            if let Some(hex) = cs.get(i) {

                if let Some((r,g,b, mut a_opt)) = parse_hex_rgba(hex) {
                    let mut stroke_rgba = RGBAColor(r,g,b,1.0);
                    if hex == "#ffffff" {
                        a_opt = Some(0.0);
                        stroke_rgba = RGBAColor(0,0,0,1.0);
                    }
                    (stroke_rgba, RGBAColor(r,g,b, a_opt.unwrap_or(0.5)))
                } else {
                    let base = Palette99::pick(i);
                    let (r0,g0,b0) = base.rgb();
                    (RGBAColor(r0,g0,b0,1.0), RGBAColor(r0,g0,b0,0.5))
                }
            } else {
                let base = Palette99::pick(i);
                let (r0,g0,b0) = base.rgb();
                (RGBAColor(r0,g0,b0,1.0), RGBAColor(r0,g0,b0,0.5))
            }
        } else {
            let base = Palette99::pick(i);
            let (r0,g0,b0) = base.rgb();
            (RGBAColor(r0,g0,b0,1.0), RGBAColor(r0,g0,b0,0.5))
        };
        picked.push((stroke_rgba, fill_rgba));
    }
    
    // Series
    for (i, (_, vals)) in series.iter().enumerate() {
        
        let (stroke_rgba, fill_rgba) = picked[i].clone();

        let mut pts: Vec<(i32, i32)> = angles
            .iter()
            .enumerate()
            .map(|(j, &ang)| pt(cx, cy, radius * clamp01(vals[j]), ang))
            .collect();
        pts.push(pts[0]);

        root.draw(&Polygon::new(pts.clone(), fill_rgba.filled()))?;
        root.draw(&PathElement::new(
            pts,
            ShapeStyle { color: stroke_rgba, filled: false, stroke_width: lineweight.unwrap_or(2) },
        ))?;
    }


    // legend swatch + layout scale with s
    let sw_w = (18.0 * s).round() as i32;
    let sw_h = (14.0 * s).round() as i32;
    let gap_sw_text = (8.0 * s).round() as i32;
    let gap_items   = (16.0 * s).round() as i32;
    let char_px     = (0.50 * legend_px as f64).round() as i32; // monospace ~0.5em per char
    
    if s > 0.0 {
        // widths using scaled char_px
        let item_widths: Vec<i32> = series.iter()
            .map(|(name, _)| sw_w + gap_sw_text + (name.chars().count() as i32) * char_px)
            .collect();

        let total_w: i32 = item_widths.iter().copied().sum::<i32>() 
                         + gap_items * ((series.len() as i32).saturating_sub(1));

        let legend_y = (cy as f64 - radius - 28.0 * s) as i32;
        let baseline_y = legend_y;
        let legend_x0 = cx - total_w / 2;

        let text_style = ("monospace", legend_px)
            .into_font()
            .color(&label_rgb)
            .into_text_style(&root)
            .pos(Pos::new(HPos::Left, VPos::Center));

        let mut x = legend_x0;
        for (i, (name, _)) in series.iter().enumerate() {
            let fill = picked[i].1;
        
            let left   = x;
            let top    = baseline_y - sw_h / 2;
            let right  = x + sw_w;
            let bottom = baseline_y + sw_h / 2;
        
            // fill
            root.draw(&Rectangle::new(
                [(left, top), (right, bottom)],
                fill.filled(),
            ))?;
        
            // black border
            let stroke = (1.0 * s).max(1.0).round() as u32;
            root.draw(&Rectangle::new(
                [(left, top), (right, bottom)],
                ShapeStyle {
                    color: picked[i].0,
                    filled: false,
                    stroke_width: stroke,
                },
            ))?;
        
            x = right + gap_sw_text;
        
            root.draw(&Text::new(name.clone(), (x, baseline_y), text_style.clone()))?;
            x += (name.chars().count() as i32) * char_px;
            if i + 1 != series.len() { x += gap_items; }
        }
    }

    // Optional extra level ring on top as reference line)
    if let Some(v) = extra_level {
        if (0.0..=1.0).contains(&v) {
            let r = radius * v.max(0.0).min(1.0);
            let mut poly: Vec<(i32, i32)> = angles.iter().map(|&ang| pt(cx, cy, r, ang)).collect();
            poly.push(poly[0]);

            // More opaque + dashed
            let style = ShapeStyle {
                color: RGBAColor(240, 187, 77, 1.0).to_rgba(),
                filled: false,
                stroke_width: 2,
            };

            root.draw(&PathElement::new(poly, style))?;
        }
    }

    root.present()?;
    Ok(())
}


fn parse_hex_rgba(s: &str) -> Option<(u8,u8,u8,Option<f64>)> {
    let h = s.strip_prefix('#').unwrap_or(s);
    let hex = h.trim();
    let to_u8 = |i| u8::from_str_radix(&hex[i..i+2], 16);
    let ok = |r: Result<u8, ParseIntError>| r.ok();

    match hex.len() {
        6 => {
            let (r,g,b) = (ok(to_u8(0))?, ok(to_u8(2))?, ok(to_u8(4))?);
            Some((r,g,b,None))
        }
        8 => {
            let (r,g,b,a) = (ok(to_u8(0))?, ok(to_u8(2))?, ok(to_u8(4))?, ok(to_u8(6))?);
            Some((r,g,b,Some(a as f64 / 255.0)))
        }
        _ => None,
    }
}

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum GridWeight { Bold, Normal, Light, Off }

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum LabelSize { ExtraLarge, Large, Normal, Small, Off }

fn grid_look(w: GridWeight) -> (RGBAColor, u32, RGBColor, RGBColor) {
    match w {
        GridWeight::Bold   => (RGBAColor(140,140,140,1.0), 2, RGBColor(30,30,30), RGBColor(120,120,120)),
        GridWeight::Normal => (RGBAColor(180,180,180,1.0), 1, RGBColor(60,60,60), RGBColor(160,160,160)),
        GridWeight::Light  => (RGBAColor(224,224,224,1.0), 1, RGBColor(90,90,90), RGBColor(200,200,200)),
        GridWeight::Off    => (RGBAColor(0,0,0,0.0),       0, RGBColor(0,0,0), RGBColor(0,0,0)), // unused
    }
}

fn label_scale(s: LabelSize) -> f64 {
    match s {
        LabelSize::ExtraLarge => 2.0,
        LabelSize::Large      => 1.75,
        LabelSize::Normal     => 1.5,
        LabelSize::Small      => 1.0,
        LabelSize::Off        => 0.0,
    }
}