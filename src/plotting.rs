// Plotting module for length distribution visualization (SVG only)
use plotters::prelude::*;
use std::error::Error;

/// Plot length distribution histogram with cumulative curve (SVG output)
pub fn plot_length_distribution(
    lengths: &[i64],
    output_path: &str,
    sample_name: &str,
) -> Result<(), Box<dyn Error>> {
    if lengths.is_empty() {
        return Err("No length data to plot".into());
    }

    // Convert lengths to kb
    let lengths_kb: Vec<f64> = lengths.iter().map(|&x| x as f64 / 1000.0).collect();

    let min_len_kb = lengths_kb.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_len_kb = lengths_kb.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Use Sturges' rule with 1.5x multiplier for histogram bins
    let n = lengths.len() as f64;
    let n_bins = ((n.log2() + 1.0) * 1.5).ceil() as usize;
    let n_bins = n_bins.max(10).min(30); // 10-30 bins

    // Calculate bin width
    let range_kb = max_len_kb - min_len_kb;
    let bin_width = if range_kb > 0.0 {
        range_kb / n_bins as f64
    } else {
        1.0
    };

    // Build histogram bins
    let mut bins: Vec<u32> = vec![0; n_bins];
    for &len_kb in &lengths_kb {
        let bin_idx = ((len_kb - min_len_kb) / bin_width) as usize;
        let bin_idx = bin_idx.min(n_bins - 1);
        bins[bin_idx] += 1;
    }

    let max_count = *bins.iter().max().unwrap_or(&1);

    // Calculate cumulative bases (in Mb) - sorted descending by length
    let mut sorted_lengths = lengths.to_vec();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a)); // descending

    let total_bases_mb: f64 = sorted_lengths.iter().sum::<i64>() as f64 / 1_000_000.0;

    // Calculate cumulative curve points (sample for performance)
    let mut cumsum: f64 = 0.0;
    let mut cumulative_points: Vec<(f64, f64)> = Vec::new();
    let sample_interval = (sorted_lengths.len() / 1000).max(1);

    for (i, &len) in sorted_lengths.iter().enumerate() {
        cumsum += len as f64 / 1_000_000.0;
        // Sample points: always include first, last, and every nth point
        if i == 0 || i == sorted_lengths.len() - 1 || i % sample_interval == 0 {
            let len_kb = len as f64 / 1000.0;
            cumulative_points.push((len_kb, cumsum));
        }
    }

    // Sort by x for proper line drawing
    cumulative_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Draw SVG
    let root = SVGBackend::new(output_path, (400, 300)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_len_kb = min_len_kb + bin_width * n_bins as f64;

    // X axis range: from 0 to actual max length
    let x_range = 0f64..max_len_kb;

    // Y axis: round up to nice number
    let y_max = if max_count >= 1000 {
        ((max_count as f64 / 1000.0).ceil() * 1000.0) as f64
    } else {
        ((max_count as f64 / 100.0).ceil() * 100.0) as f64
    };
    let y_max = y_max.max(100.0);

    // Create chart with dual Y axes
    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .margin_right(5)
        .margin_bottom(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .right_y_label_area_size(45)
        .caption(sample_name, ("sans-serif", 16).into_font().style(FontStyle::Bold))
        .build_cartesian_2d(x_range.clone(), 0f64..y_max)?
        .set_secondary_coord(x_range.clone(), 0f64..total_bases_mb * 1.1);

    // Configure primary mesh (left Y axis - Number of Reads)
    chart.configure_mesh()
        .x_label_style(("sans-serif", 14))
        .y_label_style(("sans-serif", 14))
        .x_desc("Read Length (kb)")
        .y_desc("Number of Reads")
        .axis_desc_style(("sans-serif", 15))
        .x_labels(8)
        .y_labels(8)
        .x_label_offset(0)
        .disable_mesh()
        .y_label_formatter(&|y| {
            if *y >= 1000.0 {
                format!("{}k", (*y / 1000.0) as i32)
            } else {
                format!("{}", *y as i32)
            }
        })
        .draw()?;

    // Configure secondary mesh (right Y axis - Cumulative Bases)
    chart.configure_secondary_axes()
        .y_desc("Cumulative Bases (Mb)")
        .label_style(("sans-serif", 14))
        .axis_desc_style(("sans-serif", 15))
        .y_labels(8)
        .y_label_formatter(&|y| format!("{:.0}", y))
        .draw()?;

    // Define light green color for histogram
    let bar_color = RGBColor(144, 238, 144); // lightgreen
    let bar_border = RGBColor(34, 139, 34);  // forestgreen

    // Draw histogram bars
    for (i, &count) in bins.iter().enumerate() {
        if count == 0 {
            continue;
        }
        let x0 = min_len_kb + i as f64 * bin_width;
        let x1 = x0 + bin_width;
        let count_f64 = count as f64;

        // Filled bar
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x0, 0f64), (x1, count_f64)],
            bar_color.filled(),
        )))?;

        // Bar border
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x0, 0f64), (x1, count_f64)],
            bar_border.stroke_width(1),
        )))?;
    }

    // Draw cumulative curve on secondary axis
    chart.draw_secondary_series(LineSeries::new(
        cumulative_points.iter().cloned(),
        BLACK.stroke_width(1),
    ))?;

    root.present()?;
    Ok(())
}
