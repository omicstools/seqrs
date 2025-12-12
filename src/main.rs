use bio::seq_analysis::gc::gc_content;
use clap::{Parser, Subcommand};
use indicatif::{HumanCount, ProgressBar};
use kseq::parse_path;
use modules::write_fastq;
use rand::Rng;
use regex::{bytes::RegexSet, Regex};
use std::{fs, io::BufRead, io::BufReader, io::BufWriter, process, time::Duration};
use ferris_says::say;

mod modules;
mod plotting;

const VERSION: &str = "2.0.0";

fn print_ferris_banner() {
    let message = format!(
        "seqrs v{} - Fast statistics and filtering for FASTQ/FASTA files.
Features: QC reports, length distribution plots, filtering, trimming, and more!",
        VERSION
    );
    let width = 70;
    let mut writer = BufWriter::new(std::io::stderr());
    let _ = say(&message, width, &mut writer);
}

// ============================================================================
// CLI 结构定义
// ============================================================================

#[derive(Parser)]
#[command(
    name = "seqrs",
    version = VERSION,
    author = "OMICS Tool Development Team",
    about = "Fast statistics and filtering tool for FASTQ/FASTA files",
    after_help = "EXAMPLES:
    # Statistics
    seqrs stat input.fastq                  # basic statistics
    seqrs stat -t input.fastq               # table output
    seqrs stat --qc sample1 input.fastq     # QC report
    seqrs stat -l input.fastq               # output lengths
    seqrs stat -x 0.5 input.fastq           # calculate N50

    # Summarize QC reports
    seqrs stat --sum-dir ./qc_reports                 # merge .qc_stat.txt files
    seqrs stat --sum-dir ./qc_reports -o summary.tsv  # output to TSV file
    seqrs stat --sum-dir ./qc_reports -o summary.csv  # output to CSV file

    # Filter reads
    seqrs filter -l 1000 input.fastq        # keep reads >= 1000bp
    seqrs filter -L 500 input.fastq         # keep reads <= 500bp
    seqrs filter -q 20 input.fastq          # keep reads with Q >= 20
    seqrs filter -r 'pattern' input.fastq   # filter by ID regex

    # Trim reads
    seqrs trim -f 50 input.fastq            # trim 50bp from front
    seqrs trim -t 30 input.fastq            # trim 30bp from tail

    # Sample reads
    seqrs sample -p 0.1 input.fastq         # 10% sub-sample
    seqrs sample -n 1000 input.fastq        # sample 1000 reads
"
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Statistics and QC reports
    Stat(StatArgs),
    /// Filter reads by length, quality, or ID
    Filter(FilterArgs),
    /// Trim reads from front or tail
    Trim(TrimArgs),
    /// Sub-sample reads
    Sample(SampleArgs),
}

// ============================================================================
// Stat 子命令
// ============================================================================

#[derive(Parser)]
#[command(
    about = "Statistics and QC reports",
    after_help = "EXAMPLES:
    seqrs stat input.fastq                  # basic statistics
    seqrs stat -t input.fastq               # table output
    seqrs stat -ts input.fastq              # table without header
    seqrs stat --qc sample1 input.fastq     # QC report with plots
    seqrs stat -l input.fastq               # output per-read lengths
    seqrs stat -g input.fastq               # output per-read GC content
    seqrs stat -q input.fastq               # output per-read quality scores
    seqrs stat -x 0.5 input.fastq           # calculate N50
    seqrs stat -y 20 input.fastq            # calculate Q20 yield

    # Summarize QC reports
    seqrs stat --sum-dir ./qc_reports                 # merge .qc_stat.txt to stdout
    seqrs stat --sum-dir ./qc_reports -o summary.tsv  # output as TSV
    seqrs stat --sum-dir ./qc_reports -o summary.csv  # output as CSV (auto-detect)
"
)]
struct StatArgs {
    /// Input FASTQ/FASTA file (supports gzip)
    #[arg(required_unless_present = "sum_dir")]
    input: Option<String>,

    // === Summarize Mode ===
    /// Summarize .qc_stat.txt files from directory
    #[arg(long = "sum-dir", value_name = "DIR", conflicts_with_all = ["lengths", "gc", "quality", "nx", "qyield", "qc", "table"])]
    sum_dir: Option<String>,

    // === Output Format ===
    /// Output tab-separated statistics table
    #[arg(short = 't', long = "table")]
    table: bool,

    /// Skip header line in table output
    #[arg(short = 's', long = "skip-header")]
    skip_header: bool,

    /// Generate QC report with sample name (creates .qc_stat.txt and .length_dist.svg)
    #[arg(long = "qc", value_name = "NAME")]
    qc: Option<String>,

    /// Output directory (for --qc) or output file (for --sum-dir)
    #[arg(short = 'o', long = "output", value_name = "PATH")]
    output: Option<String>,

    // === Per-read Output ===
    /// Output read lengths (one per line)
    #[arg(short = 'l', long = "lengths")]
    lengths: bool,

    /// Output GC content per read (one per line)
    #[arg(short = 'g', long = "gc")]
    gc: bool,

    /// Output mean Phred quality score per read (one per line)
    #[arg(short = 'q', long = "quality")]
    quality: bool,

    // === Metric Calculation ===
    /// Calculate NX value (e.g., 0.5 for N50, 0.9 for N90)
    #[arg(short = 'x', long = "nx", value_name = "VALUE")]
    nx: Option<f32>,

    /// Calculate percent bases >= Q score [8-60]
    #[arg(short = 'y', long = "qyield", value_name = "Q")]
    qyield: Option<u8>,
}

// ============================================================================
// Filter 子命令
// ============================================================================

#[derive(Parser)]
#[command(
    about = "Filter reads by length, quality, or ID",
    after_help = "EXAMPLES:
    seqrs filter -l 1000 input.fastq          # keep reads >= 1000bp
    seqrs filter -L 500 input.fastq           # keep reads <= 500bp
    seqrs filter -l 1000 -L 5000 input.fastq  # keep reads 1000-5000bp
    seqrs filter -q 20 input.fastq            # keep reads with Q >= 20
    seqrs filter -Q 30 input.fastq            # keep reads with Q <= 30
    seqrs filter -r 'chr1' input.fastq        # filter by ID regex
    seqrs filter -R patterns.txt input.fastq  # filter by patterns from file
    seqrs filter -r 'chr1' -v input.fastq     # exclude matching reads
"
)]
struct FilterArgs {
    /// Input FASTQ/FASTA file (supports gzip)
    #[arg(required = true)]
    input: String,

    // === Length Filter ===
    /// Minimum length filter (keep reads >= N bp)
    #[arg(short = 'l', long = "min-len", value_name = "N")]
    min_len: Option<i32>,

    /// Maximum length filter (keep reads <= N bp)
    #[arg(short = 'L', long = "max-len", value_name = "N")]
    max_len: Option<i32>,

    // === Quality Filter ===
    /// Minimum mean quality filter (keep reads with Q >= N) [8-60]
    #[arg(short = 'q', long = "min-qual", value_name = "Q")]
    min_qual: Option<i32>,

    /// Maximum mean quality filter (keep reads with Q <= N) [8-60]
    #[arg(short = 'Q', long = "max-qual", value_name = "Q")]
    max_qual: Option<i32>,

    // === ID Filter ===
    /// Filter reads by ID matching regex pattern
    #[arg(short = 'r', long = "regex", value_name = "PATTERN")]
    regex: Option<String>,

    /// Filter reads by ID matching patterns from file (one per line)
    #[arg(short = 'R', long = "regex-file", value_name = "FILE")]
    regex_file: Option<String>,

    /// Invert match (exclude matching reads)
    #[arg(short = 'v', long = "invert")]
    invert: bool,
}

// ============================================================================
// Trim 子命令
// ============================================================================

#[derive(Parser)]
#[command(
    about = "Trim reads from front or tail",
    after_help = "EXAMPLES:
    seqrs trim -f 50 input.fastq              # trim 50bp from 5' end
    seqrs trim -t 30 input.fastq              # trim 30bp from 3' end
    seqrs trim -f 50 -t 30 input.fastq        # trim both ends
    seqrs trim -f 50 -l 100 input.fastq       # trim and keep reads >= 100bp
"
)]
struct TrimArgs {
    /// Input FASTQ/FASTA file (supports gzip)
    #[arg(required = true)]
    input: String,

    /// Trim N bases from the 5' end (front)
    #[arg(short = 'f', long = "front", value_name = "N")]
    front: Option<usize>,

    /// Trim N bases from the 3' end (tail)
    #[arg(short = 't', long = "tail", value_name = "N")]
    tail: Option<usize>,

    /// Minimum length after trimming (discard shorter reads)
    #[arg(short = 'l', long = "min-len", value_name = "N")]
    min_len: Option<usize>,
}

// ============================================================================
// Sample 子命令
// ============================================================================

#[derive(Parser)]
#[command(
    about = "Sub-sample reads",
    after_help = "EXAMPLES:
    seqrs sample -p 0.1 input.fastq           # 10% sub-sample
    seqrs sample -n 1000 input.fastq          # sample 1000 reads
    seqrs sample -p 0.1 -s 42 input.fastq     # reproducible sampling
"
)]
struct SampleArgs {
    /// Input FASTQ/FASTA file (supports gzip)
    #[arg(required = true)]
    input: String,

    /// Sample by proportion [0.0-1.0]
    #[arg(short = 'p', long = "proportion", value_name = "FRAC")]
    proportion: Option<f32>,

    /// Sample by number of reads
    #[arg(short = 'n', long = "number", value_name = "N")]
    number: Option<usize>,

    /// Random seed for reproducible sampling
    #[arg(short = 's', long = "seed", value_name = "SEED")]
    seed: Option<u64>,
}

// ============================================================================
// 主函数
// ============================================================================

fn main() {
    // Check if -h or --help is in args
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print_ferris_banner();
        eprintln!();
    }

    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Stat(args)) => run_stat(args),
        Some(Commands::Filter(args)) => run_filter(args),
        Some(Commands::Trim(args)) => run_trim(args),
        Some(Commands::Sample(args)) => run_sample(args),
        None => {
            // 没有子命令时显示帮助
            print_ferris_banner();
            eprintln!();
            eprintln!("Usage: seqrs <COMMAND> [OPTIONS] <INPUT>");
            eprintln!();
            eprintln!("Commands:");
            eprintln!("  stat    Statistics and QC reports");
            eprintln!("  filter  Filter reads by length, quality, or ID");
            eprintln!("  trim    Trim reads from front or tail");
            eprintln!("  sample  Sub-sample reads");
            eprintln!();
            eprintln!("Run 'seqrs <COMMAND> --help' for more information.");
            process::exit(0);
        }
    }
}

// ============================================================================
// Stat 命令实现
// ============================================================================

fn run_stat(args: StatArgs) {
    // === Summarize mode ===
    if let Some(sum_dir) = &args.sum_dir {
        run_stat_sum(sum_dir, &args.output);
        return;
    }

    // === Normal stat mode - requires input file ===
    let input = args.input.as_ref().expect("Input file is required");
    let mut records = parse_path(input).unwrap();

    // Per-read output: lengths
    if args.lengths {
        while let Some(record) = records.iter_record().unwrap() {
            println!("{}", record.len());
        }
        process::exit(0);
    }

    // Per-read output: GC content
    if args.gc {
        while let Some(record) = records.iter_record().unwrap() {
            let seq = record.seq().as_bytes();
            println!("{}", gc_content(seq));
        }
        process::exit(0);
    }

    // Per-read output: quality scores
    if args.quality {
        while let Some(record) = records.iter_record().unwrap() {
            let mean_errorp = modules::qscore_probs(record.qual().as_bytes()) / record.seq().len() as f32;
            println!("{:.4}", -10.0 * mean_errorp.log10());
        }
        process::exit(0);
    }

    // Calculate NX value
    if let Some(nxvalue) = args.nx {
        if !(0.0..=1.0).contains(&nxvalue) {
            eprintln!("Error: NX value should be between 0.0 and 1.0");
            process::exit(1);
        }
        let mut lengths: Vec<i64> = Vec::new();
        while let Some(record) = records.iter_record().unwrap() {
            lengths.push(record.seq().len() as i64);
        }
        let nx = modules::get_nx(&mut lengths, 1.0 - nxvalue);
        let nx_label = (nxvalue * 100.0) as i16;
        println!("N{}\t{}", nx_label, nx);
        process::exit(0);
    }

    // Calculate Q yield
    if let Some(qvalue) = args.qyield {
        if !(8..=60).contains(&qvalue) {
            eprintln!("Error: Q value should be between 8 and 60");
            process::exit(1);
        }
        let mut bases: i64 = 0;
        let mut qualx: i64 = 0;
        while let Some(record) = records.iter_record().unwrap() {
            let len = record.seq().len() as i64;
            bases += len;
            qualx += modules::get_qual_bases(record.qual().as_bytes(), 33 + qvalue);
        }
        let qx = qualx as f64 / bases as f64 * 100.0;
        println!("Q{}\t{:.2}", qvalue, qx);
        process::exit(0);
    }

    // QC report generation
    if let Some(sample_name) = &args.qc {
        let outdir = args.output.as_deref().unwrap_or(".");
        let outdir_path = std::path::Path::new(outdir);
        if !outdir_path.exists() {
            std::fs::create_dir_all(outdir_path).expect("Failed to create output directory");
        }

        let stat_file = outdir_path.join(format!("{}.qc_stat.txt", sample_name));
        let plot_file = outdir_path.join(format!("{}.length_dist.svg", sample_name));

        let mut reads: i64 = 0;
        let mut bases: i64 = 0;
        let mut num_n: i32 = 0;
        let mut qual20: i64 = 0;
        let mut gc_bases: i64 = 0;
        let mut qual30: i64 = 0;
        let mut minlen: i64 = i64::MAX;
        let mut maxlen: i64 = 0;
        let mut len_vector: Vec<i64> = Vec::new();
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(Duration::from_millis(120));

        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len() as i64;
            reads += 1;
            bases += len;
            num_n += modules::get_n_bases(record.seq().as_bytes());
            qual20 += modules::get_qual_bases(record.qual().as_bytes(), 53);
            qual30 += modules::get_qual_bases(record.qual().as_bytes(), 63);
            gc_bases += modules::get_gc_bases(record.seq().as_bytes());
            minlen = len.min(minlen);
            maxlen = len.max(maxlen);
            len_vector.push(len);
            pb.set_message(format!("Processed reads: {}", HumanCount(reads as u64)));
        }

        let mean_len = modules::mean(&len_vector);
        let quart1 = modules::quartiles(&mut len_vector.clone(), 1);
        let quart2 = modules::quartiles(&mut len_vector.clone(), 2);
        let quart3 = modules::quartiles(&mut len_vector.clone(), 3);
        let n50 = modules::get_nx(&mut len_vector.clone(), 0.5);
        let q20 = qual20 as f64 / bases as f64 * 100.0;
        let q30 = qual30 as f64 / bases as f64 * 100.0;
        let gc_percent = if bases > 0 { gc_bases as f64 / bases as f64 * 100.0 } else { 0.0 };
        pb.finish_and_clear();

        let stats = modules::QcStats {
            sample_name: sample_name.to_string(),
            reads,
            bases,
            n_bases: num_n,
            min_len: minlen,
            max_len: maxlen,
            mean_len,
            q1: quart1,
            q2: quart2,
            q3: quart3,
            n50,
            q20_percent: q20,
            q30_percent: q30,
            gc_percent,
        };

        if let Err(e) = modules::write_qc_stats(&stats, &stat_file) {
            eprintln!("Failed to write QC stats: {}", e);
        }

        if let Err(e) = plotting::plot_length_distribution(&len_vector, plot_file.to_str().unwrap(), sample_name) {
            eprintln!("Failed to generate plot: {}", e);
        }

        eprintln!("QC report generated:");
        eprintln!("  Stats: {}", stat_file.display());
        eprintln!("  Plot:  {}", plot_file.display());
        process::exit(0);
    }

    // Default: table output
    let mut reads: i64 = 0;
    let mut bases: i64 = 0;
    let mut num_n: i32 = 0;
    let mut qual20: i64 = 0;
    let mut gc_bases: i64 = 0;
    let mut qual30: i64 = 0;
    let mut minlen: i64 = i64::MAX;
    let mut maxlen: i64 = 0;
    let mut len_vector: Vec<i64> = Vec::new();
    let pb = ProgressBar::new_spinner();
    pb.enable_steady_tick(Duration::from_millis(120));

    while let Some(record) = records.iter_record().unwrap() {
        let len = record.len() as i64;
        reads += 1;
        bases += len;
        num_n += modules::get_n_bases(record.seq().as_bytes());
        qual20 += modules::get_qual_bases(record.qual().as_bytes(), 53);
        qual30 += modules::get_qual_bases(record.qual().as_bytes(), 63);
        gc_bases += modules::get_gc_bases(record.seq().as_bytes());
        minlen = len.min(minlen);
        maxlen = len.max(maxlen);
        len_vector.push(len);
        pb.set_message(format!("Processed reads: {}", HumanCount(reads as u64)));
    }

    let mean_len = modules::mean(&len_vector);
    let quart1 = modules::quartiles(&mut len_vector.clone(), 1);
    let quart2 = modules::quartiles(&mut len_vector.clone(), 2);
    let quart3 = modules::quartiles(&mut len_vector.clone(), 3);
    let n50 = modules::get_nx(&mut len_vector, 0.5);
    let q20 = qual20 as f64 / bases as f64 * 100.0;
    let q30 = qual30 as f64 / bases as f64 * 100.0;
    let gc_percent = if bases > 0 { gc_bases as f64 / bases as f64 * 100.0 } else { 0.0 };
    pb.finish_and_clear();

    if !args.skip_header {
        println!("file\treads\tbases\tn_bases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent\tGC_percent");
    }
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
        input, reads, bases, num_n, minlen, maxlen, mean_len,
        quart1, quart2, quart3, n50, q20, q30, gc_percent
    );
}

/// Summarize .qc_stat.txt files from a directory
fn run_stat_sum(dir: &str, output: &Option<String>) {
    use std::io::Write;

    let dir_path = std::path::Path::new(dir);
    if !dir_path.exists() || !dir_path.is_dir() {
        eprintln!("Error: '{}' is not a valid directory", dir);
        process::exit(1);
    }

    // Find all .qc_stat.txt files
    let mut stat_files: Vec<_> = match std::fs::read_dir(dir_path) {
        Ok(entries) => entries
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().map_or(false, |ext| ext == "txt"))
            .filter(|p| p.file_name().map_or(false, |name| name.to_string_lossy().ends_with(".qc_stat.txt")))
            .collect(),
        Err(e) => {
            eprintln!("Error reading directory: {}", e);
            process::exit(1);
        }
    };

    if stat_files.is_empty() {
        eprintln!("Error: No .qc_stat.txt files found in '{}'", dir);
        process::exit(1);
    }

    // Sort files by name for consistent output
    stat_files.sort();

    // Determine format from output file extension
    let use_csv = output.as_ref().map_or(false, |p| {
        std::path::Path::new(p)
            .extension()
            .map_or(false, |ext| ext.eq_ignore_ascii_case("csv"))
    });

    // Collect all data rows
    let mut header: Option<String> = None;
    let mut data_rows: Vec<String> = Vec::new();

    for file_path in &stat_files {
        let content = match std::fs::read_to_string(file_path) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("Warning: Failed to read '{}': {}", file_path.display(), e);
                continue;
            }
        };

        let lines: Vec<&str> = content.lines().collect();
        if lines.len() < 2 {
            eprintln!("Warning: Invalid format in '{}'", file_path.display());
            continue;
        }

        // Store header from first valid file
        if header.is_none() {
            let h = if use_csv {
                lines[0].replace('\t', ",")
            } else {
                lines[0].to_string()
            };
            header = Some(h);
        }

        // Add data row (skip header line)
        for line in lines.iter().skip(1) {
            if !line.trim().is_empty() {
                let row = if use_csv {
                    line.replace('\t', ",")
                } else {
                    line.to_string()
                };
                data_rows.push(row);
            }
        }
    }

    if data_rows.is_empty() {
        eprintln!("Error: No valid data found in QC stat files");
        process::exit(1);
    }

    // Output
    let output_content = format!(
        "{}\n{}",
        header.unwrap_or_default(),
        data_rows.join("\n")
    );

    match output {
        Some(output_path) => {
            let mut file = match std::fs::File::create(output_path) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Error creating output file: {}", e);
                    process::exit(1);
                }
            };
            if let Err(e) = writeln!(file, "{}", output_content) {
                eprintln!("Error writing to file: {}", e);
                process::exit(1);
            }
            let format_name = if use_csv { "CSV" } else { "TSV" };
            eprintln!("Summary written to: {} ({})", output_path, format_name);
            eprintln!("Total samples: {}", data_rows.len());
        }
        None => {
            println!("{}", output_content);
        }
    }
}

// ============================================================================
// Filter 命令实现
// ============================================================================

fn run_filter(args: FilterArgs) {
    let mut records = parse_path(&args.input).unwrap();

    // Regex filter from string
    if let Some(pattern) = &args.regex {
        let re = Regex::new(pattern).expect("Failed to construct regex from pattern");
        while let Some(record) = records.iter_record().unwrap() {
            let readid = record.head();
            let is_match = re.is_match(readid);
            if (is_match && !args.invert) || (!is_match && args.invert) {
                write_fastq(record);
            }
        }
        process::exit(0);
    }

    // Regex filter from file
    if let Some(refilepath) = &args.regex_file {
        let refile = fs::File::open(refilepath).expect("Regex file not found");
        let re_reader = BufReader::new(refile);
        let revec: Vec<String> = re_reader.lines().map(|l| l.unwrap()).collect();
        let re_set = RegexSet::new(&revec).unwrap();

        while let Some(record) = records.iter_record().unwrap() {
            let readid = record.head().as_bytes();
            let is_match = re_set.is_match(readid);
            if (is_match && !args.invert) || (!is_match && args.invert) {
                write_fastq(record);
            }
        }
        process::exit(0);
    }

    // Length and quality filter
    while let Some(record) = records.iter_record().unwrap() {
        let seqlen = record.seq().len() as i32;

        // Check length filters
        if let Some(min_len) = args.min_len {
            if seqlen < min_len {
                continue;
            }
        }
        if let Some(max_len) = args.max_len {
            if seqlen > max_len {
                continue;
            }
        }

        // Check quality filters
        if args.min_qual.is_some() || args.max_qual.is_some() {
            let mean_errorp = modules::qscore_probs(record.qual().as_bytes()) / record.seq().len() as f32;
            let mean_qscore = -10.0 * mean_errorp.log10();

            if let Some(min_qual) = args.min_qual {
                if mean_qscore < min_qual as f32 {
                    continue;
                }
            }
            if let Some(max_qual) = args.max_qual {
                if mean_qscore > max_qual as f32 {
                    continue;
                }
            }
        }

        write_fastq(record);
    }
}

// ============================================================================
// Trim 命令实现
// ============================================================================

fn run_trim(args: TrimArgs) {
    let mut records = parse_path(&args.input).unwrap();
    let trim_front = args.front.unwrap_or(0);
    let trim_tail = args.tail.unwrap_or(0);

    while let Some(record) = records.iter_record().unwrap() {
        let seq_len = record.len();

        // Check if trimming would result in empty or too short sequence
        if trim_front + trim_tail >= seq_len {
            continue;
        }

        let trimright = seq_len - trim_tail;
        let newseq = &record.seq()[trim_front..trimright];
        let newqual = &record.qual()[trim_front..trimright];

        // Check minimum length after trimming
        if let Some(min_len) = args.min_len {
            if newseq.len() < min_len {
                continue;
            }
        }

        println!(
            "{} {}\n{}\n{}\n{}",
            "@".to_string() + record.head(),
            record.des(),
            newseq,
            "+",
            newqual
        );
    }
}

// ============================================================================
// Sample 命令实现
// ============================================================================

fn run_sample(args: SampleArgs) {
    let mut records = parse_path(&args.input).unwrap();

    // Sample by number of reads (requires two passes or reservoir sampling)
    if let Some(target_count) = args.number {
        // Use reservoir sampling for single-pass
        let mut rng = match args.seed {
            Some(seed) => rand::rngs::StdRng::seed_from_u64(seed),
            None => rand::rngs::StdRng::from_entropy(),
        };
        use rand::SeedableRng;

        let mut reservoir: Vec<(String, String, String, String)> = Vec::with_capacity(target_count);
        let mut count: usize = 0;

        while let Some(record) = records.iter_record().unwrap() {
            count += 1;
            if reservoir.len() < target_count {
                reservoir.push((
                    format!("{} {}", "@".to_string() + record.head(), record.des()),
                    record.seq().to_string(),
                    "+".to_string(),
                    record.qual().to_string(),
                ));
            } else {
                let j = rng.r#gen_range(0..count);
                if j < target_count {
                    reservoir[j] = (
                        format!("{} {}", "@".to_string() + record.head(), record.des()),
                        record.seq().to_string(),
                        "+".to_string(),
                        record.qual().to_string(),
                    );
                }
            }
        }

        for (header, seq, plus, qual) in reservoir {
            println!("{}\n{}\n{}\n{}", header, seq, plus, qual);
        }
        process::exit(0);
    }

    // Sample by proportion
    if let Some(fraction) = args.proportion {
        if !(0.0..=1.0).contains(&fraction) {
            eprintln!("Error: Sample proportion should be between 0.0 and 1.0");
            process::exit(1);
        }

        let mut rng: rand::rngs::StdRng = match args.seed {
            Some(seed) => {
                use rand::SeedableRng;
                rand::rngs::StdRng::seed_from_u64(seed)
            }
            None => {
                use rand::SeedableRng;
                rand::rngs::StdRng::from_entropy()
            }
        };

        while let Some(record) = records.iter_record().unwrap() {
            if rng.r#gen::<f32>() < fraction {
                write_fastq(record);
            }
        }
        process::exit(0);
    }

    eprintln!("Error: Please specify -p/--proportion or -n/--number for sampling");
    process::exit(1);
}
