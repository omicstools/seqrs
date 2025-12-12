# seqrs

Fast statistics and filtering tool for FASTQ/FASTA files, written in Rust.

## Features

- **Statistics**: Basic stats, QC reports with plots, per-read metrics (length, GC, quality)
- **Filtering**: By length, quality score, or read ID (regex)
- **Trimming**: Remove bases from 5' and/or 3' ends
- **Sampling**: Random sub-sampling by proportion or count
- **Summarize**: Merge multiple QC report files into one table

## Installation

### From source

```bash
git clone https://github.com/your-repo/seqrs.git
cd seqrs
cargo build --release
cargo build --release --target x86_64-unknown-linux-musl 
```

The binary will be at `target/release/seqrs`.

### Dependencies

- Rust 1.70+ (2024 edition)

## Usage

```
seqrs <COMMAND> [OPTIONS] <INPUT>

Commands:
  stat    Statistics and QC reports
  filter  Filter reads by length, quality, or ID
  trim    Trim reads from front or tail
  sample  Sub-sample reads
```

### stat - Statistics and QC Reports

```bash
# Basic statistics (table output)
seqrs stat input.fastq
seqrs stat -t input.fastq           # explicit table format
seqrs stat -ts input.fastq          # table without header

# Per-read metrics
seqrs stat -l input.fastq           # output lengths (one per line)
seqrs stat -g input.fastq           # output GC content
seqrs stat -q input.fastq           # output quality scores

# Specific metrics
seqrs stat -x 0.5 input.fastq       # calculate N50
seqrs stat -x 0.9 input.fastq       # calculate N90
seqrs stat -y 20 input.fastq        # calculate Q20 yield

# QC report with plots
seqrs stat --qc sample1 input.fastq
seqrs stat --qc sample1 -o ./reports input.fastq

# Summarize multiple QC reports
seqrs stat --sum-dir ./qc_reports                 # output to stdout (TSV)
seqrs stat --sum-dir ./qc_reports -o summary.tsv  # output to TSV file
seqrs stat --sum-dir ./qc_reports -o summary.csv  # output to CSV file (auto-detect)
```

**Options:**
| Option | Description |
|--------|-------------|
| `-t, --table` | Output tab-separated statistics table |
| `-s, --skip-header` | Skip header line in table output |
| `--qc <NAME>` | Generate QC report (creates .qc_stat.txt and .length_dist.svg) |
| `-o, --output <PATH>` | Output directory (for --qc) or output file (for --sum-dir) |
| `-l, --lengths` | Output read lengths (one per line) |
| `-g, --gc` | Output GC content per read |
| `-q, --quality` | Output mean Phred quality score per read |
| `-x, --nx <VALUE>` | Calculate NX value (e.g., 0.5 for N50) |
| `-y, --qyield <Q>` | Calculate percent bases >= Q score [8-60] |
| `--sum-dir <DIR>` | Summarize .qc_stat.txt files from directory |

### filter - Filter Reads

```bash
# Filter by length
seqrs filter -l 1000 input.fastq              # keep reads >= 1000bp
seqrs filter -L 500 input.fastq               # keep reads <= 500bp
seqrs filter -l 1000 -L 5000 input.fastq      # keep reads 1000-5000bp

# Filter by quality
seqrs filter -q 20 input.fastq                # keep reads with mean Q >= 20
seqrs filter -Q 30 input.fastq                # keep reads with mean Q <= 30

# Filter by read ID (regex)
seqrs filter -r 'chr1' input.fastq            # keep reads matching pattern
seqrs filter -R patterns.txt input.fastq      # patterns from file
seqrs filter -r 'chr1' -v input.fastq         # exclude matching reads
```

**Options:**
| Option | Description |
|--------|-------------|
| `-l, --min-len <N>` | Keep reads >= N bp |
| `-L, --max-len <N>` | Keep reads <= N bp |
| `-q, --min-qual <Q>` | Keep reads with mean quality >= Q |
| `-Q, --max-qual <Q>` | Keep reads with mean quality <= Q |
| `-r, --regex <PAT>` | Filter by ID matching regex pattern |
| `-R, --regex-file <FILE>` | Filter by patterns from file |
| `-v, --invert` | Invert match (exclude matching reads) |

### trim - Trim Reads

```bash
seqrs trim -f 50 input.fastq                  # trim 50bp from 5' end
seqrs trim -t 30 input.fastq                  # trim 30bp from 3' end
seqrs trim -f 50 -t 30 input.fastq            # trim both ends
seqrs trim -f 50 -l 100 input.fastq           # trim and keep reads >= 100bp
```

**Options:**
| Option | Description |
|--------|-------------|
| `-f, --front <N>` | Trim N bases from 5' end |
| `-t, --tail <N>` | Trim N bases from 3' end |
| `-l, --min-len <N>` | Minimum length after trimming |

### sample - Sub-sample Reads

```bash
seqrs sample -p 0.1 input.fastq               # 10% random sample
seqrs sample -n 1000 input.fastq              # sample 1000 reads
seqrs sample -p 0.1 -s 42 input.fastq         # reproducible sampling with seed
```

**Options:**
| Option | Description |
|--------|-------------|
| `-p, --proportion <FRAC>` | Sample by proportion [0.0-1.0] |
| `-n, --number <N>` | Sample by number of reads |
| `-s, --seed <SEED>` | Random seed for reproducibility |

## Examples

### Workflow: QC multiple samples and summarize

```bash
# Generate QC reports for multiple samples
for f in *.fastq.gz; do
    sample=$(basename "$f" .fastq.gz)
    seqrs stat --qc "$sample" -o ./qc_reports "$f"
done

# Summarize all QC reports into one table
seqrs stat --sum-dir ./qc_reports -o qc_summary.csv
```

### Workflow: Filter and trim

```bash
# Filter long, high-quality reads and trim adapters
seqrs filter -l 1000 -q 20 input.fastq | seqrs trim -f 50 -t 50 > clean.fastq
```

### Workflow: Random subsample for testing

```bash
# Take 10% subsample with fixed seed for reproducibility
seqrs sample -p 0.1 -s 42 large_file.fastq > subset.fastq
```

## Output Formats

### Statistics Table

```
file    reads   bases   n_bases min_len max_len mean_len    Q1  Q2  Q3  N50 Q20_percent Q30_percent GC_percent
input.fastq 10  18931   0   165 8490    1893.10 249 453 2440    5263    51.45   16.33   38.19
```

### QC Report Files

- `<name>.qc_stat.txt` - Tab-separated statistics
- `<name>.length_dist.svg` - Length distribution plot with cumulative curve

## Development

### Build

```bash
# Debug build
cargo build

# Release build (optimized)
cargo build --release

# Run tests
cargo test
```

### Project Structure

```
seqrs/
├── src/
│   ├── main.rs       # CLI and command implementations
│   ├── modules.rs    # Statistics calculation utilities
│   └── plotting.rs   # SVG plot generation
├── tests/
│   ├── cli.rs        # Integration tests
│   └── test.fastq    # Test data
├── Cargo.toml
└── README.md
```

### Dependencies

- `clap` - Command-line argument parsing
- `kseq` - FASTQ/FASTA parsing
- `bio` - Bioinformatics utilities
- `rayon` - Parallel processing
- `plotters` - SVG plotting
- `regex` - Regular expressions
- `indicatif` - Progress bars
- `rand` - Random number generation

## License

MIT

## Author

Rui Wang
