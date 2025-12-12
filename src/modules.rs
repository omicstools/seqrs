// simple helper functions for calcuting mean, quartiles etc
use rayon::prelude::*;
use kseq::record::Fastx;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;
pub fn mean(numbers: &[i64]) -> f64 {
    numbers.par_iter().sum::<i64>() as f64 / numbers.len() as f64
}

pub fn quartiles(numbers: &mut [i64], q: i8) -> i64 {
    numbers.sort_unstable();

    match q {
        1 => {
            let index = numbers.len() / 4;
            return numbers[index];
        }
        2 => {
            let index = numbers.len() / 2;
            return numbers[index];
        }
        3 => {
            // avoid having to use f64
            let index1 = numbers.len() / 4;
            let index2 = numbers.len() / 2;
            return numbers[index1 + index2];
        }
        _ => 42, //:)
    }
    // first quartile
}

pub fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.par_iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers
        .iter()
        .scan(0, |sum, i| {
            *sum += i;
            Some(*sum)
        })
        .collect::<Vec<_>>();
    let n50_index = cumsum.par_iter().position_first(|&x| x > halfsum as i64).unwrap();

    numbers[n50_index]
}

// get number of bases with q >= value
pub fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
    let mut n = 0;
    for item in q {
        if *item >= qx {
            n += 1
        }
    }
    n
}

// get number of N bases
pub fn get_n_bases(seq: &[u8]) -> i32 {
    let mut n = 0;
    for s in seq {
        if s == &78u8 || s == &110u8 {
            n += 1;
        }
    }
    n
}

// get number of GC bases (case-insensitive)
pub fn get_gc_bases(seq: &[u8]) -> i64 {
    let mut n = 0;
    for s in seq {
        // 'G' = 71, 'g' = 103, 'C' = 67, 'c' = 99
        if *s == 71u8 || *s == 103u8 || *s == 67u8 || *s == 99u8 {
            n += 1;
        }
    }
    n
}

// to get mean of q scores from a record - first convert to prob, calc mean, then back to phred
// this fn reads phred and converts to probs and returns their sum
//
// see how seqkit is doing it
// https://github.com/shenwei356/bio/blob/1886d4a9eab7315f6f445595acbdc7bc3edf0e08/seq/seq.go#L727

pub fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum = 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    qprob_sum
}

pub fn write_fastq(rec: Fastx<'_>) {
    println!(
        "{} {}\n{}\n{}\n{}", 
        "@".to_string() + rec.head(), rec.des(), 
        rec.seq(), 
        "+", 
        rec.qual()
    );
}

/// QC statistics structure
pub struct QcStats {
    pub sample_name: String,
    pub reads: i64,
    pub bases: i64,
    pub n_bases: i32,
    pub min_len: i64,
    pub max_len: i64,
    pub mean_len: f64,
    pub q1: i64,
    pub q2: i64,
    pub q3: i64,
    pub n50: i64,
    pub q20_percent: f64,
    pub q30_percent: f64,
    pub gc_percent: f64,
}

/// Write QC statistics to file
pub fn write_qc_stats(stats: &QcStats, output_path: &Path) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_path)?;

    writeln!(file, "sample\treads\tbases\tn_bases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent\tGC_percent")?;
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
        stats.sample_name,
        stats.reads,
        stats.bases,
        stats.n_bases,
        stats.min_len,
        stats.max_len,
        stats.mean_len,
        stats.q1,
        stats.q2,
        stats.q3,
        stats.n50,
        stats.q20_percent,
        stats.q30_percent,
        stats.gc_percent
    )?;

    Ok(())
}

