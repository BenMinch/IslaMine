# IslaMine

IslaMine is a Python-based tool for **detecting genomic islands (GIs)** in assembled genomes using **tetranucleotide frequency (TNF) deviation**, cumulative sum analysis, and adaptive variance-based peak detection. The pipeline is designed to work on *viral, bacterial, or microbial genomes* and emphasizes robustness to noise and assembly artifacts.

In addition to identifying candidate genomic islands, IslaMine:

* Classifies islands by confidence level
* Detects **direct repeats** flanking islands (hallmark of horizontal transfer)
* Produces publication-ready plots
* Outputs island coordinates and FASTA sequences
* Generates per-genome summary statistics

---

## Requirements

Python ≥ 3.8

Required Python packages:

```bash
pip install numpy pandas matplotlib biopython scipy
```

---

## Input

* A directory containing genome FASTA files:

  * `.fasta`, `.fa`, or `.fna`
  * One genome per file

Example:

```
genomes/
├── genome1.fasta
├── genome2.fna
└── genome3.fa
```

---

## Usage

```bash
python IslaMine.py <input_fasta_folder> <output_directory>
```

### Example

```bash
python IslaMine.py genomes/ IslaMine_output/
```

Genomes shorter than **50 kb** are automatically skipped to reduce false positives from truncated assemblies.

---

## Output Structure

For each genome, IslaMine creates a dedicated output directory:

```
IslaMine_output/
└── GenomeA/
    ├── genome_analysis_5000_islands.png
    ├── genome_analysis_2000_islands.png
    ├── genomic_islands_coordinates_5000.csv
    ├── genomic_islands_coordinates_2000.csv
    ├── GenomeA_identified_islands_5kb.fasta
    ├── GenomeA_identified_islands_2kb.fasta
```

Additionally, a **global summary file** is produced:

```
genome_summary.csv
```

---

## Output Files Explained

### Island Coordinates CSV

`genomic_islands_coordinates_<chunk>.csv`

| Column              | Description                           |
| ------------------- | ------------------------------------- |
| Start               | Island start position (bp)            |
| End                 | Island end position (bp)              |
| Peak Type           | High-Confidence / Primary / Secondary |
| Direct Repeat Found | Yes / No                              |
| Repeat Sequence     | Repeat motif (if found)               |
| Repeat Length       | Length of repeat                      |
| Upstream Position   | Genomic position of upstream repeat   |
| Downstream Position | Genomic position of downstream repeat |

---

### FASTA of Identified Islands

Each island is exported as an individual FASTA record for downstream annotation or comparative analysis.

---

### Diagnostic Plots

Plots include:

1. Pearson correlation per genome chunk
2. Cumulative sum (CUSUM)
3. First derivative of CUSUM
4. Sliding-window variance with islands highlighted

Color coding:

* **Blue**: High-confidence islands (direct repeats)
* **Red**: Primary islands
* **Orange**: Secondary islands

---

### Genome Summary CSV

`genome_summary.csv`

Includes per-genome statistics such as:

* Total genome size
* Number and total size of:

  * High-confidence islands
  * Primary islands
  * Secondary islands

---

## Interpretation Notes

* **High-confidence islands** (direct repeats detected) are strong candidates for horizontal gene transfer
* Islands near contig edges should be interpreted cautiously

---

## Recommended Downstream Analyses

IslaMine is designed to integrate cleanly with:

* Phylogenetic clustering (species-level)
* Genome length dispersion filtering
* Core gene completeness analysis
* Island presence/absence matrices
* Synteny-aware island comparison

This enables confident discrimination between:

> True biological absence vs. assembly-driven false negatives

---

## Caveats

* Requires reasonably complete genomes (>50 kb)
* Not optimized for highly fragmented metagenome bins
* Island boundaries are approximate and best used comparatively

---

## Citation

If you use IslaMine in published work, please cite appropriately or link to this repository.

---

## Contributing

Pull requests, issues, and feature suggestions are welcome.

---

## Contact

For questions, feature requests, or collaboration ideas, please open an issue or contact the repository maintainer.
