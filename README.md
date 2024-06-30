# SNP Melt

SNP Melt is a high-performance, multi-threaded C program designed to identify and report Single Nucleotide Polymorphisms (SNPs) from aligned sequencing data. It efficiently processes BAM files against a reference genome, providing detailed mutation information.

## Features

- Fast, multi-threaded processing of BAM files
- Identification of SNPs with customizable quality thresholds
- Efficient memory management for handling large datasets
- Output in CSV format for easy downstream analysis
- Detailed processing statistics

## Requirements

- GCC compiler with C99 support
- HTSlib library
- pthread library

## Installation

1. Ensure you have the required libraries installed on your system.
2. Clone this repository:

    git clone https://github.com/yourusername/snp-melt.git
    cd snp-melt"""

3. Compile the program:

    gcc -o snp_melt snp_melt.c -lhts -lpthread -lz -lm -O3

## Usage

    ./snp_melt -b <bam_file> -r <reference_file> -t <num_threads> [-m <max_snp_per_read>]

### Options:

- `-b <bam_file>`: Input BAM file (required)
- `-r <reference_file>`: Reference FASTA file (required)
- `-t <num_threads>`: Number of threads to use (required)
- `-m <max_snp_per_read>`: Estimated SNP count per alignment (default: 50)
- `--help`: Display help message

## Output

The program outputs a CSV file to stdout with the following columns:

1. `read_id`: Identifier of the read containing the SNP
2. `reference_position`: Position of the SNP in the reference genome
3. `mutated_base`: The base observed in the read
4. `reference_base`: The base in the reference genome
5. `base_quality`: Quality score of the mutated base

Additionally, the program prints processing statistics to stderr, including:

- Total number of reads processed
- Number and percentage of reads eliminated
- Total number of SNPs found and average SNPs per read

## Performance Considerations

- The program uses multi-threading to significantly speed up processing.
- Memory usage is optimized by processing the BAM file in chunks.
- The `-m` option allows for fine-tuning memory allocation based on expected SNP frequency.

## Limitations

- The program discards reads with indels or structural variants.
- Reads with soft clips are excluded from analysis.
- The program assumes a sorted and indexed BAM file.

## Contributing

Contributions to improve SNP Melt are welcome. Please feel free to submit pull requests or open issues to discuss potential enhancements.

## License

[Include your chosen license here]

## Citation

If you use SNP Melt in your research, please cite here.
