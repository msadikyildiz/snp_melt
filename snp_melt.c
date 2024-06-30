#define _GNU_SOURCE  // Ensure GNU-specific extensions are included

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <ctype.h>
#include <unistd.h>
#include <locale.h>

#define DEFAULT_MAX_SNP_PER_READ 50

// Structure for storing mutation information
typedef struct {
    char read_id[64];
    int reference_position;
    char mutated_base;
    char reference_base;
    int base_quality;
} Mutation;

// Structure for passing data to each thread
typedef struct {
    const char *bam_filename;
    const char *ref_filename;
    int thread_id;
    int start_alignment;
    int end_alignment;
    Mutation *mutations;
    int mutation_count;
    long reads_processed;
    long reads_eliminated;
    int max_mutations_per_thread;
} ThreadData;

// Function to get the total number of alignments in the BAM file
int get_num_alignments(const char *bam_filename) {
    samFile *bam_file = sam_open(bam_filename, "r");
    if (bam_file == NULL) {
        fprintf(stderr, "Failed to open BAM file: %s\n", bam_filename);
        return -1;
    }
    
    bam_hdr_t *header = sam_hdr_read(bam_file);
    if (header == NULL) {
        fprintf(stderr, "Failed to read BAM header for file: %s\n", bam_filename);
        sam_close(bam_file);
        return -1;
    }

    int num_alignments = 0;
    bam1_t *alignment = bam_init1();
    while (sam_read1(bam_file, header, alignment) >= 0) {
        num_alignments++;
    }
    
    bam_destroy1(alignment);
    bam_hdr_destroy(header);
    sam_close(bam_file);

    return num_alignments;
}

// Function to decode BAM nucleotide base
char decode_base(uint8_t base) {
    switch (base) {
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        default: return 'N';
    }
}

// Function to convert a sequence to uppercase
void to_uppercase(char *seq) {
    for (int i = 0; seq[i]; i++) {
        seq[i] = toupper(seq[i]);
    }
}

// Function to process BAM file in a thread
void *process_bam(void *arg) {
    ThreadData *data = (ThreadData *)arg;

    samFile *bam_file = sam_open(data->bam_filename, "r"); // Open BAM file for reading
    if (bam_file == NULL) {
        fprintf(stderr, "Thread %d: Failed to open BAM file: %s\n", data->thread_id, data->bam_filename);
        pthread_exit(NULL);
    }

    bam_hdr_t *header = sam_hdr_read(bam_file); // Read header
    bam1_t *alignment = bam_init1(); // Initialize an alignment structure

    hts_idx_t *idx = sam_index_load(bam_file, data->bam_filename);
    if (!idx) {
        fprintf(stderr, "Thread %d: Failed to load BAM index for file: %s\n", data->thread_id, data->bam_filename);
        pthread_exit(NULL);
    }

    // Create iterator to cover all alignments
    hts_itr_t *iter = sam_itr_queryi(idx, HTS_IDX_START, 0, HTS_IDX_NOCOOR);
    if (!iter) {
        fprintf(stderr, "Thread %d: Failed to create iterator for BAM file: %s\n", data->thread_id, data->bam_filename);
        pthread_exit(NULL);
    }

    // Load the reference sequence
    faidx_t *fai = fai_load(data->ref_filename);
    if (fai == NULL) {
        fprintf(stderr, "Thread %d: Failed to load reference sequence: %s\n", data->thread_id, data->ref_filename);
        pthread_exit(NULL);
    }

    // Skip to the start alignment for this thread
    for (int i = 0; i < data->start_alignment; i++) {
        if (sam_itr_next(bam_file, iter, alignment) < 0) {
            pthread_exit(NULL);
        }
    }

    // Process the alignments assigned to this thread
    for (int i = data->start_alignment; i < data->end_alignment; i++) {
        if (sam_itr_next(bam_file, iter, alignment) < 0) {
            break;
        }

        data->reads_processed++;
        uint32_t *cigar = bam_get_cigar(alignment);
        int num_cigar = alignment->core.n_cigar;

        // Check for indels or soft-clipped ends
        int has_indel_or_soft_clip = 0;
        for (int j = 0; j < num_cigar; ++j) {
            int cigar_op = bam_cigar_op(cigar[j]);
            if (cigar_op == BAM_CINS || cigar_op == BAM_CDEL || cigar_op == BAM_CSOFT_CLIP) {
                has_indel_or_soft_clip = 1;
                break;
            }
        }

        if (has_indel_or_soft_clip || (alignment->core.flag & BAM_FUNMAP)) {
            // Skip alignments with indels, soft clips, or unmapped reads
            data->reads_eliminated++;
            continue;
        }

        // Get reference sequence name and position
        const char *ref_name = header->target_name[alignment->core.tid];
        int ref_start = alignment->core.pos;
        int ref_end = ref_start + alignment->core.l_qseq;

        // Fetch reference sequence
        int seq_len;
        char *ref_seq = faidx_fetch_seq(fai, ref_name, ref_start, ref_end - 1, &seq_len);

        // Convert reference sequence to uppercase
        to_uppercase(ref_seq);

        // Check for SNPs
        uint8_t *seq = bam_get_seq(alignment);
        uint8_t *qual = bam_get_qual(alignment);
        int read_pos = 0;

        for (int j = 0; j < num_cigar; ++j) {
            int cigar_op = bam_cigar_op(cigar[j]);
            int cigar_len = bam_cigar_oplen(cigar[j]);

            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                for (int k = 0; k < cigar_len; ++k) {
                    int base = bam_seqi(seq, read_pos);
                    int base_quality = qual[read_pos];
                    char decoded_base = decode_base(base);
                    char ref_base = ref_seq[ref_start + read_pos - ref_start];

                    // Check for mismatch
                    if (decoded_base != ref_base && base_quality >= 20) { // Assuming quality >= 20
                        // Report SNP mutation
                        Mutation mutation;
                        snprintf(mutation.read_id, sizeof(mutation.read_id), "%s", bam_get_qname(alignment));
                        mutation.reference_position = alignment->core.pos + read_pos + 1;
                        mutation.mutated_base = decoded_base;
                        mutation.reference_base = ref_base;
                        mutation.base_quality = base_quality;

                        // Add mutation to the local list
                        if (data->mutation_count < data->max_mutations_per_thread) {
                            data->mutations[data->mutation_count++] = mutation;
                        }
                    }

                    read_pos++;
                }
            } else if (cigar_op == BAM_CINS || cigar_op == BAM_CDEL || cigar_op == BAM_CSOFT_CLIP) {
                read_pos += cigar_len;
            }
        }

        free(ref_seq);
    }

    fai_destroy(fai);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_destroy1(alignment); // Clean up alignment structure
    bam_hdr_destroy(header); // Clean up header
    sam_close(bam_file); // Close BAM file

    pthread_exit(NULL);
}

// Function to print the help message
void print_help() {
    printf("Usage: snp_melt -b <bam_file> -r <reference_file> -t <num_threads> [-m <max_snp_per_read>]\n");
    printf("Options:\n");
    printf("  -b <bam_file>         Input BAM file\n");
    printf("  -r <reference_file>   Reference FASTA file\n");
    printf("  -t <num_threads>      Number of threads to use\n");
    printf("  -m <max_snp_per_read>  Estimated SNP count per alignment (default: 100)\n");
    printf("  --help                Display this help message\n");
}

int main(int argc, char *argv[]) {
    const char *bam_filename = NULL;
    const char *ref_filename = NULL;
    int num_threads = 0;
    int estimated_snp_count_per_alignment = DEFAULT_MAX_SNP_PER_READ;

    // Set the locale to the user's default locale
    setlocale(LC_NUMERIC, "");

    // Parse command-line arguments
    int opt;
    while ((opt = getopt(argc, argv, "b:r:t:m:")) != -1) {
        switch (opt) {
            case 'b':
                bam_filename = optarg;
                break;
            case 'r':
                ref_filename = optarg;
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'm':
                estimated_snp_count_per_alignment = atoi(optarg);
                break;
            default:
                print_help();
                return 1;
        }
    }

    // Check if all required arguments are provided
    if (bam_filename == NULL || ref_filename == NULL || num_threads <= 0) {
        print_help();
        return 1;
    }

    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    // Get the total number of alignments in the BAM file
    long long int num_alignments = get_num_alignments(bam_filename);
    if (num_alignments < 0) {
        fprintf(stderr, "Failed to get the number of alignments in the BAM file: %s\n", bam_filename);
        return 1;
    }

    // Calculate mutations per thread
    long long int max_mutations_per_thread = (num_alignments * estimated_snp_count_per_alignment) / num_threads;

    // Print memory usage estimate
    fprintf(stderr, "Number of reads: %'lld\n", num_alignments);
    fprintf(stderr, "Memory allocated: %lld GB\n",
            (sizeof(Mutation) * max_mutations_per_thread * num_threads / (1024 * 1024 * 1024)));

    // Calculate alignments per thread
    int alignments_per_thread = num_alignments / num_threads;
    int remaining_alignments = num_alignments % num_threads;

    // Allocate memory for mutations per thread
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].mutations = malloc(sizeof(Mutation) * max_mutations_per_thread);
        if (thread_data[i].mutations == NULL) {
            fprintf(stderr, "Failed to allocate memory for mutations\n");
            return 1;
        }
        thread_data[i].mutation_count = 0;
        thread_data[i].reads_processed = 0;
        thread_data[i].reads_eliminated = 0;
        thread_data[i].max_mutations_per_thread = max_mutations_per_thread;
    }

    // Set up thread data
    int start_alignment = 0;
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].bam_filename = bam_filename;
        thread_data[i].ref_filename = ref_filename;
        thread_data[i].thread_id = i;
        thread_data[i].start_alignment = start_alignment;
        thread_data[i].end_alignment = start_alignment + alignments_per_thread + (i < remaining_alignments ? 1 : 0);
        start_alignment = thread_data[i].end_alignment;

        if (pthread_create(&threads[i], NULL, process_bam, (void *)&thread_data[i]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i);
            return 1;
        }
    }

    // Wait for threads to complete
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    // Write mutations to standard output
    printf("read_id,reference_position,mutated_base,reference_base,base_quality\n"); // CSV header

    long total_reads_processed = 0;
    long total_reads_eliminated = 0;
    long total_mutations = 0;

    for (int i = 0; i < num_threads; ++i) {
        total_reads_processed += thread_data[i].reads_processed;
        total_reads_eliminated += thread_data[i].reads_eliminated;
        total_mutations += thread_data[i].mutation_count;

        for (int j = 0; j < thread_data[i].mutation_count; ++j) {
            Mutation *mutation = &thread_data[i].mutations[j];
            printf("%s,%d,%c,%c,%d\n", mutation->read_id, mutation->reference_position, mutation->mutated_base, mutation->reference_base, mutation->base_quality);
        }
        free(thread_data[i].mutations);
    }

    // Print report to stderr
    double percent_eliminated = ((double)total_reads_eliminated / total_reads_processed) * 100;
    fprintf(stderr, "Reads processed:\t%'ld\n", total_reads_processed);
    fprintf(stderr, "Reads eliminated:\t%'ld\t(%.2f%%)\n", total_reads_eliminated, percent_eliminated);
    double snp_count_per_read = (double)total_mutations / (total_reads_processed-total_reads_eliminated);
    fprintf(stderr, "snps found:\t%'ld\t(%.1f snps per read)\n", total_mutations, snp_count_per_read);

    return 0;
}
