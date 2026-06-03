import gzip
import multiprocessing as mp
from itertools import islice
import time
import logging
import os
import glob
import sys

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Align and merge paired-end FASTQ reads."
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input directory containing sample folders"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output directory for merged FASTQ files"
)

parser.add_argument(
    "-t", "--threads",
    type=int,
    default=4,
    help="Number of worker processes (default: 4)"
)

# Define input and output files
args = parser.parse_args()
dir_path = args.input
output_dir_path = args.output
num_workers = args.threads

# Create output directory if needed
os.makedirs(output_dir_path, exist_ok=True)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

def read_fastq_pair(file1, file2, chunk_size=100):
    def read_chunks(f1, f2, size):
        while True:
            chunk = list(islice(zip(f1, f2), size * 4))
            if not chunk:
                break
            if len(chunk) < size * 4:
                logging.warning(f"Incomplete chunk of size {len(chunk)} detected.")
            yield chunk

    with (gzip.open(file1, 'rt') if file1.endswith('.gz') else open(file1, 'r')) as f1, \
         (gzip.open(file2, 'rt') if file2.endswith('.gz') else open(file2, 'r')) as f2:
        for chunk in read_chunks(f1, f2, chunk_size):
            paired_reads = []
            for i in range(0, len(chunk), 4):
                if i + 4 <= len(chunk):
                    header1, header2 = chunk[i]
                    seq1, seq2 = chunk[i+1]
                    plus1, plus2 = chunk[i+2]
                    qual1, qual2 = chunk[i+3]
                    paired_reads.append(((header1.strip(), seq1.strip(), qual1.strip()), (header2.strip(), seq2.strip(), qual2.strip())))
                else:
                    print(i)
                    logging.error("Incomplete read pair found.")
            if paired_reads:
                yield paired_reads

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def align_and_merge_reads(seq2, seq1): # swapping seq1 and seq2 to slide reverse complement of reverse strand over forward strand for OligoPool data
    len1, len2 = len(seq1), len(seq2)
    best_offset = -1
    best_matches = -1

    for offset in range(1, len1 + 1):
        matches = sum(1 for i in range(offset) if seq1[len1 - offset + i] == seq2[i] or seq1[len1 - offset + i] == 'N' or seq2[i] == 'N')
        if matches > best_matches:
            best_matches = matches
            best_offset = offset

    merged_sequence = seq1 + seq2[best_matches:]

    return merged_sequence if best_matches > 0 else None

def process_chunk(chunk):
    result = []
    errors = []
    for (header1, seq1, _), (header2, seq2, _) in chunk:
        header1_renamed = header1.replace("1:N:0", "2:N:0")
        if header1_renamed == header2:
            seq2_rc = reverse_complement(seq2)
            merged_sequence = align_and_merge_reads(seq1, seq2_rc)
            if merged_sequence is not None:
                result.append(f"{header1}\n{merged_sequence}\n+\n{'~' * len(merged_sequence)}\n")
            else:
                errors.append(f"Error: No alignment without mismatches for index {header1}")
        else:
            errors.append(f"Error: Headers do not match: {header1} vs {header2}")
    return result, errors

def process_paired_end_reads(file1, file2, output_file, num_workers=4):
    start_time = time.time()
    
    try:
        chunks = read_fastq_pair(file1, file2)

        with mp.Pool(num_workers) as pool, open(output_file, 'w') as out_f:
            for result, errors in pool.imap(process_chunk, chunks):
                for line in result:
                    out_f.write(line)
                for error in errors:
                    logging.error(error)

    except BrokenPipeError as e:
        logging.error(f"Broken pipe error: {e}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
    finally:
        end_time = time.time()
        logging.info(f"Processing completed in {end_time - start_time} seconds.")

# Loop through all folders in the directory
for folder in os.listdir(dir_path):
    folder_path = os.path.join(dir_path, folder)
    
    if os.path.isdir(folder_path):
        print(f"Processing folder: {folder_path}")

        # Grab files in the folder with expected filename
        file1, file2 = glob.glob(f"{folder_path}/*_1.fq") + glob.glob(f"{folder_path}/*_2.fq")

        # Extract the file names from the paths
        file_name1 = os.path.basename(file1)
        file_name2 = os.path.basename(file2)

        # Remove the '_R1.fastq' and '_R2.fastq' suffixes from the file names
        base_name1 = file_name1[:-13]
        base_name2 = file_name2[:-13]
        # Check if the base names are the same
        if base_name1 == base_name2:
            print("The files have the same base name.")
            output_file = os.path.join(output_dir_path, f"{base_name1}_paired_ends_merged.fastq")
            print(output_file)
            process_paired_end_reads(file1, file2, output_file, num_workers=num_workers)
        else:
            print("The files do not have the same base name.")
