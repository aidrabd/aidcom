import os
import sys
import csv
import itertools
import re
from collections import defaultdict

def parse_fasta(fasta_content):
    """
    Parse FASTA format string into a list of (header, sequence) tuples.
    This preserves all sequences, even if headers are duplicated.
    """
    sequences = []
    current_header = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        if line.startswith('>'):
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))

    return sequences

def write_fasta(sequences, output_file):
    """
    Write sequences dictionary to a FASTA file
    """
    with open(output_file, 'w') as f:
        for header, sequence in sequences.items():
            f.write(f">{header}\n")
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")

def extract_gene_info(header):
    """
    Extract gene symbol and protein name from FASTA header
    Returns a tuple (gene_symbol, protein_name)
    """
    gene_symbol = "N/A"
    protein_name = "N/A"
    
    gene_patterns = [
        r'gene[=:]([^\s|]+)',  # gene=SYMBOL or gene:SYMBOL
        r'\|([^|]+)\|.*gene',  # |SYMBOL| followed by gene
        r'GN=([^\s]+)',        # GN=SYMBOL
        r'gene_symbol[=:]([^\s|]+)',  # gene_symbol=SYMBOL
    ]
    
    for pattern in gene_patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            gene_symbol = match.group(1)
            break
    
    name_patterns = [
        r'protein[=:]([^|]+)',  # protein=NAME
        r'name[=:]([^|]+)',     # name=NAME
        r'description[=:]([^|]+)',  # description=NAME
        r'(?<=\s)(.+?)\s+OS=',  # Text before OS= (UniProt format)
        r'(?<=\s)(.+?)\s+\[',   # Text before [species]
        r'(?<=>)([^|]+)'        # Text right after > if nothing else matches
    ]
    
    for pattern in name_patterns:
        match = re.search(pattern, header, re.IGNORECASE)
        if match:
            protein_name = match.group(1).strip()
            break
    
    return gene_symbol, protein_name

def process_files(file_paths):
    """
    Process multiple FASTA files to find common and unique sequences
    """
    all_sequences = []
    file_names = []
    sequences_per_file = []
    duplicates_removed = []
    
    for file_path in file_paths:
        file_name = os.path.basename(file_path)
        file_names.append(file_name)
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            raw_sequences = parse_fasta(content)  # list of (header, seq)
            sequences_per_file.append(len(raw_sequences))  # total sequences including duplicates
            
            unique_seqs = {}
            duplicate_count = 0
            seq_content_set = set()
            
            for header, sequence in raw_sequences:
                if sequence in seq_content_set:
                    duplicate_count += 1
                else:
                    seq_content_set.add(sequence)
                    unique_seqs[header] = sequence
            
            duplicates_removed.append(duplicate_count)
            all_sequences.append(unique_seqs)
            
            print(f"Successfully parsed {file_name}: {len(raw_sequences)} sequences found, {duplicate_count} duplicates removed")
        except Exception as e:
            print(f"Error processing {file_name}: {e}")
            return None
    
    # Map sequence content to headers and files
    seq_to_header = defaultdict(list)
    seq_to_files = defaultdict(set)
    
    for i, sequences in enumerate(all_sequences):
        for header, sequence in sequences.items():
            seq_to_header[sequence].append((i, header))
            seq_to_files[sequence].add(i)
    
    # Common sequences across all files
    common_sequences = {}
    for sequence, file_indices in seq_to_files.items():
        if len(file_indices) == len(all_sequences):
            _, header = seq_to_header[sequence][0]
            common_sequences[header] = sequence
    
    # Pairwise common sequences
    pairwise_common = {}
    for i, j in itertools.combinations(range(len(all_sequences)), 2):
        pair_key = f"{file_names[i]}_AND_{file_names[j]}"
        pairwise_common[pair_key] = {}
        
        for sequence, file_indices in seq_to_files.items():
            if i in file_indices and j in file_indices:
                for file_idx, header in seq_to_header[sequence]:
                    if file_idx == i:
                        pairwise_common[pair_key][header] = sequence
                        break
    
    # Unique sequences per file
    unique_sequences = []
    for i in range(len(all_sequences)):
        unique_seqs = {}
        for header, sequence in all_sequences[i].items():
            if len(seq_to_files[sequence]) == 1 and i in seq_to_files[sequence]:
                unique_seqs[header] = sequence
        unique_sequences.append(unique_seqs)
    
    # Prepare statistics
    sequence_origins = {}
    for sequence, file_indices in seq_to_files.items():
        source_files = [file_names[i] for i in file_indices]
        _, header = seq_to_header[sequence][0]
        sequence_origins[header] = source_files
    
    common_gene_symbols = set()
    for header in common_sequences:
        gene_symbol, _ = extract_gene_info(header)
        if gene_symbol != "N/A":
            common_gene_symbols.add(gene_symbol)
    
    statistics = {
        'file_names': file_names,
        'sequences_per_file': sequences_per_file,
        'duplicates_removed': duplicates_removed,
        'unique_counts': [len(seqs) for seqs in unique_sequences],
        'common_all_count': len(common_sequences),
        'pairwise_common_counts': {k: len(v) for k, v in pairwise_common.items()},
        'sequence_origins': sequence_origins,
        'common_gene_count': len(common_gene_symbols)
    }
    
    return common_sequences, unique_sequences, pairwise_common, file_names, statistics

def save_results(common_sequences, unique_sequences, pairwise_common, file_names, statistics, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    common_file = os.path.join(output_dir, "common_sequences_all.fasta")
    write_fasta(common_sequences, common_file)
    print(f"Saved {len(common_sequences)} sequences common to all files in {common_file}")
    
    common_gene_info_file = os.path.join(output_dir, "common_sequences_gene_info.csv")
    save_gene_info_csv(common_sequences, common_gene_info_file, "Common sequences")
    print(f"Saved gene and protein information for common sequences in {common_gene_info_file}")
    
    for pair_key, sequences in pairwise_common.items():
        pair_file = os.path.join(output_dir, f"common_sequences_{pair_key}.fasta")
        write_fasta(sequences, pair_file)
        print(f"Saved {len(sequences)} sequences common to {pair_key} in {pair_file}")
        
        pair_gene_info_file = os.path.join(output_dir, f"common_sequences_{pair_key}_gene_info.csv")
        save_gene_info_csv(sequences, pair_gene_info_file, f"Common sequences between {pair_key}")
        print(f"Saved gene and protein information for {pair_key} common sequences in {pair_gene_info_file}")
    
    for i, unique_seqs in enumerate(unique_sequences):
        unique_file = os.path.join(output_dir, f"unique_sequences_{file_names[i]}")
        write_fasta(unique_seqs, unique_file)
        print(f"Saved {len(unique_seqs)} sequences unique to {file_names[i]} in {unique_file}")
        
        unique_gene_info_file = os.path.join(output_dir, f"unique_sequences_{file_names[i]}_gene_info.csv")
        save_gene_info_csv(unique_seqs, unique_gene_info_file, f"Unique sequences in {file_names[i]}")
        print(f"Saved gene and protein information for unique sequences in {file_names[i]} to {unique_gene_info_file}")
    
    stats_file = os.path.join(output_dir, "sequence_statistics.csv")
    save_statistics_csv(statistics, stats_file)
    print(f"Saved analysis statistics to {stats_file}")

def save_gene_info_csv(sequences, output_file, title):
    gene_symbol_count = 0
    protein_name_count = 0
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow([title])
        writer.writerow(['Sequence ID', 'Gene Symbol', 'Protein Name'])
        
        for header in sequences:
            gene_symbol, protein_name = extract_gene_info(header)
            if gene_symbol != "N/A":
                gene_symbol_count += 1
            if protein_name != "N/A":
                protein_name_count += 1
            writer.writerow([header, gene_symbol, protein_name])
        
        writer.writerow([])
        writer.writerow(['Summary:'])
        writer.writerow([f'Total sequences: {len(sequences)}'])
        if len(sequences) > 0:
            writer.writerow([f'Sequences with gene symbols: {gene_symbol_count} ({gene_symbol_count/len(sequences)*100:.1f}% of total)'])
            writer.writerow([f'Sequences with protein names: {protein_name_count} ({protein_name_count/len(sequences)*100:.1f}% of total)'])
        else:
            writer.writerow(['Sequences with gene symbols: 0 (0%)'])
            writer.writerow(['Sequences with protein names: 0 (0%)'])
        
        if gene_symbol_count == 0:
            writer.writerow(['Note: No gene symbols found in the sequence headers'])
        if protein_name_count == 0:
            writer.writerow(['Note: No protein names found in the sequence headers'])

def save_statistics_csv(statistics, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['SEQUENCE COUNTS'])
        writer.writerow(['File', 'Total Sequences', 'After Duplicates Removed', 'Unique Sequences'])
        
        for i, file_name in enumerate(statistics['file_names']):
            writer.writerow([
                file_name, 
                statistics['sequences_per_file'][i],
                statistics['sequences_per_file'][i] - statistics['duplicates_removed'][i],
                statistics['unique_counts'][i]
            ])
        
        writer.writerow([])
        
        writer.writerow(['DUPLICATES REMOVED'])
        for i, file_name in enumerate(statistics['file_names']):
            writer.writerow([file_name, statistics['duplicates_removed'][i]])
        
        writer.writerow([])
        
        writer.writerow(['COMMON SEQUENCES'])
        writer.writerow(['Files', 'Number of Common Sequences'])
        writer.writerow(['All files', statistics['common_all_count']])
        
        for pair, count in statistics['pairwise_common_counts'].items():
            writer.writerow([pair, count])
        
        writer.writerow([])

        writer.writerow(['COMMON GENES'])
        writer.writerow(['All files', statistics['common_gene_count']])
        
        writer.writerow([])
        
        writer.writerow(['SEQUENCE SOURCE SUMMARY'])
        writer.writerow(['Sequence ID', 'Present in Files'])
        
        for header, files in list(statistics['sequence_origins'].items())[:100]:
            writer.writerow([header, ', '.join(files)])
        
        if len(statistics['sequence_origins']) > 100:
            writer.writerow(['...', f'(Showing 100 of {len(statistics["sequence_origins"])} sequences)'])

def process_manual_input():
    file_contents = []
    file_names = []

    try:
        num_files = int(input("Enter the number of FASTA files to process: "))
        if num_files <= 0:
            print("Number of files must be positive. Exiting.")
            return None
    except ValueError:
        print("Invalid number. Exiting.")
        return None

    for i in range(num_files):
        file_name = input(f"Enter name for file {i + 1}: ")
        file_names.append(file_name)

        if os.path.exists(file_name):
            print(f"Reading FASTA content from existing file: {file_name}")
            try:
                with open(file_name, 'r') as f:
                    content = f.read()
                file_contents.append(content)
            except Exception as e:
                print(f"Error reading file {file_name}: {e}")
                return None
        else:
            print(f"Enter FASTA content for {file_name} (end with a line containing only 'END'):")
            lines = []
            while True:
                line = input()
                if line == "END":
                    break
                lines.append(line)
            file_contents.append('\n'.join(lines))

    if not file_names:
        print("No files provided. Exiting.")
        return None

    all_sequences = []
    sequences_per_file = []
    duplicates_removed = []

    for i, content in enumerate(file_contents):
        raw_sequences = parse_fasta(content)
        sequences_per_file.append(len(raw_sequences))

        unique_seqs = {}
        duplicate_count = 0
        seq_content_set = set()

        for header, sequence in raw_sequences:
            if sequence in seq_content_set:
                duplicate_count += 1
            else:
                seq_content_set.add(sequence)
                unique_seqs[header] = sequence

        duplicates_removed.append(duplicate_count)
        all_sequences.append(unique_seqs)

        print(f"Processed {file_names[i]}: {len(raw_sequences)} sequences found, {duplicate_count} duplicates removed")

    seq_to_header = defaultdict(list)
    seq_to_files = defaultdict(set)

    for i, sequences in enumerate(all_sequences):
        for header, sequence in sequences.items():
            seq_to_header[sequence].append((i, header))
            seq_to_files[sequence].add(i)

    common_sequences = {}
    for sequence, file_indices in seq_to_files.items():
        if len(file_indices) == len(all_sequences):
            _, header = seq_to_header[sequence][0]
            common_sequences[header] = sequence

    pairwise_common = {}
    for i, j in itertools.combinations(range(len(all_sequences)), 2):
        pair_key = f"{file_names[i]}_AND_{file_names[j]}"
        pairwise_common[pair_key] = {}

        for sequence, file_indices in seq_to_files.items():
            if i in file_indices and j in file_indices:
                for file_idx, header in seq_to_header[sequence]:
                    if file_idx == i:
                        pairwise_common[pair_key][header] = sequence
                        break

    unique_sequences = []
    for i in range(len(all_sequences)):
        unique_seqs = {}
        for header, sequence in all_sequences[i].items():
            if len(seq_to_files[sequence]) == 1 and i in seq_to_files[sequence]:
                unique_seqs[header] = sequence
        unique_sequences.append(unique_seqs)

    sequence_origins = {}
    for sequence, file_indices in seq_to_files.items():
        source_files = [file_names[i] for i in file_indices]
        _, header = seq_to_header[sequence][0]
        sequence_origins[header] = source_files

    common_gene_symbols = set()
    for header in common_sequences:
        gene_symbol, _ = extract_gene_info(header)
        if gene_symbol != "N/A":
            common_gene_symbols.add(gene_symbol)
    
    statistics = {
        'file_names': file_names,
        'sequences_per_file': sequences_per_file,
        'duplicates_removed': duplicates_removed,
        'unique_counts': [len(seqs) for seqs in unique_sequences],
        'common_all_count': len(common_sequences),
        'pairwise_common_counts': {k: len(v) for k, v in pairwise_common.items()},
        'sequence_origins': sequence_origins,
        'common_gene_count': len(common_gene_symbols)
    }
    
    return common_sequences, unique_sequences, pairwise_common, file_names, statistics

def main():
    print("Protein Sequence Analyzer")
    print("========================")
    print("This program analyzes protein sequences from FASTA files to identify common and unique sequences.")
    print("\nChoose an input method:")
    print("1. File location (provide paths to FASTA files)")
    print("2. Manual input")
    
    choice = input("\nEnter your choice (1 or 2): ")
    
    output_dir = input("Enter output directory path (default: './sequence_results'): ") or "./sequence_results"
    
    if choice == "1":
        print("\nEnter file paths (one per line, enter 'DONE' when finished):")
        print("Leave empty to use current directory as default path")
        file_paths = []
        
        while True:
            path = input()
            if path.upper() == "DONE":
                break
            elif path == "":
                current_dir = os.getcwd()
                print(f"Using current directory: {current_dir}")
                for file in os.listdir(current_dir):
                    if file.endswith('.fasta') or file.endswith('.fa'):
                        file_path = os.path.join(current_dir, file)
                        file_paths.append(file_path)
                        print(f"Found FASTA file: {file}")
                break
            elif os.path.exists(path):
                file_paths.append(path)
            else:
                print(f"Warning: File {path} not found!")
        
        if not file_paths:
            print("No valid files provided. Exiting.")
            return
        
        results = process_files(file_paths)
        if results:
            common_sequences, unique_sequences, pairwise_common, file_names, statistics = results
            print("\nSummary of sequence counts:")
            print(f"Total sequences across all files (before duplicates removed): {sum(statistics['sequences_per_file'])}")
            print(f"Sequences common to all files: {statistics['common_all_count']}")
            print(f"Unique sequences across all files: {sum(statistics['unique_counts'])}")
            save_results(common_sequences, unique_sequences, pairwise_common, file_names, statistics, output_dir)
    
    elif choice == "2":
        results = process_manual_input()
        if results:
            common_sequences, unique_sequences, pairwise_common, file_names, statistics = results
            print("\nSummary of sequence counts:")
            print(f"Total sequences across all files (before duplicates removed): {sum(statistics['sequences_per_file'])}")
            print(f"Sequences common to all files: {statistics['common_all_count']}")
            print(f"Unique sequences across all files: {sum(statistics['unique_counts'])}")
            save_results(common_sequences, unique_sequences, pairwise_common, file_names, statistics, output_dir)
    
    else:
        print("Invalid choice. Exiting.")

if __name__ == "__main__":
    main()