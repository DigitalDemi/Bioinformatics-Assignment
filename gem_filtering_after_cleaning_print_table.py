import logging
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

def find_orfs_in_protein(aa_sequence, dna_sequence=None, frame=0, is_reverse=False, min_length=30):
    """
    Find all ORFs in an amino acid sequence.
    
    Args:
        aa_sequence (str): Amino acid sequence
        dna_sequence (str, optional): Original DNA sequence for position mapping
        frame (int, optional): Reading frame (0, 1, or 2)
        is_reverse (bool, optional): Whether this is from the reverse complement
        min_length (int, optional): Minimum ORF length in amino acids
        
    Returns:
        list: List of dictionaries with ORF information
    """
    orfs = []
    i = 0
    
    while i < len(aa_sequence):
        # Find next start codon (M)
        start_aa_pos = i
        while start_aa_pos < len(aa_sequence) and aa_sequence[start_aa_pos] != 'M':
            start_aa_pos += 1
        
        if start_aa_pos >= len(aa_sequence):
            break  # No more start codons
        
        # Calculate DNA positions
        dna_start_pos = frame + (start_aa_pos * 3) if dna_sequence else -1
        
        # Continue until stop codon or end of sequence
        end_aa_pos = start_aa_pos + 1
        while end_aa_pos < len(aa_sequence) and aa_sequence[end_aa_pos] != '*':
            end_aa_pos += 1
        
        # Calculate DNA end position
        dna_end_pos = frame + (end_aa_pos * 3) + 2 if dna_sequence else -1
        
        # Adjust positions for reverse complement
        if is_reverse and dna_sequence:
            orig_dna_start = len(dna_sequence) - dna_end_pos - 1
            orig_dna_end = len(dna_sequence) - dna_start_pos - 1
            dna_start_pos, dna_end_pos = orig_dna_start, orig_dna_end
        
        # Get the ORF sequence
        orf_sequence = aa_sequence[start_aa_pos:end_aa_pos]
        if end_aa_pos < len(aa_sequence):
            orf_sequence += '*'  # Add stop codon if found
        
        # Only include ORFs that meet minimum length
        if len(orf_sequence) >= min_length:
            # Get DNA sequence of the ORF if available
            dna_sequence_of_orf = None
            if dna_sequence and dna_start_pos >= 0 and dna_end_pos >= 0:
                if not is_reverse:
                    dna_sequence_of_orf = dna_sequence[dna_start_pos:dna_end_pos+1]
                else:
                    # For reverse complement, we need to take the original segment and then reverse complement it
                    orig_segment = dna_sequence[dna_end_pos:dna_start_pos+1]
                    dna_sequence_of_orf = reverse_complement(orig_segment)
            
            orfs.append({
                'type': 'complete' if end_aa_pos < len(aa_sequence) else 'incomplete',
                'aa_start': start_aa_pos,
                'aa_end': end_aa_pos if end_aa_pos < len(aa_sequence) else len(aa_sequence) - 1,
                'dna_start': dna_start_pos,
                'dna_end': dna_end_pos,
                'length_aa': len(orf_sequence) - (1 if end_aa_pos < len(aa_sequence) else 0),
                'length_nt': (dna_end_pos - dna_start_pos + 1) if dna_sequence else -1,
                'sequence': orf_sequence,
                'dna_sequence': dna_sequence_of_orf
            })
        
        # Move to position after this ORF
        i = end_aa_pos + 1 if end_aa_pos < len(aa_sequence) else len(aa_sequence)
    
    return orfs

def reverse_complement(sequence):
    """
    Generate the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): DNA sequence (5'→3')
        
    Returns:
        str: Reverse complement of the sequence (5'→3')
    """
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    # Clean the sequence
    sequence = sequence.replace(" ", "").upper()
    
    # Generate complement and reverse it
    complement = ''.join([complement_dict.get(base, 'N') for base in sequence])
    reverse_comp = complement[::-1]  # Reverse to get 5'→3' orientation
    
    return reverse_comp

def translate_dna(sequence, frame=0):
    """
    Translate a DNA sequence into an amino acid sequence starting from a specified frame.
    
    Args:
        sequence (str): DNA sequence
        frame (int): Reading frame (0, 1, or 2)
        
    Returns:
        str: Amino acid sequence
    """
    # Dictionary of codon to amino acid mappings
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    
    protein = ""
    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:  # Make sure we have a complete codon
            amino_acid = codon_table.get(codon, 'X')  # 'X' for unknown codons
            protein += amino_acid
    
    return protein

def find_orfs_in_protein(aa_sequence, dna_sequence=None, frame=0, is_reverse=False, min_length=30):
    """
    Find all ORFs in an amino acid sequence.
    
    Args:
        aa_sequence (str): Amino acid sequence
        dna_sequence (str, optional): Original DNA sequence for position mapping
        frame (int, optional): Reading frame (0, 1, or 2)
        is_reverse (bool, optional): Whether this is from the reverse complement
        min_length (int, optional): Minimum ORF length in amino acids
        
    Returns:
        list: List of dictionaries with ORF information
    """
    logger.debug(f"Finding ORFs in {'reverse' if is_reverse else 'forward'} frame {frame}")
    orfs = []
    i = 0
    
    while i < len(aa_sequence):
        # Find next start codon (M)
        start_aa_pos = i
        while start_aa_pos < len(aa_sequence) and aa_sequence[start_aa_pos] != 'M':
            start_aa_pos += 1
        
        if start_aa_pos >= len(aa_sequence):
            break  # No more start codons
        
        # Calculate DNA positions
        dna_start_pos = frame + (start_aa_pos * 3) if dna_sequence else -1
        
        # Continue until stop codon or end of sequence
        end_aa_pos = start_aa_pos + 1
        while end_aa_pos < len(aa_sequence) and aa_sequence[end_aa_pos] != '*':
            end_aa_pos += 1
        
        # Calculate DNA end position
        dna_end_pos = frame + (end_aa_pos * 3) + 2 if dna_sequence else -1
        
        # Adjust positions for reverse complement
        if is_reverse and dna_sequence:
            orig_dna_start = len(dna_sequence) - dna_end_pos - 1
            orig_dna_end = len(dna_sequence) - dna_start_pos - 1

            dna_start_pos = orig_dna_end + 1  # Add 1 to match 1-indexing
            dna_end_pos = orig_dna_start + 1  # Add 1 to match 1-indexing

            dna_start_pos, dna_end_pos = orig_dna_start, orig_dna_end
        
        # Get the ORF sequence
        orf_sequence = aa_sequence[start_aa_pos:end_aa_pos]
        if end_aa_pos < len(aa_sequence):
            orf_sequence += '*'  # Add stop codon if found
        
        # Only include ORFs that meet minimum length
        if len(orf_sequence) >= min_length:
            # Get DNA sequence of the ORF if available
            dna_sequence_of_orf = None
            if dna_sequence and dna_start_pos >= 0 and dna_end_pos >= 0:
                if not is_reverse:
                    dna_sequence_of_orf = dna_sequence[dna_start_pos:dna_end_pos+1]
                else:
                    # For reverse complement, we need to take the original segment and then reverse complement it
                    orig_segment = dna_sequence[dna_end_pos:dna_start_pos+1]
                    dna_sequence_of_orf = reverse_complement(orig_segment)
            
            orfs.append({
                'type': 'complete' if end_aa_pos < len(aa_sequence) else 'incomplete',
                'aa_start': start_aa_pos,
                'aa_end': end_aa_pos if end_aa_pos < len(aa_sequence) else len(aa_sequence) - 1,
                'dna_start': dna_start_pos,
                'dna_end': dna_end_pos,
                'length_aa': len(orf_sequence),
                'length_nt': (dna_end_pos - dna_start_pos + 1) if dna_sequence else -1,
                'sequence': orf_sequence,
                'dna_sequence': dna_sequence_of_orf
            })
        
        # Move to position after this ORF
        i = end_aa_pos + 1 if end_aa_pos < len(aa_sequence) else len(aa_sequence)
    
    logger.debug(f"Found {len(orfs)} ORFs of minimum length {min_length}")
    return orfs

def find_orfs_six_frames(dna_sequence, min_length=30, check_operons=False):
    """
    Find ORFs in all six reading frames of a DNA sequence.
    
    Args:
        dna_sequence (str): DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        
    Returns:
        dict: Dictionary with ORFs for all six frames
    """
    # Clean the sequence
    dna_sequence = dna_sequence.replace(" ", "").upper()
    
    # Get the reverse complement
    reverse_comp = reverse_complement(dna_sequence)
    
    # Store results for all frames
    results = {}
    
    # Store all unique ORFs across all frames
    all_orfs = []
    seen_aa_sequences = set()
    
    # Forward frames (5'→3')
    for frame in range(3):
        frame_name = f"5'→3' Frame {frame+1}"
        
        # Translate DNA to amino acids
        aa_sequence = translate_dna(dna_sequence, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, False, min_length)
        
        # Filter out ORFs with sequences we've seen before across all frames
        filtered_orfs = []
        for orf in orfs:
            if orf['sequence'] not in seen_aa_sequences:
                seen_aa_sequences.add(orf['sequence'])
                filtered_orfs.append(orf)
                all_orfs.append(orf)
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': filtered_orfs
        }
    
    # Reverse complement frames (3'→5')
    for frame in range(3):
        frame_name = f"3'→5' Frame {frame+1}"
        
        # Translate reverse complement to amino acids
        aa_sequence = translate_dna(reverse_comp, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, True, min_length)
        
        # Filter out ORFs with sequences we've seen before across all frames
        filtered_orfs = []
        for orf in orfs:
            if orf['sequence'] not in seen_aa_sequences:
                seen_aa_sequences.add(orf['sequence'])
                filtered_orfs.append(orf)
                all_orfs.append(orf)
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': filtered_orfs
        }
    
    return results

def print_orf_results(results, verbose=False):
    """
    Print ORF results in a tabular format.
    
    Args:
        results (dict): Results from find_orfs_six_frames
        verbose (bool, optional): Whether to include detailed information
    """
    # Collect all ORFs from all frames
    all_orfs = []
    orf_counter = 1
    
    for frame_name, frame_data in results.items():
        strand = "+" if frame_name.startswith("5'") else "-"
        frame_num = frame_name[-1]  # Extract frame number
        
        for orf in frame_data['orfs']:
            all_orfs.append({
                'label': f"ORF{orf_counter}",
                'strand': strand,
                'frame': frame_num,
                'start': orf['dna_start'] + 1,  # Convert to 1-based indexing
                'stop': orf['dna_end'] + 1,     # Convert to 1-based indexing
                'length_nt': orf['length_nt'],
                'length_aa': orf['length_aa'],
                'sequence': orf['sequence']
            })
            orf_counter += 1
    
    # Print header
    print("\nOpen Reading Frames (ORFs)")
    print("=========================\n")
    print(f"Total ORFs found: {len(all_orfs)}\n")
    
    # Print table header
    header = f"{'Label':<10}{'Strand':<10}{'Frame':<10}{'Start':<10}{'Stop':<10}{'Length (nt | aa)':<20}"
    print(header)
    print("-" * len(header))
    
    # Sort ORFs by length (descending)
    all_orfs.sort(key=lambda x: x['length_aa'], reverse=True)
    
    # Print each ORF
    for orf in all_orfs:
        print(f"{orf['label']:<10}{orf['strand']:<10}{orf['frame']:<10}{orf['start']:<10}{orf['stop']:<10}{orf['length_nt']} | {orf['length_aa']:<15}")
    
    # Print sequences if verbose
    if verbose:
        print("\nORF Sequences:")
        print("=============\n")
        for orf in all_orfs:
            print(f"{orf['label']} ({orf['length_aa']} aa):")
            # Print in blocks of 60 with position markers
            for i in range(0, len(orf['sequence']), 60):
                block = orf['sequence'][i:i+60]
                print(f"  {i+1:4d}: {block}")

def analyze_from_file(file_path, min_length=30, verbose=False, check_operons=False):
    """
    Analyze a DNA sequence from a file.
    
    Args:
        file_path (str): Path to file containing DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        verbose (bool, optional): Whether to include detailed information
    """
    try:
        # Read the sequence from file
        with open(file_path, 'r') as file:
            sequence = ""
            for line in file:
                if line.startswith('>'):  # Skip FASTA header
                    continue
                sequence += line.strip()
        
        if not sequence:
            print(f"Error: No valid sequence found in {file_path}")
            return
        
        # Perform ORF analysis
        results = find_orfs_six_frames(sequence, min_length, False)
        
        # Print results
        print_orf_results(results, verbose)
            
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"Error: {str(e)}")



def main():
    parser = argparse.ArgumentParser(description='Open Reading Frame (ORF) Finder')
    
    parser.add_argument('filepath', help='Path to file containing DNA sequence')
    parser.add_argument('--min-length', type=int, default=10, 
                      help='Minimum ORF length in amino acids (default: 10)')
    parser.add_argument('--verbose', action='store_true',
                      help='Show amino acid sequences for each ORF')
    
    args = parser.parse_args()
    
    try:
        analyze_from_file(args.filepath, args.min_length, args.verbose)
    except Exception as e:
        print(f"Error: {str(e)}")


if __name__ == "__main__":
    main()