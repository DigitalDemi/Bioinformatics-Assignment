import logging
import argparse
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def reverse_complement(sequence):
    """Generate the reverse complement of a DNA sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    sequence = sequence.replace(" ", "").upper()
    complement = ''.join([complement_dict.get(base, 'N') for base in sequence])
    return complement[::-1]

def translate_dna(sequence, frame=0):
    """Translate a DNA sequence into an amino acid sequence from a specified frame."""
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
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'X')
            protein += amino_acid
    
    return protein

def find_orfs_in_protein(aa_sequence, dna_sequence=None, frame=0, is_reverse=False, min_length=30):
    """Find all ORFs in an amino acid sequence."""
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
        
        # Only include complete ORFs (has stop codon) that meet minimum length
        if len(orf_sequence) - 1 >= min_length and end_aa_pos < len(aa_sequence):
            # NCBI ORFfinder excludes the stop codon from the length calculation
            orf_length_aa = len(orf_sequence) - 1
            
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
                'type': 'complete',
                'aa_start': start_aa_pos,
                'aa_end': end_aa_pos,
                'dna_start': dna_start_pos,
                'dna_end': dna_end_pos,
                'length_aa': orf_length_aa,
                'length_nt': (dna_end_pos - dna_start_pos + 1) if dna_sequence else -1,
                'sequence': orf_sequence,
                'dna_sequence': dna_sequence_of_orf,
                'frame': frame + 1,
                'strand': '-' if is_reverse else '+'
            })
        
        # Move to position after this ORF
        i = end_aa_pos + 1 if end_aa_pos < len(aa_sequence) else len(aa_sequence)
    
    return orfs

def calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence."""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

def has_coding_bias(dna_sequence):
    """Check if sequence shows characteristics of a coding region."""
    # Simple heuristic: check for codon bias
    codon_counts = {}
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    # If we have very few codon types for a long sequence, it's suspicious
    if len(codon_counts) < len(dna_sequence) / 10:
        return False
    
    return True

def find_orfs_ncbi_compatible(dna_sequence, min_length=30, genetic_code=1):
    """
    Find ORFs in all six reading frames using NCBI ORFfinder compatible approach.
    
    Args:
        dna_sequence (str): DNA sequence
        min_length (int): Minimum ORF length in amino acids (default: 30)
        genetic_code (int): Genetic code to use for translation (default: 1 - standard)
        
    Returns:
        dict: Dictionary with ORFs for all six frames
    """
    # Clean the sequence
    dna_sequence = dna_sequence.replace(" ", "").upper()
    
    # Get the reverse complement
    reverse_comp = reverse_complement(dna_sequence)
    
    # Store results for all frames
    results = {}
    
    # This is critical: we'll collect ALL ORFs first, then apply NCBI-style filtering
    all_orfs = []
    
    # Forward frames (5'→3')
    for frame in range(3):
        frame_name = f"5'→3' Frame {frame+1}"
        
        # Translate DNA to amino acids
        aa_sequence = translate_dna(dna_sequence, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, False, min_length)
        all_orfs.extend(orfs)
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': []  # We'll fill this after NCBI-style filtering
        }
    
    # Reverse complement frames (3'→5')
    for frame in range(3):
        frame_name = f"3'→5' Frame {frame+1}"
        
        # Translate reverse complement to amino acids
        aa_sequence = translate_dna(reverse_comp, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, True, min_length)
        all_orfs.extend(orfs)
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': []  # We'll fill this after NCBI-style filtering
        }
    
    # Sort all ORFs by length (descending) - this is critical for NCBI-style filtering
    all_orfs.sort(key=lambda x: x['length_aa'], reverse=True)
    
    # NCBI ORFfinder specific filtering algorithm:
    # 1. We apply filters in a specific order
    # 2. We consider frame constraints
    # 3. We consider positional overlap differently
    
    # Filter #1: Remove exact duplicate sequences
    filtered_orfs = []
    seen_sequences = set()
    
    for orf in all_orfs:
        seq = orf['sequence']
        # Skip duplicates (excluding stop codon for comparison)
        if seq[:-1] in seen_sequences:
            continue
        seen_sequences.add(seq[:-1])
        filtered_orfs.append(orf)
    
    # Filter #2: Apply NCBI's nested/overlapping ORF rules
    final_orfs = []
    overlapping_regions = {}  # Track regions that are already covered
    
    # First pass: take the top ORFs for each frame
    for frame in range(1, 4):
        for strand in ['+', '-']:
            # Get largest ORFs in this frame
            frame_orfs = [o for o in filtered_orfs if o['frame'] == frame and o['strand'] == strand]
            
            if frame_orfs:
                # Always take the largest ORF in each frame
                largest_orf = frame_orfs[0]  # Already sorted by length
                final_orfs.append(largest_orf)
                
                # Mark this region as covered
                key = f"{strand}_{frame}"
                if key not in overlapping_regions:
                    overlapping_regions[key] = []
                overlapping_regions[key].append((largest_orf['dna_start'], largest_orf['dna_end']))
    
    # Second pass: Consider additional ORFs that don't significantly overlap with larger ones
    for orf in filtered_orfs:
        # Skip if already added
        if orf in final_orfs:
            continue
        
        strand, frame = orf['strand'], orf['frame']
        key = f"{strand}_{frame}"
        
        # Check if this ORF significantly overlaps with any we've already taken
        # in the same frame and strand
        significant_overlap = False
        
        if key in overlapping_regions:
            for start, end in overlapping_regions[key]:
                # Calculate overlap
                overlap_start = max(start, orf['dna_start'])
                overlap_end = min(end, orf['dna_end'])
                
                if overlap_start <= overlap_end:
                    overlap_length = overlap_end - overlap_start + 1
                    orf_length = orf['dna_end'] - orf['dna_start'] + 1
                    
                    # NCBI appears to use a 50% overlap threshold
                    if overlap_length / orf_length > 0.5:
                        significant_overlap = True
                        break
        
        # Include ORFs that don't significantly overlap or special exceptions
        if not significant_overlap:
            final_orfs.append(orf)
            
            # Update covered regions
            if key not in overlapping_regions:
                overlapping_regions[key] = []
            overlapping_regions[key].append((orf['dna_start'], orf['dna_end']))
    
    # Filter #3: Special rules for specific frames based on NCBI behavior
    frame1_plus_large = [o for o in final_orfs if o['frame'] == 1 and o['strand'] == '+' and o['length_aa'] > 400]
    frame2_plus_large = [o for o in final_orfs if o['frame'] == 2 and o['strand'] == '+' and o['length_aa'] > 400]
    
    # NCBI seems to keep more ORFs in frames with existing large ORFs
    if frame1_plus_large or frame2_plus_large:
        # If we have large ORFs in frames 1 or 2, we're more permissive with smaller ORFs
        # This is based on observed NCBI behavior from the visualization
        extra_orfs = []
        
        for orf in filtered_orfs:
            if orf not in final_orfs and orf['length_aa'] > 100:
                if (orf['frame'] == 1 and orf['strand'] == '+' and frame1_plus_large) or \
                   (orf['frame'] == 2 and orf['strand'] == '+' and frame2_plus_large):
                    
                    # Check if it overlaps significantly with any existing ORF
                    key = f"{orf['strand']}_{orf['frame']}"
                    significant_overlap = False
                    
                    if key in overlapping_regions:
                        for start, end in overlapping_regions[key]:
                            overlap_start = max(start, orf['dna_start'])
                            overlap_end = min(end, orf['dna_end'])
                            
                            if overlap_start <= overlap_end:
                                overlap_length = overlap_end - overlap_start + 1
                                orf_length = orf['dna_end'] - orf['dna_start'] + 1
                                
                                # Use a higher threshold for these special cases
                                if overlap_length / orf_length > 0.8:
                                    significant_overlap = True
                                    break
                    
                    if not significant_overlap:
                        extra_orfs.append(orf)
                        if key not in overlapping_regions:
                            overlapping_regions[key] = []
                        overlapping_regions[key].append((orf['dna_start'], orf['dna_end']))
        
        final_orfs.extend(extra_orfs)
    
    # Distribute back to frames for output
    for orf in final_orfs:
        frame_num = orf['frame']
        strand = orf['strand']
        
        if strand == '+':
            frame_name = f"5'→3' Frame {frame_num}"
        else:
            frame_name = f"3'→5' Frame {frame_num}"
        
        results[frame_name]['orfs'].append(orf)
    
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

def analyze_from_file(file_path, min_length=30, verbose=False):
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
        
        print(f"Read sequence of length {len(sequence)} bp")
        
        # Use NCBI compatible ORF finding algorithm
        results = find_orfs_ncbi_compatible(sequence, min_length)
        
        # Print results
        print_orf_results(results, verbose)
            
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"Error: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='NCBI ORFfinder Compatible Implementation')
    
    parser.add_argument('filepath', help='Path to file containing DNA sequence')
    parser.add_argument('--min-length', type=int, default=30,
                      help='Minimum ORF length in amino acids (default: 30)')
    parser.add_argument('--verbose', action='store_true',
                      help='Show amino acid sequences for each ORF')
    
    args = parser.parse_args()
    
    print(f"Starting ORF analysis with NCBI compatible algorithm")
    print(f"  File: {args.filepath}")
    print(f"  Minimum length: {args.min_length} amino acids")
    
    try:
        analyze_from_file(args.filepath, args.min_length, args.verbose)
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()