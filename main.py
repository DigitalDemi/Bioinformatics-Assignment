import argparse
import sys
import urllib.request
import urllib.parse
import json
import time
import re
from collections import Counter

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
                'length_aa': len(orf_sequence) - 1,
                'length_nt': (dna_end_pos - dna_start_pos + 1) if dna_sequence else -1,
                'sequence': orf_sequence,
                'dna_sequence': dna_sequence_of_orf
            })
        
        # Move to position after this ORF
        i = end_aa_pos + 1 if end_aa_pos < len(aa_sequence) else len(aa_sequence)
    
    return orfs

def calculate_codon_bias(dna_sequence):
    """
    Calculate codon usage bias for a DNA sequence.
    
    Args:
        dna_sequence (str): DNA sequence
        
    Returns:
        dict: Dictionary with codon usage statistics
    """
    if not dna_sequence:
        return None
    
    # Make sure the sequence length is a multiple of 3
    dna_sequence = dna_sequence[:len(dna_sequence) - (len(dna_sequence) % 3)]
    
    # Count codons
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    codon_counts = Counter(codons)
    
    # Group by amino acid
    aa_table = {
        'I': ['ATA', 'ATC', 'ATT'],
        'M': ['ATG'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'N': ['AAC', 'AAT'],
        'K': ['AAA', 'AAG'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'P': ['CCA', 'CCC', 'CCG', 'CCT'],
        'H': ['CAC', 'CAT'],
        'Q': ['CAA', 'CAG'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'G': ['GGA', 'GGC', 'GGG', 'GGT'],
        'F': ['TTC', 'TTT'],
        'Y': ['TAC', 'TAT'],
        'C': ['TGC', 'TGT'],
        'W': ['TGG'],
        '*': ['TAA', 'TAG', 'TGA']
    }
    
    # Calculate relative usage within each amino acid group
    codon_bias = {}
    for aa, codons_list in aa_table.items():
        total_aa = sum(codon_counts[codon] for codon in codons_list)
        if total_aa > 0:
            codon_bias[aa] = {codon: codon_counts[codon] / total_aa for codon in codons_list}
    
    return {
        'codon_counts': dict(codon_counts),
        'codon_bias': codon_bias
    }

def is_prokaryotic_operon(sequence, orf_data):
    """
    Check if an ORF might be part of a prokaryotic operon based on spacing and orientation.
    
    Args:
        sequence (str): Full DNA sequence
        orf_data (dict): ORF information dictionary
        
    Returns:
        bool: True if the ORF might be part of an operon
    """
    # Look for nearby ORFs in the same orientation
    # This is a simplified heuristic - real operon detection requires more context
    
    # Criteria:
    # 1. ORF is complete (has stop codon)
    # 2. Has appropriate spacing to other ORFs (typically -4 to +15 bp)
    # 3. In same orientation as neighboring genes
    
    if orf_data['type'] != 'complete':
        return False
    
    return True  # Simplified for this example, would need genome context for real evaluation

def blast_sequence(dna_sequence, program="blastn", database="nt", max_hits=5):
    """
    Perform a BLAST search on a DNA sequence.
    Note: This function performs a remote BLAST search using NCBI's API.
    
    Args:
        dna_sequence (str): DNA sequence to BLAST
        program (str): BLAST program to use
        database (str): Database to search against
        max_hits (int): Maximum number of hits to return
        
    Returns:
        dict: BLAST results
    """
    print("Simulating BLAST search (not actually connecting to NCBI)...")
    print(f"This would search the {database} database using {program}")
    print(f"Sequence length: {len(dna_sequence)} bp")
    
    # In a real implementation, this would connect to NCBI's BLAST API
    # Since we can't make actual API calls in this context, we'll simulate results
    
    mock_results = {
        "hits": [
            {
                "id": "mock_hit_1",
                "description": "Hypothetical protein [Escherichia coli]",
                "evalue": 1e-30,
                "identity": 95.0,
                "query_cover": 98.0
            },
            {
                "id": "mock_hit_2",
                "description": "Predicted transcriptional regulator [Bacillus subtilis]",
                "evalue": 1e-15,
                "identity": 70.0,
                "query_cover": 85.0
            }
        ],
        "status": "Success",
        "message": "BLAST search completed successfully (simulated)"
    }
    
    return mock_results

def is_likely_real_orf(orf_data, blast_results=None, min_length=30, check_codon_bias=True):
    """
    Determine if an ORF is likely a real coding sequence.
    
    Args:
        orf_data (dict): ORF information dictionary
        blast_results (dict, optional): BLAST search results
        min_length (int): Minimum length in amino acids
        check_codon_bias (bool): Whether to check codon bias
        
    Returns:
        tuple: (is_real, reasons)
    """
    is_real = True
    reasons = []
    
    # Check length
    if orf_data['length_aa'] < min_length:
        is_real = False
        reasons.append(f"ORF is too short ({orf_data['length_aa']} aa, minimum {min_length})")
    
    # Check if it's complete
    if orf_data['type'] != 'complete':
        reasons.append("ORF is incomplete (no stop codon)")
        # Not disqualifying, but worth noting
    
    # Check BLAST results if available
    if blast_results and 'hits' in blast_results and blast_results['hits']:
        # If there are good BLAST hits, this is strong evidence for a real ORF
        best_hit = blast_results['hits'][0]
        if best_hit['evalue'] < 1e-10 and best_hit['identity'] > 70:
            reasons.append(f"Strong BLAST hit: {best_hit['description']} (E-value: {best_hit['evalue']})")
        else:
            reasons.append("No strong BLAST hits")
    
    # In a full implementation, we would also check:
    # - Codon bias compared to known genes in the organism
    # - Presence of ribosome binding sites
    # - GC content in the third position
    # - Whether it fits with operon structure
    
    return is_real, reasons

def find_orfs_six_frames(dna_sequence, min_length=30, check_operons=False, run_blast=False):
    """
    Find ORFs in all six reading frames of a DNA sequence.
    
    Args:
        dna_sequence (str): DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        check_operons (bool, optional): Whether to check for operon structures
        run_blast (bool, optional): Whether to run BLAST on found ORFs
        
    Returns:
        dict: Dictionary with ORFs for all six frames
    """
    # Clean the sequence
    dna_sequence = dna_sequence.replace(" ", "").upper()
    
    # Get the reverse complement
    reverse_comp = reverse_complement(dna_sequence)
    
    # Store results for all frames
    results = {}
    
    # Forward frames (5'→3')
    for frame in range(3):
        frame_name = f"5'→3' Frame {frame+1}"
        
        # Translate DNA to amino acids
        aa_sequence = translate_dna(dna_sequence, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, False, min_length)
        
        # Additional analysis for each ORF
        for i, orf in enumerate(orfs):
            # Calculate codon bias if DNA sequence is available
            if orf['dna_sequence']:
                orf['codon_bias'] = calculate_codon_bias(orf['dna_sequence'])
            
            # Check if part of prokaryotic operon
            if check_operons:
                orf['potential_operon'] = is_prokaryotic_operon(dna_sequence, orf)
            
            # Run BLAST if requested
            if run_blast and orf['dna_sequence']:
                orf['blast_results'] = blast_sequence(orf['dna_sequence'])
            
            # Evaluate if likely a real ORF
            blast_results = orf.get('blast_results') if run_blast else None
            orf['is_real'], orf['assessment_reasons'] = is_likely_real_orf(
                orf, blast_results, min_length
            )
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': orfs
        }
    
    # Reverse complement frames (3'→5')
    for frame in range(3):
        frame_name = f"3'→5' Frame {frame+1}"
        
        # Translate reverse complement to amino acids
        aa_sequence = translate_dna(reverse_comp, frame)
        
        # Find ORFs in this frame
        orfs = find_orfs_in_protein(aa_sequence, dna_sequence, frame, True, min_length)
        
        # Additional analysis for each ORF
        for i, orf in enumerate(orfs):
            # Calculate codon bias if DNA sequence is available
            if orf['dna_sequence']:
                orf['codon_bias'] = calculate_codon_bias(orf['dna_sequence'])
            
            # Check if part of prokaryotic operon
            if check_operons:
                orf['potential_operon'] = is_prokaryotic_operon(dna_sequence, orf)
            
            # Run BLAST if requested
            if run_blast and orf['dna_sequence']:
                orf['blast_results'] = blast_sequence(orf['dna_sequence'])
            
            # Evaluate if likely a real ORF
            blast_results = orf.get('blast_results') if run_blast else None
            orf['is_real'], orf['assessment_reasons'] = is_likely_real_orf(
                orf, blast_results, min_length
            )
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': orfs
        }
    
    return results

def print_orf_results(results, verbose=False, real_only=False):
    """
    Print ORF results in a readable format.
    
    Args:
        results (dict): Results from find_orfs_six_frames
        verbose (bool, optional): Whether to include detailed information
        real_only (bool, optional): Whether to only show likely real ORFs
    """
    print("\nOpen Reading Frames (ORFs) Analysis")
    print("==================================\n")
    
    # Count total and real ORFs
    total_orfs = sum(len(frame_data['orfs']) for frame_data in results.values())
    real_orfs = sum(
        sum(1 for orf in frame_data['orfs'] if orf.get('is_real', True))
        for frame_data in results.values()
    )
    
    print(f"Total ORFs found: {total_orfs}")
    if real_only:
        print(f"Likely real ORFs: {real_orfs}")
    print()
    
    for frame_name, frame_data in results.items():
        aa_sequence = frame_data['aa_sequence']
        orfs = frame_data['orfs']
        
        # Filter to only real ORFs if requested
        if real_only:
            orfs = [orf for orf in orfs if orf.get('is_real', True)]
        
        print(f"\n{frame_name}:")
        print("-" * len(frame_name))
        
        # Print ORFs
        if orfs:
            print(f"\nORFs found: {len(orfs)}")
            for i, orf in enumerate(orfs):
                print(f"\nORF {i+1}:")
                print(f"  Type: {orf['type']}")
                print(f"  Amino Acid Positions: {orf['aa_start']+1}-{orf['aa_end']+1}")
                correction_factor = 1  # adjust as needed
                print(f"  DNA Positions: {orf['dna_start']+correction_factor}-{orf['dna_end']+correction_factor}")
                print(f"  Length: {orf['length_aa']} amino acids, {orf['length_nt']} nucleotides")
                
                if 'is_real' in orf:
                    print(f"  Assessment: {'Likely real' if orf['is_real'] else 'Likely false positive'}")
                    if 'assessment_reasons' in orf and orf['assessment_reasons']:
                        print(f"  Reasons: {', '.join(orf['assessment_reasons'])}")
                
                # Print sequence in blocks
                aa_seq = orf['sequence']
                print(f"  Amino Acid Sequence: {aa_seq[:60]}" + ("..." if len(aa_seq) > 60 else ""))
                
                if verbose:
                    if orf['dna_sequence']:
                        dna_seq = orf['dna_sequence']
                        print(f"  DNA Sequence: {dna_seq[:60]}" + ("..." if len(dna_seq) > 60 else ""))
                    
                    if 'codon_bias' in orf and orf['codon_bias']:
                        print("  Codon Usage:")
                        for aa, codons in list(orf['codon_bias']['codon_bias'].items())[:3]:  # Show first 3 AAs
                            print(f"    {aa}: " + ", ".join(f"{codon}:{usage:.2f}" for codon, usage in codons.items()))
                        if len(orf['codon_bias']['codon_bias']) > 3:
                            print("    ...")
                    
                    if 'blast_results' in orf and orf['blast_results']:
                        print("  BLAST Results:")
                        for hit in orf['blast_results']['hits'][:2]:  # Show top 2 hits
                            print(f"    {hit['description']} (E-value: {hit['evalue']}, Identity: {hit['identity']}%)")
        else:
            print("\nNo ORFs found")

def analyze_sequence(sequence, min_length=30, verbose=False, check_operons=False, 
                   run_blast=False, real_only=False):
    """
    Analyze a DNA sequence for ORFs in all six reading frames.
    
    Args:
        sequence (str): DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        verbose (bool, optional): Whether to include detailed information
        check_operons (bool, optional): Whether to check for operon structures
        run_blast (bool, optional): Whether to run BLAST on found ORFs
        real_only (bool, optional): Whether to only show likely real ORFs
    """
    # Clean the sequence
    sequence = sequence.replace(" ", "").upper()
    
    print("\nSequence Information:")
    print("====================")
    print(f"Length: {len(sequence)} nucleotides")
    
    # Count nucleotides
    counts = {base: sequence.count(base) for base in 'ATGC'}
    for base, count in counts.items():
        percentage = (count / len(sequence)) * 100
        print(f"{base}: {count} ({percentage:.1f}%)")
    
    print(f"GC Content: {((counts['G'] + counts['C']) / len(sequence)) * 100:.1f}%")
    
    # Find ORFs in all six frames
    results = find_orfs_six_frames(sequence, min_length, check_operons, run_blast)
    
    # Print results
    print_orf_results(results, verbose, real_only)

def analyze_from_file(file_path, min_length=30, verbose=False, check_operons=False, 
                    run_blast=False, real_only=False):
    """
    Analyze a DNA sequence from a file.
    
    Args:
        file_path (str): Path to file containing DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        verbose (bool, optional): Whether to include detailed information
        check_operons (bool, optional): Whether to check for operon structures
        run_blast (bool, optional): Whether to run BLAST on found ORFs
        real_only (bool, optional): Whether to only show likely real ORFs
    """
    try:
        with open(file_path, 'r') as file:
            sequence = ""
            for line in file:
                # Skip header lines in FASTA format
                if line.startswith('>'):
                    continue
                sequence += line.strip()
            
            if not sequence:
                print(f"Error: No valid sequence found in {file_path}")
                return
            
            analyze_sequence(sequence, min_length, verbose, check_operons, run_blast, real_only)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        print(f"Error reading file: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Enhanced ORF Finder with False Positive Elimination')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze a DNA sequence for ORFs')
    analyze_parser.add_argument('sequence', help='DNA sequence to analyze')
    analyze_parser.add_argument('--min-length', type=int, default=30, 
                               help='Minimum ORF length in amino acids (default: 30)')
    analyze_parser.add_argument('--verbose', action='store_true',
                               help='Include detailed information in output')
    analyze_parser.add_argument('--check-operons', action='store_true',
                               help='Check for potential prokaryotic operon structures')
    analyze_parser.add_argument('--blast', action='store_true',
                               help='Run BLAST on found ORFs (simulation)')
    analyze_parser.add_argument('--real-only', action='store_true',
                               help='Only show likely real ORFs')
    
    # File command
    file_parser = subparsers.add_parser('file', help='Analyze DNA from a file')
    file_parser.add_argument('filepath', help='Path to file containing DNA sequence')
    file_parser.add_argument('--min-length', type=int, default=30, 
                            help='Minimum ORF length in amino acids (default: 30)')
    file_parser.add_argument('--verbose', action='store_true',
                            help='Include detailed information in output')
    file_parser.add_argument('--check-operons', action='store_true',
                            help='Check for potential prokaryotic operon structures')
    file_parser.add_argument('--blast', action='store_true',
                            help='Run BLAST on found ORFs (simulation)')
    file_parser.add_argument('--real-only', action='store_true',
                            help='Only show likely real ORFs')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        analyze_sequence(args.sequence, args.min_length, args.verbose, 
                        args.check_operons, args.blast, args.real_only)
    elif args.command == 'file':
        analyze_from_file(args.filepath, args.min_length, args.verbose, 
                         args.check_operons, args.blast, args.real_only)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()