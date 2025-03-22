import argparse
from collections import Counter
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

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
    logger.debug(f"Translating DNA sequence of length {len(sequence)} from frame {frame}")
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
    
    logger.debug(f"Translation complete. Amino acid sequence length: {len(protein)}")
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
        check_operons (bool, optional): Whether to check for operon structures
        
    Returns:
        dict: Dictionary with ORFs for all six frames
    """
    logger.info(f"Analyzing sequence of length {len(dna_sequence)} for ORFs in all six frames")
    # Clean the sequence
    dna_sequence = dna_sequence.replace(" ", "").upper()
    
    # Get the reverse complement
    reverse_comp = reverse_complement(dna_sequence)
    logger.debug("Generated reverse complement sequence")
    
    # Store results for all frames
    results = {}
    
    # Store all unique ORFs across all frames
    all_orfs = []
    seen_aa_sequences = set()
    
    # Forward frames (5'→3')
    for frame in range(3):
        frame_name = f"5'→3' Frame {frame+1}"
        logger.info(f"Analyzing {frame_name}")
        
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
        
        logger.info(f"Found {len(filtered_orfs)} unique ORFs in {frame_name}")
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': filtered_orfs
        }
    
    # Reverse complement frames (3'→5')
    for frame in range(3):
        frame_name = f"3'→5' Frame {frame+1}"
        logger.info(f"Analyzing {frame_name}")
        
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
        
        logger.info(f"Found {len(filtered_orfs)} unique ORFs in {frame_name}")
        
        results[frame_name] = {
            'aa_sequence': aa_sequence,
            'orfs': filtered_orfs
        }
    
    total_orfs = len(all_orfs)
    logger.info(f"Analysis complete. Total unique ORFs found across all frames: {total_orfs}")
    return results

def print_orf_results(results, verbose=False):
    """
    Print ORF results in a readable format.
    
    Args:
        results (dict): Results from find_orfs_six_frames
        verbose (bool, optional): Whether to include detailed information
    """
    print("\nOpen Reading Frames (ORFs) Analysis")
    print("==================================\n")
    
    # Count total ORFs
    total_orfs = sum(len(frame_data['orfs']) for frame_data in results.values())
    
    print(f"Total ORFs found: {total_orfs}")
    print()
    
    for frame_name, frame_data in results.items():
        aa_sequence = frame_data['aa_sequence']
        orfs = frame_data['orfs']
        
        print(f"\n{frame_name}:")
        print("-" * len(frame_name))
        
        # Print complete amino acid sequence for the frame
        print(f"\nComplete Amino Acid Sequence:")
        # Print in blocks of 60 characters with position markers
        for i in range(0, len(aa_sequence), 60):
            block = aa_sequence[i:i+60]
            print(f"{i+1:4d}: {block}")
        
        # Print ORFs
        if orfs:
            print(f"\nORFs found: {len(orfs)}")
            for i, orf in enumerate(orfs):
                print(f"\nORF {i+1}:")
                print(f"  Type: {orf['type']}")
                print(f"  Amino Acid Positions: {orf['aa_start']+1}-{orf['aa_end']+1}")
                
                # Handle DNA positions (may be -1 for some frames)
                dna_start = orf['dna_start']
                dna_end = orf['dna_end']
                if dna_start >= 0 and dna_end >= 0:
                    correction_factor = 1  # adjust as needed for 1-based indexing
                    print(f"  DNA Positions: {dna_start+correction_factor}-{dna_end+correction_factor}")
                    print(f"  Length: {orf['length_aa']} amino acids, {orf['length_nt']} nucleotides")
                else:
                    print(f"  Length: {orf['length_aa']} amino acids")
                
                # Print complete amino acid sequence with formatting
                aa_seq = orf['sequence']
                print(f"  Amino Acid Sequence:")
                # Print in blocks of 60 with position markers relative to the ORF
                for j in range(0, len(aa_seq), 60):
                    block = aa_seq[j:j+60]
                    print(f"    {j+1:4d}: {block}")
                
                # Print DNA sequence if verbose and available
                if verbose and orf.get('dna_sequence'):
                    dna_seq = orf['dna_sequence']
                    print(f"  DNA Sequence:")
                    # Print DNA sequence in blocks of 60
                    for j in range(0, len(dna_seq), 60):
                        block = dna_seq[j:j+60]
                        print(f"    {j+1:4d}: {block}")
        else:
            print("\nNo ORFs found")

def analyze_from_file(file_path, min_length=30, verbose=False, check_operons=False):
    """
    Analyze a DNA sequence from a file.
    
    Args:
        file_path (str): Path to file containing DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        verbose (bool, optional): Whether to include detailed information
        check_operons (bool, optional): Whether to check for operon structures
    """
    logger.info(f"Reading sequence from file: {file_path}")
    try:
        with open(file_path, 'r') as file:
            sequence = ""
            is_fasta = False
            for line in file:
                # Skip header lines in FASTA format
                if line.startswith('>'):
                    is_fasta = True
                    logger.info(f"Detected FASTA format: {line.strip()}")
                    continue
                sequence += line.strip()
            
            if not sequence:
                logger.error(f"No valid sequence found in {file_path}")
                print(f"Error: No valid sequence found in {file_path}")
                return
            
            logger.info(f"Successfully read sequence of length {len(sequence)}" + 
                      (f" in FASTA format" if is_fasta else ""))
            
            # Perform ORF analysis
            logger.debug("Starting ORF analysis...")
            results = find_orfs_six_frames(sequence, min_length, check_operons)
            logger.debug(f"Analysis complete with {sum(len(frame_data['orfs']) for frame_data in results.values())} total ORFs")
            
            # Print results
            logger.debug("Printing results...")
            print_orf_results(results, verbose)
            logger.debug("Finished printing results")
            
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        logger.error(f"Error reading file: {str(e)}", exc_info=True)
        print(f"Error reading file: {str(e)}")


def main():
    parser = argparse.ArgumentParser(description='Enhanced ORF Finder')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # File command
    file_parser = subparsers.add_parser('file', help='Analyze DNA from a file')
    file_parser.add_argument('filepath', help='Path to file containing DNA sequence')
    file_parser.add_argument('--min-length', type=int, default=30, 
                            help='Minimum ORF length in amino acids (default: 30)')
    file_parser.add_argument('--verbose', action='store_true',
                            help='Include detailed information in output')
    file_parser.add_argument('--check-operons', action='store_true',
                            help='Check for potential prokaryotic operon structures')
    file_parser.add_argument('--debug', action='store_true',
                            help='Enable debug logging')
    
    # Test command
    test_parser = subparsers.add_parser('test', help='Run snapshot tests on a directory of DNA sequences')
    test_parser.add_argument('test_dir', help='Directory containing test FASTA files')
    test_parser.add_argument('--snapshot-dir', default='./snapshots', 
                            help='Directory to store/load snapshots (default: ./snapshots)')
    test_parser.add_argument('--min-length', type=int, default=30,
                            help='Minimum ORF length in amino acids (default: 30)')
    test_parser.add_argument('--debug', action='store_true',
                            help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Enable debug logging if requested
    if hasattr(args, 'debug') and args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug logging enabled")
    
    logger.info("Starting ORF analysis")
   
    if args.command == 'file':
        try:
            logger.debug(f"Analyzing file: {args.filepath} with min_length={args.min_length}")
            analyze_from_file(args.filepath, args.min_length, args.verbose, 
                             args.check_operons)
            logger.debug("Analysis completed successfully")
        except Exception as e:
            logger.error(f"Error during analysis: {str(e)}", exc_info=True)
            print(f"Error during analysis: {str(e)}")
    else:
        parser.print_help()
        logger.info("No command provided, displaying help")

if __name__ == "__main__":
    main()