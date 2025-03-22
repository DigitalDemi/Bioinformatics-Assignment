import argparse
import sys
import re
import logging
import os
import io
import difflib
import hashlib
from datetime import datetime

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

def sequence_similarity(seq1, seq2, threshold=0.8):
    """
    Check if two sequences are similar beyond a threshold.
    Returns True if sequences are highly similar.
    """
    if len(seq1) != len(seq2):
        return False
    
    matches = sum(a == b for a, b in zip(seq1, seq2))
    similarity = matches / len(seq1)
    return similarity >= threshold

def find_orfs_six_frames(dna_sequence, min_length=30):
    """
    Find ORFs in all six reading frames of a DNA sequence.
    
    Args:
        dna_sequence (str): DNA sequence
        min_length (int, optional): Minimum ORF length in amino acids
        
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
            # Check if we've seen this sequence or a highly similar one
            is_similar = False
            for seen_seq in seen_aa_sequences:
                if len(orf['sequence']) == len(seen_seq) and sequence_similarity(orf['sequence'], seen_seq):
                    is_similar = True
                    break
            
            if not is_similar:
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
            # Check if we've seen this sequence or a highly similar one
            is_similar = False
            for seen_seq in seen_aa_sequences:
                if len(orf['sequence']) == len(seen_seq) and sequence_similarity(orf['sequence'], seen_seq):
                    is_similar = True
                    break
            
            if not is_similar:
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
    
    logger.debug(f"Preparing to print {len(all_orfs)} ORFs")
    
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
        logger.debug("Printing detailed ORF sequences (verbose mode)")
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
            
            # Clean the sequence
            sequence = sequence.replace(" ", "").upper()
            
            # Find ORFs in all six frames
            results = find_orfs_six_frames(sequence, min_length)
            
            # Print results
            print_orf_results(results, verbose)
            logger.info("Analysis completed successfully")
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        print(f"Error: File {file_path} not found.")
    except Exception as e:
        logger.error(f"Error reading file: {str(e)}", exc_info=True)
        print(f"Error reading file: {str(e)}")

def run_snapshot_tests(test_data_dir, snapshot_dir, min_length=30):
    """
    Run snapshot tests on a collection of DNA sequences.
    
    Args:
        test_data_dir (str): Directory containing test FASTA files
        snapshot_dir (str): Directory to store/load snapshots
        min_length (int): Minimum ORF length in amino acids
    """
    # Create directories if they don't exist
    os.makedirs(snapshot_dir, exist_ok=True)
    
    # Find all test files
    if not os.path.exists(test_data_dir):
        logger.error(f"Test data directory not found: {test_data_dir}")
        print(f"Error: Test data directory not found: {test_data_dir}")
        return
        
    test_files = [f for f in os.listdir(test_data_dir) if f.lower().endswith(('.fasta', '.fa'))]
    
    if not test_files:
        logger.error(f"No test files found in {test_data_dir}")
        print(f"Error: No test files found in {test_data_dir}")
        return
        
    # Track test results
    passed = 0
    failed = 0
    first_run = 0
    
    print(f"\nRunning snapshot tests from {test_data_dir}")
    print(f"Snapshot directory: {snapshot_dir}")
    print(f"Found {len(test_files)} test files")
    print("=" * 60)
    
    # Record start time
    start_time = datetime.now()
    
    for test_file in test_files:
        test_path = os.path.join(test_data_dir, test_file)
        print(f"\nTesting {test_file}...", end="")
        sys.stdout.flush()
        
        # Generate a unique snapshot name based on file content and parameters
        with open(test_path, 'rb') as f:
            content = f.read()
            file_hash = hashlib.md5(content).hexdigest()[:10]
        
        param_hash = hashlib.md5(f"min_length={min_length}".encode()).hexdigest()[:6]
        snapshot_name = f"{os.path.splitext(test_file)[0]}_{file_hash}_{param_hash}.snapshot"
        snapshot_path = os.path.join(snapshot_dir, snapshot_name)
        
        # Capture stdout to get the output
        original_stdout = sys.stdout
        captured_output = io.StringIO()
        sys.stdout = captured_output
        
        try:
            # Run analysis
            analyze_from_file(test_path, min_length, verbose=False)
            # Get captured output
            sys.stdout = original_stdout
            current_output = captured_output.getvalue()
            
            # Check if snapshot exists
            if os.path.exists(snapshot_path):
                # Compare with existing snapshot
                with open(snapshot_path, 'r') as f:
                    expected_output = f.read()
                
                if current_output == expected_output:
                    print(" ✅ Passed")
                    passed += 1
                else:
                    print(" ❌ Failed")
                    failed += 1
                    
                    # Show diff
                    print("\nDifferences found:")
                    diff = difflib.unified_diff(
                        expected_output.splitlines(),
                        current_output.splitlines(),
                        fromfile='Expected',
                        tofile='Actual',
                        lineterm=''
                    )
                    # Show only first 10 lines of diff to avoid excessive output
                    diff_lines = list(diff)
                    for line in diff_lines[:10]:
                        print(line)
                    if len(diff_lines) > 10:
                        print(f"...and {len(diff_lines) - 10} more differences")
            else:
                # First run, create snapshot
                with open(snapshot_path, 'w') as f:
                    f.write(current_output)
                print(" 📸 Created new snapshot")
                first_run += 1
                
        except Exception as e:
            sys.stdout = original_stdout
            print(f" ❌ Error: {str(e)}")
            logger.error(f"Error running test on {test_file}: {str(e)}", exc_info=True)
            failed += 1
    
    # Calculate run time
    total_time = (datetime.now() - start_time).total_seconds()
    
    # Print summary
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed, {first_run} new snapshots")
    print(f"Total time: {total_time:.2f} seconds")
    
    # Return test results
    return {
        'passed': passed,
        'failed': failed,
        'new_snapshots': first_run,
        'total': passed + failed + first_run,
        'runtime': total_time
    }

def main():
    parser = argparse.ArgumentParser(description='ORF Finder')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # File command
    file_parser = subparsers.add_parser('file', help='Analyze DNA from a file')
    file_parser.add_argument('filepath', help='Path to file containing DNA sequence')
    file_parser.add_argument('--min-length', type=int, default=10, 
                           help='Minimum ORF length in amino acids (default: 10)')
    file_parser.add_argument('--verbose', action='store_true',
                           help='Show amino acid sequences for each ORF')
    file_parser.add_argument('--debug', action='store_true',
                           help='Enable debug logging')
    
    # Test command
    test_parser = subparsers.add_parser('test', help='Run snapshot tests on a directory of DNA sequences')
    test_parser.add_argument('test_dir', help='Directory containing test FASTA files')
    test_parser.add_argument('--snapshot-dir', default='./snapshots', 
                           help='Directory to store/load snapshots (default: ./snapshots)')
    test_parser.add_argument('--min-length', type=int, default=10,
                           help='Minimum ORF length in amino acids (default: 10)')
    test_parser.add_argument('--debug', action='store_true',
                           help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Enable debug logging if requested
    if hasattr(args, 'debug') and args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug logging enabled")
    
    logger.info("Starting ORF analysis")
    
    if args.command == 'file':
        analyze_from_file(args.filepath, args.min_length, args.verbose)
    elif args.command == 'test':
        run_snapshot_tests(args.test_dir, args.snapshot_dir, args.min_length)
    else:
        parser.print_help()
        logger.info("No command provided, displaying help")

if __name__ == "__main__":
    main()