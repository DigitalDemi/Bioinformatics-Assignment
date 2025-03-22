import argparse
import sys
import urllib.request
import urllib.parse
import json
import time
import re
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
    With debug statements to trace reverse complement Frame 3 specifically.
    """
    # print(f"DEBUG: Finding ORFs in {'reverse' if is_reverse else 'forward'} frame {frame+1}")
    
    # IMPORTANT: Check if this is Frame 3 reverse complement
    # is_frame3_reverse = is_reverse and frame == 2
    
    # if is_frame3_reverse:
    #     print(f"DEBUG FRAME3-REV: AA sequence length: {len(aa_sequence)}")
    #     print(f"DEBUG FRAME3-REV: DNA sequence length: {len(dna_sequence)}")
    #     # Print first 20 AA to verify translation
    #     print(f"DEBUG FRAME3-REV: First 20 AA: {aa_sequence[:20]}")
    
    orfs = []
    i = 0
    
    while i < len(aa_sequence):
        # Find next start codon (M)
        start_aa_pos = i
        while start_aa_pos < len(aa_sequence) and aa_sequence[start_aa_pos] != 'M':
            start_aa_pos += 1
        
        if start_aa_pos >= len(aa_sequence):
            break  # No more start codons
        
        # Continue until stop codon or end of sequence
        end_aa_pos = start_aa_pos + 1
        while end_aa_pos < len(aa_sequence) and aa_sequence[end_aa_pos] != '*':
            end_aa_pos += 1
        
        # Calculate DNA positions
        dna_start_pos = frame + (start_aa_pos * 3) if dna_sequence else -1
        dna_end_pos = frame + (end_aa_pos * 3) + 2 if dna_sequence else -1
        
        # if is_frame3_reverse:
        #     print(f"DEBUG FRAME3-REV: Found ORF at AA positions {start_aa_pos}-{end_aa_pos}")
        #     print(f"DEBUG FRAME3-REV: Initial DNA positions: {dna_start_pos}-{dna_end_pos}")
        
        # Adjust positions for reverse complement
        if is_reverse and dna_sequence:
            orig_dna_start = len(dna_sequence) - dna_end_pos - 1
            orig_dna_end = len(dna_sequence) - dna_start_pos - 1
            
            # if is_frame3_reverse:
            #     print(f"DEBUG FRAME3-REV: Reverse position calculation:")
            #     print(f"DEBUG FRAME3-REV:   len(dna) = {len(dna_sequence)}")
            #     print(f"DEBUG FRAME3-REV:   orig_start = {len(dna_sequence)} - {dna_end_pos} - 1 = {orig_dna_start}")
            #     print(f"DEBUG FRAME3-REV:   orig_end = {len(dna_sequence)} - {dna_start_pos} - 1 = {orig_dna_end}")
                
            #     # Check if this matches ORF7 coordinates (405-310)
            #     if (abs(orig_dna_start - 310) <= 1) and (abs(orig_dna_end - 405) <= 1):
            #         print(f"DEBUG FRAME3-REV: FOUND MATCH FOR ORF7!")
            
            # dna_start_pos, dna_end_pos = orig_dna_start, orig_dna_end
            
            # if is_frame3_reverse:
            #     print(f"DEBUG FRAME3-REV: Final DNA positions: {dna_start_pos}-{dna_end_pos}")
        
        # Get the ORF sequence
        orf_sequence = aa_sequence[start_aa_pos:end_aa_pos]
        if end_aa_pos < len(aa_sequence):
            orf_sequence += '*'  # Add stop codon if found
        
        # if is_frame3_reverse:
        #     print(f"DEBUG FRAME3-REV: ORF sequence: {orf_sequence}")
        #     print(f"DEBUG FRAME3-REV: Length: {len(orf_sequence)} AA")
        
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
            
            # if is_frame3_reverse:
            #     print(f"DEBUG FRAME3-REV: Added ORF to results ({dna_start_pos}-{dna_end_pos})")
        # else:
        #     if is_frame3_reverse:
        #         print(f"DEBUG FRAME3-REV: ORF too short ({len(orf_sequence)} < {min_length}), skipping")
        
        # Move to position after this ORF
        i = end_aa_pos + 1 if end_aa_pos < len(aa_sequence) else len(aa_sequence)
    
    # if is_frame3_reverse:
    #     print(f"DEBUG FRAME3-REV: Total ORFs found in Frame 3 reverse: {len(orfs)}")
    #     for idx, orf in enumerate(orfs):
    #         print(f"DEBUG FRAME3-REV:   ORF {idx+1}: {orf['dna_start']}-{orf['dna_end']}, {orf['length_aa']} AA")
    
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
        'I': ['ATA', 'AT', 'ATT'],
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
    logger.info(f"Simulating BLAST search using {program} against {database}")
    logger.info(f"Sequence length: {len(dna_sequence)} bp")
    
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
        
        # Additional analysis for each ORF
        for i, orf in enumerate(filtered_orfs):
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
        
        # Additional analysis for each ORF
        for i, orf in enumerate(filtered_orfs):
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
            'orfs': filtered_orfs
        }
    
    total_orfs = len(all_orfs)
    logger.info(f"Analysis complete. Total unique ORFs found across all frames: {total_orfs}")
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
                correction_factor = 1  # adjust as needed
                print(f"  DNA Positions: {orf['dna_start']+correction_factor}-{orf['dna_end']+correction_factor}")
                print(f"  Length: {orf['length_aa']} amino acids, {orf['length_nt']} nucleotides")
                
                if 'is_real' in orf:
                    print(f"  Assessment: {'Likely real' if orf['is_real'] else 'Likely false positive'}")
                    if 'assessment_reasons' in orf and orf['assessment_reasons']:
                        print(f"  Reasons: {', '.join(orf['assessment_reasons'])}")
                
                # Print complete amino acid sequence with formatting
                aa_seq = orf['sequence']
                print(f"  Amino Acid Sequence:")
                # Print in blocks of 60 with position markers relative to the ORF
                for j in range(0, len(aa_seq), 60):
                    block = aa_seq[j:j+60]
                    print(f"    {j+1:4d}: {block}")
                
                if verbose:
                    if orf['dna_sequence']:
                        dna_seq = orf['dna_sequence']
                        print(f"  DNA Sequence:")
                        # Print DNA sequence in blocks of 60
                        for j in range(0, len(dna_seq), 60):
                            block = dna_seq[j:j+60]
                            print(f"    {j+1:4d}: {block}")
                    
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
    
    logger.info(f"Starting analysis of sequence with length {len(sequence)}")
    
    print("\nSequence Information:")
    print("====================")
    print(f"Length: {len(sequence)} nucleotides")
    
    # Count nucleotides
    counts = {base: sequence.count(base) for base in 'ATGC'}
    for base, count in counts.items():
        percentage = (count / len(sequence)) * 100
        print(f"{base}: {count} ({percentage:.1f}%)")
    
    gc_content = ((counts['G'] + counts['C']) / len(sequence)) * 100
    print(f"GC Content: {gc_content:.1f}%")
    logger.info(f"Sequence GC content: {gc_content:.1f}%")
    
    # Find ORFs in all six frames
    results = find_orfs_six_frames(sequence, min_length, check_operons, run_blast)
    
    # Print results
    print_orf_results(results, verbose, real_only)
    logger.info("Analysis completed successfully")

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
            
            analyze_sequence(sequence, min_length, verbose, check_operons, run_blast, real_only)
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
    import os
    import io
    import sys
    import difflib
    import hashlib
    from datetime import datetime
    
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
            analyze_from_file(test_path, min_length, verbose=False, check_operons=False, 
                             run_blast=False, real_only=False)
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
    analyze_parser.add_argument('--debug', action='store_true',
                               help='Enable debug logging')
    
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
    
    if args.command == 'analyze':
        analyze_sequence(args.sequence, args.min_length, args.verbose, 
                        args.check_operons, args.blast, args.real_only)
    elif args.command == 'file':
        analyze_from_file(args.filepath, args.min_length, args.verbose, 
                         args.check_operons, args.blast, args.real_only)
    elif args.command == 'test':
        run_snapshot_tests(args.test_dir, args.snapshot_dir, args.min_length)
    else:
        parser.print_help()
        logger.info("No command provided, displaying help")

if __name__ == "__main__":
    main()