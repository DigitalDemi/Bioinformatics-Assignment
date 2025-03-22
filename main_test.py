import sys
import os
import csv

csv_rows = []

compliment_nucleotide = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

CodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def read_dna_seq(file):
    """
    Reads a DNA sequence from a FASTA file.
    Skips the first line (header) and returns the sequence as a string.
    """
    with open(file, "r") as f:
        f.readline()  # Skip the header
        return f.read().replace("\n", "")

def compliment_strand(dna_seq):
    """
    Generates the complementary (reverse complement) strand of a given DNA sequence.
    """
    return "".join(compliment_nucleotide[n] for n in dna_seq[::-1])

def read_codon(dna_seq, reading_frame_no=0):
    """
    Reads codons from a DNA sequence starting at a given reading frame.
    Returns a list of tuples (position, codon) where position is the 1-indexed nucleotide index.
    """
    return [(i + 1, dna_seq[i:i+3]) for i in range(reading_frame_no, len(dna_seq), 3) if len(dna_seq[i:i+3]) == 3]

def get_possible_orfs(aminos, positions, nested=False):
    """
    Identifies possible Open Reading Frames (ORFs) from an amino acid sequence.
    Returns a list of tuples:
      (start_index, stop_index, start_nt, stop_nt, protein_sequence)
    where start_index and stop_index are 0-indexed positions in the amino acid list,
    start_nt is the nucleotide position (first nucleotide of the start codon),
    and stop_nt is the nucleotide position (third nucleotide of the stop codon).
    This function does not handle nested ORFs
    """
    start_amino = "M"  # Start codon
    stop_amino = "*"   # Stop codon
    orfs = []
    
    if nested:
        start_positions = [i for i, a in enumerate(aminos) if a == start_amino]
        stop_positions = [i for i, a in enumerate(aminos) if a == stop_amino]
        
        # Remove start positions that are before the current start and current stop position (remove nested ORFs).
        for i, start in enumerate(start_positions):
            for stop in stop_positions:
                if stop > start:
                    start_positions = start_positions[i:]
                    break



        for start in start_positions:
            for stop in stop_positions:
                if stop > start:
                    start_nt = positions[start]      # first base of start codon
                    stop_nt = positions[stop] + 2      # third base of stop codon
                    # Exclude the stop codon from the protein sequence.
                    protein_seq = "".join(aminos[start:stop+1])
                    
                    if len(protein_seq) > 10:
                        orfs.append((start, stop, start_nt, stop_nt, protein_seq))
                    
                    break  # use the first stop codon after the start
    else:
        start = None
        for i, a in enumerate(aminos):
            if a == start_amino and start is None:
                start = i  # Set the first start codon
            elif a == stop_amino and start is not None:
                start_nt = positions[start]  # First base of start codon
                stop_nt = positions[i] + 2   # Third base of stop codon
                protein_seq = "".join(aminos[start:i+1])  # Exclude stop codon
                if len(protein_seq) > 10:
                    orfs.append((start+1, i+1, start_nt, stop_nt, protein_seq))
                
                start = None  # Reset start to avoid nested ORFs

    return orfs

def complement_to_original_pos(pos, seq_length):
    """
    Converts a nucleotide position from complement coordinates (1-indexed from the 5' end of the complement)
    to the corresponding position in the original sequence.
    """
    return seq_length - pos + 1

def main(dna_seq_file):
    # Read primary DNA sequence and determine its length.
    dna_seq = read_dna_seq(dna_seq_file)
    seq_length = len(dna_seq)
    
    # Generate the reverse complement.
    comp_dna_seq = compliment_strand(dna_seq)
    
    with open("log.txt", "w") as log_file:
        log_file.write(f"\nPrimary DNA Sequence (5' to 3'):\n{dna_seq}\n\n")
        
        log_file.write("Primary Strand ORFs:\n")
        # Process all three reading frames for the primary strand.
        for frame in [0, 1, 2]:
            codons = read_codon(dna_seq, frame)
            if not codons:
                continue
            actual_codons = [c[1] for c in codons]
            # Translate each codon; skip those not in our CodonTable.
            amino_acids = [CodonTable[c] for c in actual_codons if c in CodonTable]
            positions = [c[0] for c in codons]
            orfs = get_possible_orfs(amino_acids, positions, nested=False)
            
            log_file.write(f"\nReading Frame {frame+1}:\n")
            log_file.write("Codons:\n" + " ".join(actual_codons) + "\n")
            log_file.write("Amino Acid Sequence:\n" + "".join(amino_acids) + "\n\n")
            log_file.write(f"Open Reading Frames:    Primary Strand {frame+1}\n")
            for start, stop, start_nt, stop_nt, seq in orfs:
                
                # Skip ORFs that are too short.
                if len(seq) < 10:
                    continue
                
                # For primary strand, convert 0-indexed AA positions to 1-indexed.
                aa_start = start + 1
                aa_stop = stop  # since the stop codon is not included in the protein sequence
                aa_length = stop - start
                nt_length = stop_nt - start_nt + 1
                log_file.write(f"Start (AA): {aa_start}, Stop (AA): {aa_stop}, Length: {aa_length} aa\n")
                log_file.write(f"Start (NT): {start_nt}, Stop (NT): {stop_nt}, Length: {nt_length} nt\n")
                log_file.write(f"AA Sequence: {seq}\n\n")
        
        log_file.write("\n\n\n\n\n")
        log_file.write(f"\nComplementary DNA Sequence (3' to 5'):\n{comp_dna_seq[::-1]}\n\n")

        log_file.write("\nComplementary Strand ORFs:\n")
        # Process all three reading frames for the complementary strand.
        for frame in [0, 1, 2]:
            codons = read_codon(comp_dna_seq, frame)
            if not codons:
                continue
            actual_codons = [c[1] for c in codons]
            amino_acids = [CodonTable[c] for c in actual_codons if c in CodonTable]
            positions = [c[0] for c in codons]
            orfs = get_possible_orfs(amino_acids, positions, nested=False)
            
            log_file.write(f"\nReading Frame {frame+1} (complement):\n")
            log_file.write("Codons:\n" + " ".join(actual_codons[::-1]) + "\n")
            log_file.write("Amino Acid Sequence:\n" + "".join(amino_acids[::-1]) + "\n")
            log_file.write(f"Open Reading Frames:    Complimentary Strand {frame+1}\n")
            total_codons = len(codons)  # total codons in this reading frame on the complement
            
            for start, stop, comp_start_nt, comp_stop_nt, seq in orfs:
                # Skip ORFs that are too short.
                if len(seq) < 10:
                    continue
                # Convert nucleotide coordinates from the complement to original coordinates.
                orig_stop_nt = complement_to_original_pos(comp_stop_nt, seq_length)
                orig_start_nt  = complement_to_original_pos(comp_start_nt, seq_length)
                nt_length = abs(orig_stop_nt - orig_start_nt) + 1
                # Reverse-map amino acid positions.
                # The codon at index 0 in the complement corresponds to the last codon in the original orientation.
                aa_start = total_codons - start
                aa_stop  = total_codons - stop
                aa_length = (total_codons - start) - (total_codons - stop + 1)
                # Alternatively, you can compute aa_length as stop - start.
                aa_length = stop - start
                log_file.write(f"Start (AA): {aa_start}, Stop (AA): {aa_stop}, Length: {aa_length} aa\n")
                log_file.write(f"Start (NT): {orig_start_nt}, Stop (NT): {orig_stop_nt}, Length: {nt_length} nt\n")
                log_file.write(f"AA Sequence: {seq[::-1]}\n\n")


def save_orfs_to_csv(orfs, reading_frame, strand, seq_length, full_aa_length, output_file):
    """
    Saves ORF data to a CSV file.
    """
    for start_aa, stop_aa, start_nt, stop_nt, aa_sequence in orfs:
        if strand == "-":
            # Reverse the start and stop for the complementary strand. 3' to 5'
            start_aa = full_aa_length - start_aa + 1
            stop_aa = full_aa_length - stop_aa + 1
            start_nt = complement_to_original_pos(start_nt, seq_length)
            stop_nt = complement_to_original_pos(stop_nt, seq_length)

        # Use abs to ensure positive values.
        aa_length = abs(stop_aa - start_aa)  # AA length
        nt_length = abs(stop_nt - start_nt) + 1  # NT length

        row = [start_aa, stop_aa, aa_length, start_nt, stop_nt, nt_length, reading_frame, strand, aa_sequence]
        csv_rows.append(row)

def write_to_csv(output_file):
    with open(output_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        for row in csv_rows:
            writer.writerow(row)


def save_orf_csv(dna_seq_file, output_csv):
    dna_seq = read_dna_seq(dna_seq_file)
    seq_length = len(dna_seq)
    comp_dna_seq = compliment_strand(dna_seq)
    
    # Create CSV and write header
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Start (AA)", "Stop (AA)", "AA Length", "Start (NT)", "Stop (NT)", "NT Length", "Reading Frame", "Strand", "Amino Acid Sequence"])
    
    # Process all three reading frames for the primary strand
    for frame in [0, 1, 2]:
        codons = read_codon(dna_seq, frame)
        if not codons:
            continue
        amino_acids = [CodonTable[c[1]] for c in codons if c[1] in CodonTable]
        full_aa_length = len(amino_acids)
        positions = [c[0] for c in codons]
        orfs = get_possible_orfs(amino_acids, positions, nested=False)
        save_orfs_to_csv(orfs, frame + 1, "+", seq_length, full_aa_length, output_csv)
    
    # Process all three reading frames for the complementary strand
    for frame in [0, 1, 2]:
        codons = read_codon(comp_dna_seq, frame)
        if not codons:
            continue
        amino_acids = [CodonTable[c[1]] for c in codons if c[1] in CodonTable]
        positions = [c[0] for c in codons]
        orfs = get_possible_orfs(amino_acids, positions, nested=False)
        save_orfs_to_csv(orfs, frame + 1, "-", seq_length, full_aa_length, output_csv)

    csv_rows.sort(key=lambda x: x[5], reverse=True)

    write_to_csv(output_csv)


    print(f"ORF data saved to {output_csv}")

# This program outputs two files logs with more in depth information about ORF and codons and a csv with a summary of the ORF
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Invalid arguments. Please run the program as: python assignment.py <dna_seq_file>")
        sys.exit(1)
    dna_seq_file = sys.argv[1]
    if not os.path.exists(dna_seq_file):
        print("File does not exist")
        sys.exit(1)
    if not (dna_seq_file.endswith(".fasta") or dna_seq_file.endswith(".fa")):
        print("File is not a .fasta file")
        sys.exit(1)
    main(dna_seq_file)

    save_orf_csv(dna_seq_file, "orf_data.csv")