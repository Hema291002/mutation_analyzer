"""
Mutation Impact Analyzer - Fixed Version
Correctly handles stop codons in mutation classification
"""

from Bio.Seq import Seq
from Bio import SeqIO
import requests
import json
from datetime import datetime

# Genetic code dictionary with stop codons explicitly marked
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * represents STOP
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',  # * represents STOP
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Stop codons set for easy checking
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Amino acid names (including STOP)
AA_NAMES = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic acid', 'E': 'Glutamic acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine',
    '*': 'STOP'
}


def is_stop_codon(codon):
    """
    Check if a codon is a stop codon.
    
    Args:
        codon (str): Three-nucleotide DNA sequence
        
    Returns:
        bool: True if codon is a stop codon (TAA, TAG, TGA)
    """
    return codon.upper() in STOP_CODONS


def translate_codon(codon):
    """
    Translate a DNA codon to its corresponding amino acid.
    Now properly handles stop codons.
    
    Args:
        codon (str): Three-nucleotide DNA sequence
        
    Returns:
        str: Single-letter amino acid code (or '*' for stop codon)
        
    Raises:
        ValueError: If codon is not exactly 3 nucleotides or contains invalid characters
    """
    if len(codon) != 3:
        raise ValueError(f"Codon must be exactly 3 nucleotides, got {len(codon)}")
    
    codon = codon.upper()
    
    # Validate nucleotides
    valid_nucleotides = set('ATCG')
    if not all(n in valid_nucleotides for n in codon):
        raise ValueError(f"Codon contains invalid nucleotides: {codon}")
    
    if codon not in GENETIC_CODE:
        raise ValueError(f"Unknown codon: {codon}")
    
    return GENETIC_CODE[codon]


def validate_sequence(sequence):
    """
    Validate DNA sequence.
    Now allows stop codons in the sequence.
    
    Args:
        sequence (str): DNA sequence to validate
        
    Returns:
        tuple: (is_valid (bool), error_message (str or None))
    """
    # Remove whitespace
    sequence = sequence.replace(" ", "").replace("\n", "").replace("\r", "")
    
    # Check if empty
    if not sequence:
        return False, "Sequence is empty"
    
    # Check for valid nucleotides only
    valid_nucleotides = set('ATCGatcg')
    if not all(n in valid_nucleotides for n in sequence):
        invalid_chars = set(sequence) - valid_nucleotides
        return False, f"Sequence contains invalid characters: {invalid_chars}"
    
    # Check if length is divisible by 3
    if len(sequence) % 3 != 0:
        return False, f"Sequence length ({len(sequence)}) is not divisible by 3. Incomplete codons detected."
    
    return True, None


def classify_mutation(original_codon, mutated_codon):
    """
    Classify mutation type with proper stop codon handling.
    
    Args:
        original_codon (str): Original three-nucleotide sequence
        mutated_codon (str): Mutated three-nucleotide sequence
        
    Returns:
        dict: Classification information including type, amino acids, and description
    """
    original_codon = original_codon.upper()
    mutated_codon = mutated_codon.upper()
    
    # Translate codons
    original_aa = translate_codon(original_codon)
    mutated_aa = translate_codon(mutated_codon)
    
    # Get amino acid names
    original_aa_name = AA_NAMES.get(original_aa, "Unknown")
    mutated_aa_name = AA_NAMES.get(mutated_aa, "Unknown")
    
    # Classify mutation type
    if original_aa == mutated_aa:
        # Silent mutation
        mutation_type = "Silent"
        description = f"Synonymous mutation - both codons code for {original_aa_name}"
        impact = "LOW"
        
    elif is_stop_codon(mutated_codon):
        # Nonsense mutation (creates stop codon)
        mutation_type = "Nonsense"
        description = f"Creates premature stop codon - truncates protein"
        impact = "VERY HIGH"
        
    elif is_stop_codon(original_codon):
        # Stop-loss mutation (removes stop codon)
        mutation_type = "Stop-loss"
        description = f"Removes stop codon - extends protein with {mutated_aa_name}"
        impact = "HIGH"
        
    else:
        # Missense mutation
        mutation_type = "Missense"
        description = f"Changes {original_aa_name} to {mutated_aa_name}"
        impact = "MODERATE"
    
    return {
        'type': mutation_type,
        'original_codon': original_codon,
        'mutated_codon': mutated_codon,
        'original_aa': original_aa,
        'mutated_aa': mutated_aa,
        'original_aa_name': original_aa_name,
        'mutated_aa_name': mutated_aa_name,
        'description': description,
        'impact': impact
    }


def calculate_pathogenicity_score(mutation_info, gene_name):
    """
    Calculate pathogenicity score (0-5) based on mutation characteristics.
    
    Args:
        mutation_info (dict): Mutation classification information
        gene_name (str): Name of the gene
        
    Returns:
        dict: Pathogenicity assessment including score, level, and recommendations
    """
    # Critical genes that warrant higher risk scores
    critical_genes = {
        'TP53', 'BRCA1', 'BRCA2', 'PTEN', 'APC', 'RB1', 'VHL',
        'MLH1', 'MSH2', 'MSH6', 'PMS2', 'KRAS', 'BRAF', 'EGFR'
    }
    
    # Base score from mutation type
    base_scores = {
        'Silent': 0,
        'Missense': 3,
        'Nonsense': 5,
        'Stop-loss': 4
    }
    
    mutation_type = mutation_info['type']
    score = base_scores.get(mutation_type, 3)
    
    # Adjust score for gene importance
    if gene_name.upper() in critical_genes:
        if mutation_type == 'Missense':
            score = 4  # Upgrade missense in critical genes
        # Nonsense and Stop-loss already at max or near-max
    
    # Determine risk level
    risk_levels = {
        0: "NO RISK",
        1: "VERY LOW",
        2: "LOW",
        3: "MODERATE",
        4: "HIGH",
        5: "VERY HIGH"
    }
    
    risk_level = risk_levels.get(score, "UNKNOWN")
    
    # Determine classification
    if score >= 4:
        classification = "LIKELY PATHOGENIC"
    elif score == 3:
        classification = "UNCERTAIN SIGNIFICANCE"
    elif score >= 1:
        classification = "LIKELY BENIGN"
    else:
        classification = "BENIGN"
    
    # Generate recommendations
    recommendations = []
    if score >= 4:
        recommendations = [
            "Genetic counseling strongly recommended",
            "Clinical confirmatory testing advised",
            "Family history evaluation warranted",
            "Enhanced screening protocol",
            "Consider genetic testing for family members"
        ]
    elif score == 3:
        recommendations = [
            "Consider genetic counseling",
            "Monitor for additional clinical evidence",
            "Family history may provide context"
        ]
    elif score >= 1:
        recommendations = [
            "Routine monitoring",
            "Low clinical concern"
        ]
    else:
        recommendations = [
            "No action required",
            "Silent mutation - no protein change"
        ]
    
    return {
        'score': score,
        'risk_level': risk_level,
        'classification': classification,
        'recommendations': recommendations,
        'is_critical_gene': gene_name.upper() in critical_genes
    }


def query_clinvar(gene_name, amino_acid_change):
    """
    Query ClinVar database for variant information.
    This is a simplified version - in production, use proper API.
    
    Args:
        gene_name (str): Gene symbol
        amino_acid_change (str): Amino acid change notation (e.g., R175H)
        
    Returns:
        dict: ClinVar information if available
    """
    try:
        # Note: This is a placeholder. Real implementation would use ClinVar API
        # https://www.ncbi.nlm.nih.gov/clinvar/
        
        print("Querying ClinVar database...")
        
        # Simulated response for demonstration
        known_variants = {
            'TP53': {
                'R175H': {
                    'significance': 'Pathogenic',
                    'submissions': 146,
                    'review_status': 'Reviewed by expert panel'
                }
            }
        }
        
        if gene_name in known_variants and amino_acid_change in known_variants[gene_name]:
            return {
                'found': True,
                'data': known_variants[gene_name][amino_acid_change]
            }
        
        return {'found': False, 'message': 'Variant not found in ClinVar'}
        
    except Exception as e:
        return {'found': False, 'error': str(e)}


def query_cosmic(gene_name):
    """
    Query COSMIC database for gene information.
    This is a simplified version.
    
    Args:
        gene_name (str): Gene symbol
        
    Returns:
        dict: COSMIC information if available
    """
    try:
        print("Querying COSMIC database...")
        
        # Simulated cancer gene information
        cancer_genes = {
            'TP53': {
                'role': 'Tumor suppressor',
                'cancers': ['Li-Fraumeni syndrome', 'Colorectal cancer', 'Breast cancer', 'Lung cancer'],
                'mutation_frequency': 'Very high (>50% of cancers)'
            },
            'KRAS': {
                'role': 'Oncogene',
                'cancers': ['Pancreatic cancer (90%)', 'Colorectal cancer (40%)', 'Lung cancer (30%)'],
                'mutation_frequency': 'High (25-30% of cancers)'
            },
            'BRAF': {
                'role': 'Oncogene',
                'cancers': ['Melanoma (50%)', 'Colorectal cancer (10%)', 'Thyroid cancer (40%)'],
                'mutation_frequency': 'Moderate'
            }
        }
        
        if gene_name in cancer_genes:
            return {
                'found': True,
                'data': cancer_genes[gene_name]
            }
        
        return {'found': False, 'message': 'Gene not found in COSMIC'}
        
    except Exception as e:
        return {'found': False, 'error': str(e)}


def generate_report(gene_name, position, original_nt, mutated_nt, 
                   mutation_info, pathogenicity, clinvar_data, cosmic_data):
    """
    Generate comprehensive text report.
    
    Args:
        gene_name (str): Gene name
        position (int): Mutation position
        original_nt (str): Original nucleotide
        mutated_nt (str): Mutated nucleotide
        mutation_info (dict): Mutation classification
        pathogenicity (dict): Pathogenicity assessment
        clinvar_data (dict): ClinVar query results
        cosmic_data (dict): COSMIC query results
        
    Returns:
        str: Formatted report text
    """
    report = []
    report.append("=" * 70)
    report.append("MUTATION ANALYSIS REPORT")
    report.append("=" * 70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # Basic Information
    report.append("BASIC INFORMATION")
    report.append("-" * 70)
    report.append(f"Gene: {gene_name}")
    report.append(f"Position: {position}")
    report.append(f"Original Nucleotide: {original_nt}")
    report.append(f"Mutated Nucleotide: {mutated_nt}")
    report.append("")
    
    # Mutation Classification
    report.append("MUTATION CLASSIFICATION")
    report.append("-" * 70)
    report.append(f"Type: {mutation_info['type']}")
    report.append(f"Original Codon: {mutation_info['original_codon']} ({mutation_info['original_aa_name']})")
    report.append(f"Mutated Codon: {mutation_info['mutated_codon']} ({mutation_info['mutated_aa_name']})")
    
    if mutation_info['type'] != 'Silent':
        codon_number = (position - 1) // 3 + 1
        aa_notation = f"{mutation_info['original_aa']}{codon_number}{mutation_info['mutated_aa']}"
        report.append(f"Amino Acid Change: {aa_notation}")
    
    report.append(f"Description: {mutation_info['description']}")
    report.append(f"Impact: {mutation_info['impact']}")
    report.append("")
    
    # Pathogenicity Assessment
    report.append("PATHOGENICITY ASSESSMENT")
    report.append("-" * 70)
    report.append(f"Risk Score: {pathogenicity['score']}/5")
    report.append(f"Risk Level: {pathogenicity['risk_level']}")
    report.append(f"Classification: {pathogenicity['classification']}")
    
    if pathogenicity['is_critical_gene']:
        report.append(f"Note: {gene_name} is a critical disease/cancer gene")
    
    if pathogenicity['score'] >= 4:
        report.append("\nWARNING: This mutation may significantly impact protein function")
    
    report.append("")
    
    # Clinical Recommendations
    report.append("CLINICAL RECOMMENDATIONS")
    report.append("-" * 70)
    for i, rec in enumerate(pathogenicity['recommendations'], 1):
        report.append(f"{i}. {rec}")
    report.append("")
    
    # Database Information
    report.append("DATABASE INFORMATION")
    report.append("-" * 70)
    
    # ClinVar
    if clinvar_data.get('found'):
        data = clinvar_data['data']
        report.append(f"ClinVar: {data['significance']}")
        report.append(f"  Submissions: {data['submissions']}")
        report.append(f"  Review Status: {data['review_status']}")
    else:
        report.append("ClinVar: Not found or not queried")
    
    report.append("")
    
    # COSMIC
    if cosmic_data.get('found'):
        data = cosmic_data['data']
        report.append(f"COSMIC Gene Information:")
        report.append(f"  Role: {data['role']}")
        report.append(f"  Mutation Frequency: {data['mutation_frequency']}")
        report.append(f"  Associated Cancers:")
        for cancer in data['cancers']:
            report.append(f"    - {cancer}")
    else:
        report.append("COSMIC: Not found or not queried")
    
    report.append("")
    report.append("=" * 70)
    report.append("DISCLAIMER")
    report.append("-" * 70)
    report.append("This analysis is for research and educational purposes only.")
    report.append("NOT intended for clinical diagnosis or medical decisions.")
    report.append("Always consult qualified genetic professionals for clinical interpretation.")
    report.append("=" * 70)
    
    return "\n".join(report)


def load_fasta(filename):
    """
    Load DNA sequence from FASTA file.
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        tuple: (sequence (str), description (str))
    """
    try:
        record = SeqIO.read(filename, "fasta")
        return str(record.seq), record.description
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {filename}")
    except Exception as e:
        raise Exception(f"Error reading FASTA file: {e}")


def analyze_mutation(sequence, gene_name, position, new_nucleotide):
    """
    Main analysis function.
    
    Args:
        sequence (str): DNA sequence
        gene_name (str): Gene name
        position (int): Position to mutate (1-indexed)
        new_nucleotide (str): New nucleotide (A/T/C/G)
        
    Returns:
        dict: Complete analysis results
    """
    # Clean and validate sequence
    sequence = sequence.replace(" ", "").replace("\n", "").replace("\r", "").upper()
    is_valid, error_msg = validate_sequence(sequence)
    
    if not is_valid:
        raise ValueError(f"Invalid sequence: {error_msg}")
    
    # Validate position
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} is out of range (1-{len(sequence)})")
    
    # Validate nucleotide
    new_nucleotide = new_nucleotide.upper()
    if new_nucleotide not in 'ATCG':
        raise ValueError(f"Invalid nucleotide: {new_nucleotide}. Must be A, T, C, or G")
    
    # Get original nucleotide
    original_nt = sequence[position - 1]
    
    # Check if mutation actually changes anything
    if original_nt == new_nucleotide:
        raise ValueError(f"Position {position} already contains nucleotide {new_nucleotide}")
    
    # Determine codon boundaries
    codon_start = ((position - 1) // 3) * 3
    codon_number = codon_start // 3 + 1
    
    # Extract original codon
    original_codon = sequence[codon_start:codon_start + 3]
    
    # Create mutated codon
    position_in_codon = (position - 1) % 3
    mutated_codon = list(original_codon)
    mutated_codon[position_in_codon] = new_nucleotide
    mutated_codon = ''.join(mutated_codon)
    
    # Classify mutation
    mutation_info = classify_mutation(original_codon, mutated_codon)
    
    # Calculate pathogenicity
    pathogenicity = calculate_pathogenicity_score(mutation_info, gene_name)
    
    # Query databases
    if mutation_info['type'] != 'Silent':
        aa_change = f"{mutation_info['original_aa']}{codon_number}{mutation_info['mutated_aa']}"
        clinvar_data = query_clinvar(gene_name, aa_change)
    else:
        clinvar_data = {'found': False, 'message': 'Silent mutations not typically in ClinVar'}
    
    cosmic_data = query_cosmic(gene_name)
    
    # Generate report
    report_text = generate_report(
        gene_name, position, original_nt, new_nucleotide,
        mutation_info, pathogenicity, clinvar_data, cosmic_data
    )
    
    return {
        'sequence': sequence,
        'gene_name': gene_name,
        'position': position,
        'original_nt': original_nt,
        'mutated_nt': new_nucleotide,
        'codon_number': codon_number,
        'mutation_info': mutation_info,
        'pathogenicity': pathogenicity,
        'clinvar_data': clinvar_data,
        'cosmic_data': cosmic_data,
        'report_text': report_text
    }


def main():
    """
    Main interactive function with user-friendly error handling.
    """
    print("=" * 70)
    print("MUTATION IMPACT ANALYZER")
    print("Fixed version with proper stop codon recognition")
    print("=" * 70)
    print()
    
    # Input method selection
    print("Choose input method:")
    print("1. Load from FASTA file")
    print("2. Paste sequence directly")
    print()
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    sequence = None
    gene_name = None
    
    if choice == '1':
        # Load from file
        filename = input("Enter FASTA filename (in sequences/ folder): ").strip()
        filepath = f"sequences/{filename}"
        
        try:
            sequence, description = load_fasta(filepath)
            print(f"\nLoaded: {description}")
            print(f"Sequence length: {len(sequence)} nucleotides")
            
            # Extract gene name from description if possible
            gene_name = input("Enter gene name: ").strip()
            
        except Exception as e:
            print(f"Error: {e}")
            return
    
    elif choice == '2':
        # Direct input with validation loop
        gene_name = input("Enter gene name (e.g., TP53, BRCA1): ").strip()
        
        while True:
            print("Enter DNA sequence (only A, T, C, G):")
            sequence = input().strip()
            
            # Clean sequence
            sequence = sequence.replace(" ", "").replace("\n", "").replace("\r", "").upper()
            
            # Validate sequence
            is_valid, error_msg = validate_sequence(sequence)
            
            if not is_valid:
                print(f"\nError: {error_msg}")
                
                # Check if it's an incomplete codon issue
                if "not divisible by 3" in error_msg:
                    seq_len = len(sequence)
                    remainder = seq_len % 3
                    needed = 3 - remainder if remainder != 0 else 0
                    
                    print(f"Your sequence has {seq_len} nucleotides.")
                    print(f"You need to add {needed} more nucleotide(s) to complete the last codon.")
                    print(f"Or remove {remainder} nucleotide(s) from the end.")
                
                retry = input("\nWould you like to re-enter the sequence? (y/n): ").strip().lower()
                if retry != 'y':
                    print("Exiting...")
                    return
                continue
            else:
                print(f"\nSequence validated successfully!")
                print(f"Length: {len(sequence)} nucleotides ({len(sequence)//3} codons)")
                break
    
    else:
        print("Invalid choice")
        return
    
    # Get mutation details with position validation loop
    print()
    
    while True:
        try:
            position = int(input(f"Enter position to mutate (1-{len(sequence)}): ").strip())
            
            # Check if position is in valid range
            if position < 1 or position > len(sequence):
                print(f"Error: Position must be between 1 and {len(sequence)}")
                continue
            
            # Show current nucleotide at that position
            current_nt = sequence[position - 1]
            codon_start = ((position - 1) // 3) * 3
            codon_number = codon_start // 3 + 1
            current_codon = sequence[codon_start:codon_start + 3]
            position_in_codon = (position - 1) % 3
            
            print(f"\nPosition {position} information:")
            print(f"  Current nucleotide: {current_nt}")
            print(f"  Current codon: {current_codon} (codon #{codon_number})")
            print(f"  Position in codon: {position_in_codon + 1} of 3")
            
            try:
                current_aa = translate_codon(current_codon)
                aa_name = AA_NAMES.get(current_aa, "Unknown")
                print(f"  Amino acid: {aa_name} ({current_aa})")
            except:
                print(f"  Amino acid: Unable to translate")
            
            break
            
        except ValueError:
            print("Error: Please enter a valid number")
            continue
    
    # Get new nucleotide with validation
    while True:
        new_nucleotide = input("\nEnter new nucleotide (A/T/C/G): ").strip().upper()
        
        if new_nucleotide not in 'ATCG':
            print(f"Error: Invalid nucleotide '{new_nucleotide}'. Must be A, T, C, or G")
            continue
        
        # Check if it's the same as current
        if new_nucleotide == current_nt:
            print(f"Error: Position {position} already contains nucleotide {new_nucleotide}")
            retry = input("Would you like to choose a different nucleotide? (y/n): ").strip().lower()
            if retry != 'y':
                print("Exiting...")
                return
            continue
        
        break
    
    # Show what the mutation will do
    mutated_codon_list = list(current_codon)
    mutated_codon_list[position_in_codon] = new_nucleotide
    mutated_codon = ''.join(mutated_codon_list)
    
    print(f"\nMutation preview:")
    print(f"  {current_codon} → {mutated_codon}")
    
    try:
        original_aa = translate_codon(current_codon)
        mutated_aa = translate_codon(mutated_codon)
        original_aa_name = AA_NAMES.get(original_aa, "Unknown")
        mutated_aa_name = AA_NAMES.get(mutated_aa, "Unknown")
        print(f"  {original_aa_name} ({original_aa}) → {mutated_aa_name} ({mutated_aa})")
    except:
        pass
    
    confirm = input("\nProceed with analysis? (y/n): ").strip().lower()
    if confirm != 'y':
        print("Analysis cancelled.")
        return
    
    # Perform analysis
    try:
        print("\nAnalyzing mutation...")
        results = analyze_mutation(sequence, gene_name, position, new_nucleotide)
        
        # Display report
        print("\n" + results['report_text'])
        
        # Save report to file
        output_filename = f"output/{gene_name}_P{position}{new_nucleotide}_report.txt"
        
        import os
        os.makedirs("output", exist_ok=True)
        
        with open(output_filename, 'w') as f:
            f.write(results['report_text'])
        
        print(f"\nReport saved to: {output_filename}")
        
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()