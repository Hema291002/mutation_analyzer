"""
Mutation Impact Analyzer - Main Module
Analyzes DNA mutations and their impact on protein function
Generates comprehensive reports and visualizations
"""

import re
from pathogenicity import PathogenicityPredictor


# Genetic Code Dictionary
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
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

# Amino Acid Properties
AA_PROPERTIES = {
    'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'],
    'polar': ['S', 'T', 'C', 'Y', 'N', 'Q'],
    'charged_positive': ['K', 'R', 'H'],
    'charged_negative': ['D', 'E'],
    'special': ['G'],
    'stop': ['*']
}

# COSMIC Database (simplified)
COSMIC_DB = {
    'TP53': {
        'R175H': {'count': 450, 'cancer_types': ['Breast', 'Lung', 'Colorectal'], 'role': 'DNA binding'},
        'R248Q': {'count': 380, 'cancer_types': ['Breast', 'Ovarian'], 'role': 'DNA binding'},
        'R273H': {'count': 420, 'cancer_types': ['Colorectal', 'Lung'], 'role': 'DNA binding'},
    },
    'BRCA1': {
        'C61G': {'count': 120, 'cancer_types': ['Breast', 'Ovarian'], 'role': 'RING domain'},
        'M1652I': {'count': 85, 'cancer_types': ['Breast'], 'role': 'BRCT domain'},
    },
    'EGFR': {
        'L858R': {'count': 890, 'cancer_types': ['Lung'], 'role': 'Kinase domain'},
        'T790M': {'count': 650, 'cancer_types': ['Lung'], 'role': 'Kinase domain'},
    }
}


def validate_sequence(sequence):
    """
    Validate DNA sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequence:
        return False, "Empty sequence"
    
    if not re.match('^[ATCG]+$', sequence.upper()):
        invalid_chars = set(sequence.upper()) - set('ATCG')
        return False, f"Invalid characters found: {', '.join(invalid_chars)}"
    
    if len(sequence) % 3 != 0:
        return False, f"Sequence length ({len(sequence)}) is not divisible by 3"
    
    return True, "Valid"


def load_fasta(filepath):
    """
    Load sequence from FASTA file
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        tuple: (sequence, description)
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            raise ValueError("Empty file")
        
        description = lines[0].strip().lstrip('>')
        sequence = ''.join(line.strip() for line in lines[1:] if line.strip())
        sequence = sequence.upper()
        
        is_valid, error = validate_sequence(sequence)
        if not is_valid:
            raise ValueError(f"Invalid sequence in FASTA file: {error}")
        
        return sequence, description
    
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    except Exception as e:
        raise Exception(f"Error reading FASTA file: {e}")


def translate_codon(codon):
    """Translate DNA codon to amino acid"""
    return GENETIC_CODE.get(codon.upper(), 'X')


def get_aa_property(aa):
    """Get amino acid property category"""
    for prop, amino_acids in AA_PROPERTIES.items():
        if aa in amino_acids:
            return prop
    return 'unknown'


def classify_mutation(original_aa, mutated_aa, original_codon, mutated_codon):
    """
    Classify mutation type and impact
    
    Returns:
        dict: Classification information
    """
    if original_aa == mutated_aa:
        return {
            'type': 'Synonymous',
            'impact': 'Silent mutation - no amino acid change',
            'severity': 'benign'
        }
    
    if mutated_aa == '*':
        return {
            'type': 'Nonsense',
            'impact': 'Premature stop codon - protein truncation',
            'severity': 'high'
        }
    
    if original_aa == '*':
        return {
            'type': 'Nonstop',
            'impact': 'Stop codon lost - extended protein',
            'severity': 'high'
        }
    
    # Missense mutation analysis
    orig_prop = get_aa_property(original_aa)
    mut_prop = get_aa_property(mutated_aa)
    
    if orig_prop == mut_prop:
        return {
            'type': 'Missense (Conservative)',
            'impact': f'Similar amino acid properties ({orig_prop})',
            'severity': 'moderate'
        }
    else:
        return {
            'type': 'Missense (Non-conservative)',
            'impact': f'Changed from {orig_prop} to {mut_prop}',
            'severity': 'high'
        }


def check_cosmic_database(gene_name, mutation_notation):
    """
    Check if mutation is in COSMIC database
    
    Args:
        gene_name: Gene name
        mutation_notation: e.g., 'R175H'
        
    Returns:
        dict: COSMIC data or empty dict
    """
    gene_data = COSMIC_DB.get(gene_name.upper(), {})
    mutation_data = gene_data.get(mutation_notation, None)
    
    if mutation_data:
        return {
            'found': True,
            'data': mutation_data
        }
    
    return {'found': False}


def analyze_mutation(sequence, gene_name, position, new_nucleotide):
    """
    Main analysis function
    
    Args:
        sequence: Original DNA sequence
        gene_name: Name of the gene
        position: 1-indexed position to mutate
        new_nucleotide: New nucleotide (A/T/C/G)
        
    Returns:
        dict: Complete analysis results
    """
    # Validate inputs
    is_valid, error = validate_sequence(sequence)
    if not is_valid:
        raise ValueError(f"Invalid sequence: {error}")
    
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} out of range (1-{len(sequence)})")
    
    if new_nucleotide.upper() not in 'ATCG':
        raise ValueError(f"Invalid nucleotide: {new_nucleotide}")
    
    # Convert to 0-indexed
    pos_idx = position - 1
    
    # Get original nucleotide
    original_nucleotide = sequence[pos_idx]
    
    if original_nucleotide == new_nucleotide.upper():
        raise ValueError(f"Position {position} already contains {new_nucleotide}")
    
    # Determine codon position
    codon_number = (pos_idx // 3) + 1
    codon_start = (pos_idx // 3) * 3
    position_in_codon = (pos_idx % 3) + 1
    
    # Get original codon
    original_codon = sequence[codon_start:codon_start + 3]
    original_aa = translate_codon(original_codon)
    
    # Create mutated codon
    mutated_codon = list(original_codon)
    mutated_codon[position_in_codon - 1] = new_nucleotide.upper()
    mutated_codon = ''.join(mutated_codon)
    mutated_aa = translate_codon(mutated_codon)
    
    # Classify mutation
    classification = classify_mutation(original_aa, mutated_aa, original_codon, mutated_codon)
    
    # Create mutation notation
    mutation_notation = f"{original_aa}{codon_number}{mutated_aa}"
    
    # Check COSMIC database
    cosmic_data = check_cosmic_database(gene_name, mutation_notation)
    
    # Pathogenicity prediction
    predictor = PathogenicityPredictor()
    pathogenicity = predictor.predict(
        original_aa=original_aa,
        mutated_aa=mutated_aa,
        position=codon_number,
        gene_name=gene_name,
        cosmic_found=cosmic_data['found']
    )
    
    # Generate report
    report = generate_report(
        gene_name=gene_name,
        sequence_length=len(sequence),
        position=position,
        codon_number=codon_number,
        position_in_codon=position_in_codon,
        original_nucleotide=original_nucleotide,
        new_nucleotide=new_nucleotide,
        original_codon=original_codon,
        mutated_codon=mutated_codon,
        original_aa=original_aa,
        mutated_aa=mutated_aa,
        mutation_notation=mutation_notation,
        classification=classification,
        cosmic_data=cosmic_data,
        pathogenicity=pathogenicity
    )
    
    return {
        'gene_name': gene_name,
        'sequence_length': len(sequence),
        'position': position,
        'codon_number': codon_number,
        'mutation_info': {
            'original_nucleotide': original_nucleotide,
            'new_nucleotide': new_nucleotide,
            'original_codon': original_codon,
            'mutated_codon': mutated_codon,
            'original_aa': original_aa,
            'mutated_aa': mutated_aa,
            'notation': mutation_notation,
            'type': classification['type']
        },
        'classification': classification,
        'cosmic_data': cosmic_data,
        'pathogenicity': pathogenicity,
        'report_text': report
    }


def generate_report(gene_name, sequence_length, position, codon_number, position_in_codon,
                   original_nucleotide, new_nucleotide, original_codon, mutated_codon,
                   original_aa, mutated_aa, mutation_notation, classification,
                   cosmic_data, pathogenicity):
    """Generate formatted text report"""
    
    report = []
    report.append("=" * 80)
    report.append(f"MUTATION ANALYSIS REPORT: {gene_name}")
    report.append("=" * 80)
    report.append("")
    
    # Basic Information
    report.append("SEQUENCE INFORMATION:")
    report.append(f"  Gene:               {gene_name}")
    report.append(f"  Sequence Length:    {sequence_length} nucleotides ({sequence_length//3} codons)")
    report.append(f"  Mutation Position:  {position}")
    report.append(f"  Codon Number:       {codon_number}")
    report.append(f"  Position in Codon:  {position_in_codon}")
    report.append("")
    
    # Mutation Details
    report.append("MUTATION DETAILS:")
    report.append(f"  Original Nucleotide: {original_nucleotide}")
    report.append(f"  Mutated Nucleotide:  {new_nucleotide}")
    report.append(f"  Original Codon:      {original_codon}")
    report.append(f"  Mutated Codon:       {mutated_codon}")
    report.append(f"  Original Amino Acid: {original_aa}")
    report.append(f"  Mutated Amino Acid:  {mutated_aa}")
    report.append(f"  Notation:            {mutation_notation}")
    report.append("")
    
    # Classification
    report.append("MUTATION CLASSIFICATION:")
    report.append(f"  Type:     {classification['type']}")
    report.append(f"  Impact:   {classification['impact']}")
    report.append(f"  Severity: {classification['severity'].upper()}")
    report.append("")
    
    # COSMIC Database
    report.append("COSMIC DATABASE CHECK:")
    if cosmic_data['found']:
        data = cosmic_data['data']
        report.append(f"  ✓ FOUND in COSMIC database")
        report.append(f"  Occurrences:   {data['count']}")
        report.append(f"  Cancer Types:  {', '.join(data['cancer_types'])}")
        report.append(f"  Protein Role:  {data['role']}")
    else:
        report.append(f"  ✗ Not found in COSMIC database")
    report.append("")
    
    # Pathogenicity
    report.append("PATHOGENICITY PREDICTION:")
    report.append(f"  Score:          {pathogenicity['score']:.2f}/10")
    report.append(f"  Risk Level:     {pathogenicity['risk_level']}")
    report.append(f"  Classification: {pathogenicity['classification']}")
    report.append("")
    report.append("  Factors Contributing to Score:")
    for factor in pathogenicity['factors']:
        report.append(f"    • {factor}")
    report.append("")
    
    # Clinical Recommendations
    report.append("CLINICAL RECOMMENDATIONS:")
    for i, rec in enumerate(pathogenicity['recommendations'], 1):
        report.append(f"  {i}. {rec}")
    report.append("")
    
    report.append("=" * 80)
    report.append("END OF REPORT")
    report.append("=" * 80)
    
    return '\n'.join(report)


def main():
    """Main function - runs the complete analysis with HTML generation"""
    from visualize_results import create_comprehensive_dashboard
    import os
    
    print("=" * 70)
    print("MUTATION IMPACT ANALYZER")
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
        filename = input("Enter FASTA filename (e.g., TP53.fasta): ").strip()
        
        # Clean up the filename - remove any leading path separators
        filename = filename.replace("sequences/", "").strip()
        filepath = f"sequences/{filename}"
        
        print(f"\nLooking for file: {filepath}")
        
        try:
            sequence, description = load_fasta(filepath)
            print(f"\n✅ Loaded: {description}")
            print(f"Sequence length: {len(sequence)} nucleotides")
            gene_name = input("Enter gene name (e.g., TP53): ").strip()
        except FileNotFoundError:
            print(f"\n❌ Error: File not found at {filepath}")
            print(f"\nAvailable files in sequences/ folder:")
            try:
                files = os.listdir("sequences/")
                fasta_files = [f for f in files if f.endswith(('.fasta', '.fa', '.fna'))]
                if fasta_files:
                    for f in fasta_files:
                        print(f"   - {f}")
                else:
                    print("   (No FASTA files found)")
            except:
                print("   (Could not list files)")
            return
        except Exception as e:
            print(f"\n❌ Error: {e}")
            return
    
    elif choice == '2':
        gene_name = input("Enter gene name (e.g., TP53, BRCA1): ").strip()
        
        while True:
            print("\nEnter DNA sequence (only A, T, C, G):")
            sequence = input().strip()
            sequence = sequence.replace(" ", "").replace("\n", "").replace("\r", "").upper()
            
            is_valid, error_msg = validate_sequence(sequence)
            
            if not is_valid:
                print(f"\n❌ Error: {error_msg}")
                retry = input("\nWould you like to re-enter the sequence? (y/n): ").strip().lower()
                if retry != 'y':
                    return
                continue
            else:
                print(f"\n✅ Sequence validated!")
                print(f"Length: {len(sequence)} nucleotides ({len(sequence)//3} codons)")
                break
    else:
        print("❌ Invalid choice")
        return
    
    # Get mutation position
    print()
    print("-" * 70)
    while True:
        try:
            position = int(input(f"Enter position to mutate (1-{len(sequence)}): ").strip())
            
            if position < 1 or position > len(sequence):
                print(f"❌ Error: Position must be between 1 and {len(sequence)}")
                continue
            
            current_nt = sequence[position - 1]
            codon_start = ((position - 1) // 3) * 3
            current_codon = sequence[codon_start:codon_start + 3]
            
            print(f"\n📍 Current nucleotide at position {position}: {current_nt}")
            print(f"📍 Current codon: {current_codon}")
            break
        except ValueError:
            print("❌ Error: Please enter a valid number")
            continue
    
    # Get new nucleotide
    current_nt = sequence[position - 1]
    while True:
        new_nucleotide = input("\nEnter new nucleotide (A/T/C/G): ").strip().upper()
        
        if new_nucleotide not in 'ATCG':
            print(f"❌ Error: Invalid nucleotide. Must be A, T, C, or G")
            continue
        
        if new_nucleotide == current_nt:
            print(f"❌ Error: Position already contains {new_nucleotide}")
            retry = input("Choose a different nucleotide? (y/n): ").strip().lower()
            if retry != 'y':
                return
            continue
        
        break
    
    # Perform analysis
    try:
        print("\n" + "="*70)
        print("🔬 ANALYZING MUTATION...")
        print("="*70)
        
        results = analyze_mutation(sequence, gene_name, position, new_nucleotide)
        
        # Display text report
        print("\n" + results['report_text'])
        
        # Create output directory
        os.makedirs("output", exist_ok=True)
        
        # Generate filename
        mutation_label = f"{results['mutation_info']['original_aa']}{results['codon_number']}{results['mutation_info']['mutated_aa']}"
        
        # Save text report
        text_filename = f"output/{gene_name}_{mutation_label}_report.txt"
        with open(text_filename, 'w') as f:
            f.write(results['report_text'])
        print(f"\n📄 Text report saved: {text_filename}")
        
        # Generate HTML dashboard
        print("\n" + "="*70)
        print("🎨 GENERATING INTERACTIVE HTML DASHBOARD...")
        print("="*70)
        
        # Create mutated sequence
        mutated_sequence = list(sequence)
        mutated_sequence[position - 1] = new_nucleotide
        mutated_sequence = ''.join(mutated_sequence)
        
        # Prepare visualization data
        viz_classification = {
            'notation': mutation_label,
            'type': results['mutation_info']['type'],
            'original_codon': results['mutation_info']['original_codon'],
            'mutated_codon': results['mutation_info']['mutated_codon'],
            'ref_aa': results['mutation_info']['original_aa'],
            'alt_aa': results['mutation_info']['mutated_aa'],
            'codon_number': results['codon_number']
        }
        
        viz_pathogenicity = {
            'score': results['pathogenicity']['score'],
            'risk': results['pathogenicity']['risk_level'],
            'prediction': results['pathogenicity']['classification'],
            'recommendation': results['pathogenicity']['recommendations'][0] if results['pathogenicity']['recommendations'] else 'See full report'
        }
        
        # Extract pathway data
        pathways = []
        if results['cosmic_data'].get('found'):
            cosmic_data = results['cosmic_data'].get('data', {})
            if 'role' in cosmic_data:
                pathways.append((cosmic_data['role'], 'COSMIC'))
        
        # Create dashboard
        try:
            html_file = create_comprehensive_dashboard(
                gene_name=gene_name,
                sequence=sequence,
                mutated_seq=mutated_sequence,
                position=position - 1,
                classification=viz_classification,
                pathogenicity=viz_pathogenicity,
                cosmic=results['cosmic_data'],
                omim={},
                pathways=pathways
            )
            
            print("\n" + "="*70)
            print("✅ ANALYSIS COMPLETE!")
            print("="*70)
            print(f"\n📊 Files generated:")
            print(f"   1. 📄 Text Report:    {text_filename}")
            print(f"   2. 🌐 HTML Dashboard: {html_file}")
            print(f"\n💡 Open HTML in browser: google-chrome {html_file}")
            print("="*70)
            
        except Exception as viz_error:
            print(f"\n⚠️  Warning: Could not create HTML dashboard: {viz_error}")
            print(f"   Text report saved: {text_filename}")
            
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️  Analysis cancelled by user.")
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
