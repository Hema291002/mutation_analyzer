#!/usr/bin/env python3
"""
Enhanced Mutation Impact Analyzer - Standalone Version
Complete integration with disease search and pathway mapping
"""

from Bio import SeqIO
from Bio.Seq import Seq
import requests
import time
import os

# ============================================================
# CORE FUNCTIONS - Sequence Analysis
# ============================================================

def load_fasta(file_path):
    """Load FASTA file and return sequence information"""
    try:
        record = SeqIO.read(file_path, "fasta")
        return {
            'id': record.id,
            'description': record.description,
            'sequence': str(record.seq).upper()
        }
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return None
    except Exception as e:
        print(f"Error loading FASTA: {e}")
        return None


def validate_sequence(seq_string):
    """Validate that sequence contains only valid DNA nucleotides"""
    seq_string = seq_string.upper().replace(" ", "").replace("\n", "")
    valid_nucleotides = set('ATCGUN')
    invalid_chars = set(seq_string) - valid_nucleotides
    if invalid_chars:
        return False, f"Invalid nucleotides found: {invalid_chars}"
    return True, "Valid DNA sequence"


def parse_direct_input(seq_string):
    """Accept sequence pasted directly by user"""
    seq_string = seq_string.upper().replace(" ", "").replace("\n", "")
    is_valid, message = validate_sequence(seq_string)
    if is_valid:
        return seq_string, message
    else:
        return None, message


def apply_mutation(sequence, position, new_nucleotide):
    """Apply SNP mutation to sequence"""
    if position < 0 or position >= len(sequence):
        return None, f"Position out of range. Valid range: 0-{len(sequence)-1}"
    if new_nucleotide.upper() not in 'ATCG':
        return None, "Invalid nucleotide. Must be A, T, C, or G"
    
    new_nucleotide = new_nucleotide.upper()
    original_nt = sequence[position]
    if original_nt == new_nucleotide:
        return None, f"No mutation - position {position} already has {new_nucleotide}"
    
    seq_list = list(sequence)
    seq_list[position] = new_nucleotide
    mutated_seq = ''.join(seq_list)
    return mutated_seq, original_nt


def display_mutation(original_seq, mutated_seq, position, original_nt, new_nt):
    """Display the mutation in a readable format"""
    start = max(0, position - 10)
    end = min(len(original_seq), position + 11)
    
    print(f"\n{'='*70}")
    print(f"MUTATION APPLIED")
    print(f"{'='*70}")
    print(f"Position: {position} (0-based) / {position + 1} (1-based)")
    print(f"Change: {original_nt} → {new_nt}")
    print(f"\nContext (position {start}-{end}):")
    print(f"Original: {original_seq[start:end]}")
    print(f"          {' ' * (position - start)}^")
    print(f"Mutated:  {mutated_seq[start:end]}")
    print(f"{'='*70}\n")


def translate_sequence(dna_seq):
    """Translate DNA sequence to protein"""
    trim_length = len(dna_seq) - (len(dna_seq) % 3)
    seq_trimmed = dna_seq[:trim_length]
    seq = Seq(seq_trimmed)
    return str(seq.translate())


def classify_mutation(original_seq, mutated_seq, position):
    """Classify mutation type by comparing proteins"""
    original_protein = translate_sequence(original_seq)
    mutated_protein = translate_sequence(mutated_seq)
    
    codon_number = (position // 3) + 1
    codon_start = (codon_number - 1) * 3
    original_codon = original_seq[codon_start:codon_start + 3]
    mutated_codon = mutated_seq[codon_start:codon_start + 3]
    
    codon_index = codon_number - 1
    if codon_index < len(original_protein) and codon_index < len(mutated_protein):
        ref_aa = original_protein[codon_index]
        alt_aa = mutated_protein[codon_index]
        
        if ref_aa == alt_aa:
            mutation_type = "Silent/Synonymous"
            description = "No change in amino acid - mutation has no effect on protein"
        elif alt_aa == "*":
            mutation_type = "Nonsense"
            description = "Creates a STOP codon - protein is truncated"
        elif ref_aa == "*":
            mutation_type = "Stop-loss"
            description = "Removes STOP codon - protein is extended"
        else:
            mutation_type = "Missense"
            description = "Changes amino acid - may affect protein function"
        
        return {
            'type': mutation_type,
            'description': description,
            'codon_number': codon_number,
            'aa_position': codon_number,
            'original_codon': original_codon,
            'mutated_codon': mutated_codon,
            'ref_aa': ref_aa,
            'alt_aa': alt_aa,
            'notation': f"{ref_aa}{codon_number}{alt_aa}"
        }
    else:
        return {
            'type': 'Unknown',
            'description': 'Mutation near sequence boundary',
            'notation': 'N/A'
        }


def display_classification(classification):
    """Display mutation classification results"""
    print("\n" + "="*70)
    print("MUTATION CLASSIFICATION")
    print("="*70)
    
    print(f"\n🔬 Mutation Type: {classification['type']}")
    print(f"   {classification['description']}")
    
    print(f"\n📊 Details:")
    print(f"   Codon #{classification['codon_number']}")
    print(f"   Original codon: {classification['original_codon']}")
    print(f"   Mutated codon:  {classification['mutated_codon']}")
    print(f"   Original amino acid: {classification['ref_aa']}")
    print(f"   New amino acid:      {classification['alt_aa']}")
    print(f"   Standard notation: {classification['notation']}")
    
    if classification['type'] == 'Missense':
        print("⚠️  MISSENSE mutation detected!")
        print("    → Changes protein sequence")
        print("    → May affect function")
    elif classification['type'] == 'Nonsense':
        print("🛑 NONSENSE mutation detected!")
        print("    → Protein truncated, likely pathogenic")
    elif classification['type'] == 'Silent/Synonymous':
        print("✓ SILENT mutation detected")
        print("    → No change to protein sequence")
    
    print("="*70 + "\n")


# ============================================================
# ENHANCED FEATURES - Disease & Pathway Analysis
# ============================================================

def search_clinvar(gene_name, aa_change):
    """Search ClinVar for disease associations"""
    print(f"\n{'='*70}")
    print("🔍 SEARCHING CLINVAR DATABASE")
    print('='*70)
    print(f"Query: {gene_name} {aa_change}")
    
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_term = f"{gene_name}[gene] AND {aa_change}"
        search_url = f"{base_url}esearch.fcgi?db=clinvar&term={search_term}&retmode=json"
        
        response = requests.get(search_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            id_list = data.get('esearchresult', {}).get('idlist', [])
            
            if id_list:
                print(f"✓ Found {len(id_list)} ClinVar entries for this variant")
                print(f"  → This mutation has been clinically reported")
                return {'found': True, 'count': len(id_list), 'ids': id_list[:5]}
            else:
                print("ℹ️  No exact matches in ClinVar database")
                print("  → This specific variant may be novel")
                return {'found': False, 'count': 0}
        else:
            print(f"⚠️  ClinVar API temporarily unavailable (Status: {response.status_code})")
            return {'found': False, 'error': 'API error'}
            
    except Exception as e:
        print(f"⚠️  Could not connect to ClinVar: {e}")
        return {'found': False, 'error': str(e)}


def search_cosmic(gene_name):
    """Check COSMIC cancer database"""
    print(f"\n{'='*70}")
    print("🔍 CHECKING COSMIC CANCER DATABASE")
    print('='*70)
    
    cosmic_genes = {
        'TP53': {
            'cancer_types': ['Lung cancer', 'Colorectal cancer', 'Breast cancer', 'Ovarian cancer', 'Pancreatic cancer'],
            'mutation_frequency': 'Very High (>50% in many cancers)',
            'role': 'Tumor suppressor - Guardian of the genome',
            'hotspots': ['R175H', 'R248Q', 'R248W', 'R273H', 'R273C', 'R282W']
        },
        'BRCA1': {
            'cancer_types': ['Breast cancer', 'Ovarian cancer', 'Prostate cancer'],
            'mutation_frequency': 'High in hereditary cancers',
            'role': 'DNA repair - Homologous recombination',
            'hotspots': ['Throughout protein - truncating mutations common']
        },
        'BRCA2': {
            'cancer_types': ['Breast cancer', 'Ovarian cancer', 'Pancreatic cancer', 'Prostate cancer'],
            'mutation_frequency': 'High in hereditary cancers',
            'role': 'DNA repair - Homologous recombination',
            'hotspots': ['Throughout protein - truncating mutations common']
        },
        'EGFR': {
            'cancer_types': ['Lung cancer (NSCLC)', 'Glioblastoma', 'Colorectal cancer'],
            'mutation_frequency': 'High in lung adenocarcinoma (~15%)',
            'role': 'Growth factor receptor - Cell proliferation',
            'hotspots': ['L858R', 'T790M', 'Exon 19 deletions', 'G719X']
        },
        'KRAS': {
            'cancer_types': ['Pancreatic cancer', 'Colorectal cancer', 'Lung cancer'],
            'mutation_frequency': 'Very High (90% pancreatic, 40% colorectal)',
            'role': 'GTPase - Cell signaling',
            'hotspots': ['G12D', 'G12V', 'G12C', 'G13D', 'Q61H']
        },
        'PIK3CA': {
            'cancer_types': ['Breast cancer', 'Colorectal cancer', 'Endometrial cancer'],
            'mutation_frequency': 'High (~30% breast cancer)',
            'role': 'Kinase - PI3K/AKT pathway',
            'hotspots': ['E545K', 'H1047R', 'E542K']
        },
        'BRAF': {
            'cancer_types': ['Melanoma', 'Colorectal cancer', 'Thyroid cancer'],
            'mutation_frequency': 'Very High (50% melanoma)',
            'role': 'Kinase - MAPK pathway',
            'hotspots': ['V600E', 'V600K']
        }
    }
    
    if gene_name in cosmic_genes:
        info = cosmic_genes[gene_name]
        print(f"✓ {gene_name} is a well-characterized cancer gene")
        print(f"\n📊 Cancer Profile:")
        print(f"   Role: {info['role']}")
        print(f"   Mutation frequency: {info['mutation_frequency']}")
        print(f"\n🎗️  Associated cancers:")
        for cancer in info['cancer_types']:
            print(f"   • {cancer}")
        print(f"\n🔥 Known hotspot mutations:")
        print(f"   {', '.join(info['hotspots'][:5])}")
        return info
    else:
        print(f"ℹ️  {gene_name} not in common cancer gene database")
        print("  → May still be relevant - check literature")
        return None


def search_omim(gene_name):
    """Search OMIM for genetic diseases"""
    print(f"\n{'='*70}")
    print("🔍 SEARCHING OMIM GENETIC DISEASE DATABASE")
    print('='*70)
    
    omim_genes = {
        'TP53': {
            'diseases': ['Li-Fraumeni syndrome (LFS)', 'Various cancers'],
            'inheritance': 'Autosomal dominant',
            'omim_id': '191170',
            'description': 'Predisposition to multiple cancer types at young age'
        },
        'BRCA1': {
            'diseases': ['Hereditary breast and ovarian cancer syndrome (HBOC)'],
            'inheritance': 'Autosomal dominant',
            'omim_id': '113705',
            'description': 'High lifetime risk of breast and ovarian cancer'
        },
        'BRCA2': {
            'diseases': ['Hereditary breast and ovarian cancer syndrome (HBOC)', 'Fanconi anemia'],
            'inheritance': 'Autosomal dominant / recessive (FA)',
            'omim_id': '600185',
            'description': 'Cancer predisposition and DNA repair disorder'
        },
        'CFTR': {
            'diseases': ['Cystic fibrosis', 'Congenital bilateral absence of vas deferens'],
            'inheritance': 'Autosomal recessive',
            'omim_id': '602421',
            'description': 'Chloride channel dysfunction affecting lungs and pancreas'
        },
        'HBB': {
            'diseases': ['Sickle cell anemia', 'Beta-thalassemia'],
            'inheritance': 'Autosomal recessive',
            'omim_id': '141900',
            'description': 'Hemoglobin disorders causing anemia'
        },
        'DMD': {
            'diseases': ['Duchenne muscular dystrophy', 'Becker muscular dystrophy'],
            'inheritance': 'X-linked recessive',
            'omim_id': '300377',
            'description': 'Progressive muscle degeneration'
        }
    }
    
    if gene_name in omim_genes:
        info = omim_genes[gene_name]
        print(f"✓ {gene_name} associated with genetic disease(s)")
        print(f"\n🏥 Disease Information:")
        print(f"   OMIM ID: {info['omim_id']}")
        print(f"   Inheritance: {info['inheritance']}")
        print(f"\n💊 Associated diseases:")
        for disease in info['diseases']:
            print(f"   • {disease}")
        print(f"\n📝 Description:")
        print(f"   {info['description']}")
        return info
    else:
        print(f"ℹ️  No common Mendelian diseases found for {gene_name}")
        print("  → Gene may be involved in complex/multifactorial conditions")
        return None


def map_to_pathways(gene_name):
    """Map gene to biological pathways"""
    print(f"\n{'='*70}")
    print("🔍 MAPPING TO BIOLOGICAL PATHWAYS")
    print('='*70)
    
    pathway_db = {
        'TP53': {
            'pathways': [
                ('p53 signaling pathway', 'KEGG:04115'),
                ('Cell cycle checkpoint control', 'Reactome:R-HSA-69620'),
                ('Apoptosis signaling', 'Reactome:R-HSA-109581'),
                ('DNA damage response', 'Reactome:R-HSA-5693532'),
                ('Cellular senescence', 'KEGG:04218')
            ],
            'summary': 'Master regulator of cell fate decisions'
        },
        'BRCA1': {
            'pathways': [
                ('Homologous recombination repair', 'KEGG:03440'),
                ('BRCA1-mediated DNA repair', 'Reactome:R-HSA-5685942'),
                ('Cell cycle checkpoint', 'Reactome:R-HSA-69620'),
                ('Fanconi anemia pathway', 'KEGG:03460')
            ],
            'summary': 'Critical DNA repair and genome stability'
        },
        'BRCA2': {
            'pathways': [
                ('Homologous recombination', 'KEGG:03440'),
                ('DNA double-strand break repair', 'Reactome:R-HSA-5693538'),
                ('Fanconi anemia pathway', 'KEGG:03460')
            ],
            'summary': 'Essential for DNA repair by homologous recombination'
        },
        'EGFR': {
            'pathways': [
                ('EGFR tyrosine kinase inhibitor resistance', 'KEGG:01521'),
                ('MAPK/ERK signaling pathway', 'KEGG:04010'),
                ('PI3K-Akt signaling', 'KEGG:04151'),
                ('Cell proliferation', 'Reactome:R-HSA-69278')
            ],
            'summary': 'Growth factor signaling and cell proliferation'
        },
        'KRAS': {
            'pathways': [
                ('RAS signaling pathway', 'KEGG:04014'),
                ('MAPK/ERK pathway', 'KEGG:04010'),
                ('PI3K-Akt signaling', 'KEGG:04151'),
                ('Regulation of actin cytoskeleton', 'KEGG:04810')
            ],
            'summary': 'Oncogenic signaling hub'
        }
    }
    
    if gene_name in pathway_db:
        data = pathway_db[gene_name]
        print(f"✓ {gene_name}: {data['summary']}")
        print(f"\n🧬 Biological pathways ({len(data['pathways'])} found):")
        for i, (pathway, pathway_id) in enumerate(data['pathways'], 1):
            print(f"   {i}. {pathway}")
            print(f"      ID: {pathway_id}")
        return data['pathways']
    else:
        print(f"ℹ️  Limited pathway data for {gene_name}")
        return []


def predict_pathogenicity(classification, gene_name):
    """Predict pathogenicity with detailed scoring"""
    print(f"\n{'='*70}")
    print("🔍 PATHOGENICITY PREDICTION")
    print('='*70)
    
    score = 0
    factors = []
    
    # Factor 1: Mutation type (0-3 points)
    if classification['type'] == 'Nonsense':
        score += 3
        factors.append("⚠️  Nonsense mutation → Truncated protein (HIGH IMPACT)")
    elif classification['type'] == 'Missense':
        score += 2
        factors.append("⚠️  Missense mutation → Altered amino acid (MODERATE IMPACT)")
    elif classification['type'] == 'Stop-loss':
        score += 3
        factors.append("⚠️  Stop-loss → Extended protein (HIGH IMPACT)")
    elif classification['type'] == 'Silent/Synonymous':
        score += 0
        factors.append("✓ Silent mutation → No protein change (LOW IMPACT)")
    
    # Factor 2: Gene importance (0-2 points)
    critical_genes = {
        'TP53': 2, 'BRCA1': 2, 'BRCA2': 2, 'EGFR': 2, 'KRAS': 2,
        'PIK3CA': 2, 'BRAF': 2, 'CFTR': 2, 'DMD': 2
    }
    
    if gene_name in critical_genes:
        gene_score = critical_genes[gene_name]
        score += gene_score
        factors.append(f"⚠️  {gene_name} is a critical disease/cancer gene")
    else:
        factors.append(f"ℹ️  {gene_name} - gene significance unknown")
    
    # Determine classification
    if score >= 4:
        prediction = "LIKELY PATHOGENIC"
        risk = "HIGH"
        color = "🔴"
        recommendation = "Genetic counseling and clinical testing STRONGLY recommended"
    elif score >= 2:
        prediction = "UNCERTAIN SIGNIFICANCE (VUS)"
        risk = "MODERATE"
        color = "🟡"
        recommendation = "Further investigation recommended - consult specialist"
    else:
        prediction = "LIKELY BENIGN"
        risk = "LOW"
        color = "🟢"
        recommendation = "Routine monitoring - low clinical concern"
    
    print(f"\n{color} PREDICTION: {prediction}")
    print(f"{'='*70}")
    print(f"Risk Level: {risk}")
    print(f"Pathogenicity Score: {score}/5")
    print(f"\n📋 Contributing Factors:")
    for factor in factors:
        print(f"   {factor}")
    
    print(f"\n💡 CLINICAL RECOMMENDATION:")
    print(f"   {recommendation}")
    
    if risk == "HIGH":
        print(f"\n⚠️  HIGH RISK ACTIONS:")
        print(f"   1. Discuss with genetic counselor")
        print(f"   2. Consider confirmatory testing")
        print(f"   3. Family history evaluation")
        print(f"   4. Cancer screening protocol")
    
    return {
        'prediction': prediction,
        'risk': risk,
        'score': score,
        'factors': factors,
        'recommendation': recommendation
    }


def generate_final_report(gene_name, classification, clinvar, cosmic, omim, pathways, pathogenicity):
    """Generate comprehensive final report"""
    print(f"\n{'='*70}")
    print(" " * 15 + "COMPREHENSIVE ANALYSIS REPORT")
    print('='*70)
    
    print(f"\n📋 MUTATION SUMMARY:")
    print(f"   Gene: {gene_name}")
    print(f"   Variant: {classification['notation']}")
    print(f"   Type: {classification['type']}")
    print(f"   DNA Change: {classification['original_codon']} → {classification['mutated_codon']}")
    print(f"   Protein Change: {classification['ref_aa']} → {classification['alt_aa']}")
    
    print(f"\n🔬 MOLECULAR IMPACT:")
    print(f"   Codon Position: {classification['codon_number']}")
    print(f"   Amino Acid Position: {classification['aa_position']}")
    print(f"   Effect: {classification['description']}")
    
    print(f"\n⚕️  CLINICAL SIGNIFICANCE:")
    print(f"   Classification: {pathogenicity['prediction']}")
    print(f"   Risk Level: {pathogenicity['risk']}")
    print(f"   Confidence Score: {pathogenicity['score']}/5")
    
    if clinvar and clinvar.get('found'):
        print(f"\n📚 CLINVAR STATUS:")
        print(f"   ✓ {clinvar['count']} clinical report(s) found")
        print(f"   → Previously observed in patients")
    
    if cosmic:
        print(f"\n🎗️  CANCER ASSOCIATIONS:")
        for cancer in cosmic.get('cancer_types', [])[:3]:
            print(f"   • {cancer}")
    
    if omim:
        print(f"\n🏥 GENETIC DISEASE:")
        for disease in omim.get('diseases', []):
            print(f"   • {disease}")
    
    if pathways:
        print(f"\n🧬 AFFECTED PATHWAYS:")
        for pathway, pathway_id in pathways[:3]:
            print(f"   • {pathway}")
    
    print(f"\n💊 RECOMMENDATIONS:")
    print(f"   {pathogenicity['recommendation']}")
    
    print(f"\n📅 Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print('='*70)
    
    # Save to file
    try:
        os.makedirs('output', exist_ok=True)
        filename = f"output/{gene_name}_{classification['notation']}_report.txt"
        
        with open(filename, 'w') as f:
            f.write(f"{'='*70}\n")
            f.write(f"MUTATION ANALYSIS REPORT\n")
            f.write(f"{'='*70}\n\n")
            f.write(f"Gene: {gene_name}\n")
            f.write(f"Variant: {classification['notation']}\n")
            f.write(f"Classification: {pathogenicity['prediction']}\n")
            f.write(f"Risk: {pathogenicity['risk']}\n")
            f.write(f"Score: {pathogenicity['score']}/5\n\n")
            f.write(f"Analysis performed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        print(f"\n✓ Report saved: {filename}")
    except Exception as e:
        print(f"\n⚠️  Could not save report: {e}")


# ============================================================
# MAIN PROGRAM
# ============================================================

def main():
    print("="*70)
    print(" " * 10 + "ENHANCED MUTATION IMPACT ANALYZER")
    print(" " * 8 + "with Disease Association & Pathway Mapping")
    print("="*70)
    
    # Step 1: Load sequence
    choice = input("\nChoose input method:\n1. Load FASTA file\n2. Paste sequence directly\nEnter choice (1/2): ")
    
    if choice == '1':
        file_path = input("Enter FASTA file path: ").strip()
        fasta_data = load_fasta(file_path)
        if not fasta_data:
            print("Failed to load file. Exiting.")
            return
        sequence = fasta_data['sequence']
        gene_name = input("Enter gene name (e.g., TP53): ").strip().upper()
        print(f"\n✓ Loaded: {fasta_data['id']} ({len(sequence)} bp)")
    else:
        seq_input = input("\nPaste your DNA sequence: ").strip()
        sequence, msg = parse_direct_input(seq_input)
        if not sequence:
            print(f"✗ {msg}. Exiting.")
            return
        gene_name = input("Enter gene name: ").strip().upper()
        print(f"✓ {msg} ({len(sequence)} bp)")
    
    print(f"Total codons: {len(sequence)//3}")
    
    # Step 2: Get mutation
    while True:
        try:
            position = int(input(f"\nEnter nucleotide position to mutate (1-{len(sequence)}): ")) - 1
            if 0 <= position < len(sequence):
                break
            print(f"✗ Position out of range.")
        except ValueError:
            print("✗ Please enter a valid number")
    
    current_nt = sequence[position]
    while True:
        new_nt = input(f"Enter new nucleotide to replace {current_nt} (A/T/C/G): ").strip().upper()
        if new_nt in 'ATCG':
            if new_nt != current_nt:
                break
            print(f"⚠️  Already {new_nt}. Choose different.")
        else:
            print("✗ Invalid nucleotide")
    
    # Step 3: Apply & classify
    mutated_seq, _ = apply_mutation(sequence, position, new_nt)
    if not mutated_seq:
        print("Mutation failed")
        return
    
    display_mutation(sequence, mutated_seq, position, current_nt, new_nt)
    classification = classify_mutation(sequence, mutated_seq, position)
    display_classification(classification)
    
    # Step 4: Enhanced analysis
    print("\n" + "="*70)
    print("PERFORMING ADVANCED ANALYSIS...")
    print("="*70)
    print("⏳ Searching medical databases and mapping pathways...")
    print("   This may take 10-15 seconds...\n")
    
    clinvar_results = search_clinvar(gene_name, classification['notation'])
    time.sleep(1)  # Rate limiting
    
    cosmic_results = search_cosmic(gene_name)
    omim_results = search_omim(gene_name)
    pathways = map_to_pathways(gene_name)
    pathogenicity = predict_pathogenicity(classification, gene_name)
    
    # Step 5: Final report
    generate_final_report(gene_name, classification, clinvar_results, 
                         cosmic_results, omim_results, pathways, pathogenicity)
    
    # Step 6: Generate visualizations
    viz_choice = input("\n🎨 Generate interactive visualizations? (y/n): ").strip().lower()
    if viz_choice == 'y':
        try:
            from visualize_results import create_comprehensive_dashboard
            print("\n⏳ Creating interactive dashboard...")
            dashboard_file = create_comprehensive_dashboard(
                gene_name, sequence, mutated_seq, position,
                classification, pathogenicity, cosmic_results, omim_results, pathways
            )
            print(f"\n✨ SUCCESS! Open {dashboard_file} in your browser!")
        except ImportError:
            print("\n⚠️  Visualization module not found. Make sure visualize_results.py is in the same directory.")
        except Exception as e:
            print(f"\n⚠️  Could not generate visualizations: {e}")
    
    print("\n✅ ANALYSIS COMPLETE!")
    print("   Check the 'output' folder for all reports and visualizations")


if __name__ == "__main__":
    main()