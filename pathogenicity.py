"""
Pathogenicity Assessment Module
Evaluates the clinical significance of mutations
"""

# Amino acid properties
AMINO_ACID_PROPERTIES = {
    'nonpolar': ['G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W'],
    'polar': ['S', 'T', 'C', 'Y', 'N', 'Q'],
    'acidic': ['D', 'E'],
    'basic': ['K', 'R', 'H']
}

# Mutation type severity
MUTATION_SEVERITY = {
    'Silent/Synonymous': {
        'base_score': 0,
        'severity': 'Benign',
        'color': '🟢',
        'clinical_action': 'No action needed'
    },
    'Missense': {
        'base_score': 2,
        'severity': 'Variable',
        'color': '🟡',
        'clinical_action': 'Further evaluation recommended'
    },
    'Nonsense': {
        'base_score': 4,
        'severity': 'Pathogenic',
        'color': '🔴',
        'clinical_action': 'Genetic counseling strongly recommended'
    },
    'Stop-loss': {
        'base_score': 3,
        'severity': 'Likely Pathogenic',
        'color': '🟠',
        'clinical_action': 'Genetic counseling recommended'
    }
}


def get_aa_property(amino_acid):
    """Get the property category of an amino acid"""
    for prop, aas in AMINO_ACID_PROPERTIES.items():
        if amino_acid in aas:
            return prop
    return 'unknown'


def assess_biochemical_change(ref_aa, alt_aa):
    """
    Assess the biochemical impact of amino acid substitution
    
    Returns:
        dict: Impact assessment with score and description
    """
    ref_prop = get_aa_property(ref_aa)
    alt_prop = get_aa_property(alt_aa)
    
    # Same property = conservative change
    if ref_prop == alt_prop:
        return {
            'conservative': True,
            'impact_score': 1,
            'description': f'Conservative change within {ref_prop} group',
            'severity_modifier': 0
        }
    else:
        return {
            'conservative': False,
            'impact_score': 3,
            'description': f'Non-conservative change: {ref_prop} → {alt_prop}',
            'severity_modifier': 1
        }


def calculate_pathogenicity_score(mutation_type, ref_aa, alt_aa, is_hotspot=False):
    """
    Calculate comprehensive pathogenicity score
    
    Args:
        mutation_type (str): Type of mutation
        ref_aa (str): Reference amino acid
        alt_aa (str): Alternate amino acid
        is_hotspot (bool): Whether mutation is in known hotspot
    
    Returns:
        dict: Pathogenicity assessment
    """
    # Get base severity
    severity_info = MUTATION_SEVERITY.get(mutation_type, {
        'base_score': 1,
        'severity': 'Unknown',
        'color': '⚪',
        'clinical_action': 'Uncertain significance'
    })
    
    score = severity_info['base_score']
    
    # Enhanced assessment for missense mutations
    if mutation_type == 'Missense':
        biochem = assess_biochemical_change(ref_aa, alt_aa)
        score += biochem['severity_modifier']
        
        # Hotspot bonus
        if is_hotspot:
            score += 1
            severity_info['severity'] = 'Likely Pathogenic'
        elif biochem['conservative']:
            severity_info['severity'] = 'Likely Benign'
        else:
            severity_info['severity'] = 'Uncertain Significance'
    
    # Cap score at 4
    score = min(score, 4)
    
    # Determine color and interpretation
    if score == 0:
        color = '🟢'
        interpretation = 'BENIGN'
    elif score == 1:
        color = '🟢'
        interpretation = 'LIKELY BENIGN'
    elif score == 2:
        color = '🟡'
        interpretation = 'UNCERTAIN SIGNIFICANCE'
    elif score == 3:
        color = '🟠'
        interpretation = 'LIKELY PATHOGENIC'
    else:  # score == 4
        color = '🔴'
        interpretation = 'PATHOGENIC'
    
    return {
        'score': score,
        'severity': severity_info['severity'],
        'color': color,
        'interpretation': interpretation,
        'clinical_action': severity_info['clinical_action'],
        'base_score': severity_info['base_score']
    }


def get_clinical_recommendations(score, gene_name, mutation_type):
    """
    Provide clinical recommendations based on pathogenicity
    
    Args:
        score (int): Pathogenicity score
        gene_name (str): Gene name
        mutation_type (str): Type of mutation
    
    Returns:
        list: List of recommendations
    """
    recommendations = []
    
    if score >= 3:
        recommendations.append("🔍 Recommend confirmatory testing with independent method")
        recommendations.append("👨‍⚕️ Genetic counseling strongly advised")
        recommendations.append("👨‍👩‍👧‍👦 Consider family testing if hereditary condition")
        
        # Gene-specific recommendations
        if gene_name in ['BRCA1', 'BRCA2']:
            recommendations.append("📋 Consider enhanced cancer screening protocols")
        elif gene_name == 'TP53':
            recommendations.append("📋 Li-Fraumeni syndrome surveillance protocol may be indicated")
        elif gene_name == 'HBB':
            recommendations.append("🩸 Hematology consultation recommended")
    
    elif score == 2:
        recommendations.append("📊 Additional evidence needed for clinical interpretation")
        recommendations.append("🔬 Consider segregation analysis in family")
        recommendations.append("📚 Search for variant in clinical databases (ClinVar, HGMD)")
    
    else:  # score < 2
        recommendations.append("✅ Likely not clinically significant")
        recommendations.append("📝 Document for future reference")
    
    return recommendations