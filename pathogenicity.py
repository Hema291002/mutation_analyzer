"""
Pathogenicity Prediction Module
Predicts the pathogenicity of mutations based on multiple factors
"""


class PathogenicityPredictor:
    """
    Predicts pathogenicity of mutations using rule-based scoring
    """
    
    def __init__(self):
        """Initialize the predictor"""
        self.weights = {
            'mutation_type': 3.0,
            'aa_property_change': 2.5,
            'position': 1.5,
            'cosmic': 2.0,
            'conservation': 1.0
        }
        
        # Amino acid properties for similarity scoring
        self.aa_properties = {
            'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'],
            'polar': ['S', 'T', 'C', 'Y', 'N', 'Q'],
            'charged_positive': ['K', 'R', 'H'],
            'charged_negative': ['D', 'E'],
            'special': ['G'],
            'stop': ['*']
        }
    
    def get_aa_property(self, aa):
        """Get the property category of an amino acid"""
        for prop, amino_acids in self.aa_properties.items():
            if aa in amino_acids:
                return prop
        return 'unknown'
    
    def predict(self, original_aa, mutated_aa, position, gene_name='', cosmic_found=False):
        """
        Predict pathogenicity of a mutation
        
        Args:
            original_aa: Original amino acid (single letter)
            mutated_aa: Mutated amino acid (single letter)
            position: Position in the protein (codon number)
            gene_name: Name of the gene
            cosmic_found: Whether mutation is found in COSMIC database
            
        Returns:
            dict: Pathogenicity prediction results
        """
        score = 0.0
        factors = []
        
        # 1. Mutation Type Scoring
        if mutated_aa == '*':
            # Nonsense mutation - creates stop codon
            score += 4.0
            factors.append("Nonsense mutation (premature stop codon) - HIGH IMPACT")
        elif original_aa == '*':
            # Nonstop mutation
            score += 3.5
            factors.append("Nonstop mutation (stop codon lost) - HIGH IMPACT")
        elif original_aa == mutated_aa:
            # Synonymous - no change
            score += 0.0
            factors.append("Synonymous mutation (no amino acid change) - BENIGN")
        else:
            # Missense mutation
            orig_prop = self.get_aa_property(original_aa)
            mut_prop = self.get_aa_property(mutated_aa)
            
            if orig_prop == mut_prop:
                # Conservative substitution
                score += 1.5
                factors.append(f"Conservative missense ({orig_prop} to {mut_prop}) - MODERATE")
            else:
                # Non-conservative substitution
                score += 3.0
                factors.append(f"Non-conservative missense ({orig_prop} to {mut_prop}) - HIGH")
                
                # Extra penalty for specific property changes
                if (orig_prop in ['charged_positive', 'charged_negative'] and 
                    mut_prop == 'hydrophobic'):
                    score += 0.5
                    factors.append("Charge to hydrophobic change - additional penalty")
        
        # 2. COSMIC Database Factor
        if cosmic_found:
            score += 2.0
            factors.append("Found in COSMIC cancer database - SIGNIFICANT")
        
        # 3. Position-based scoring (simplified - domains typically critical)
        # Known critical domains for common cancer genes
        critical_domains = {
            'TP53': [(100, 300)],  # DNA-binding domain
            'BRCA1': [(1, 100), (1560, 1863)],  # RING and BRCT domains
            'EGFR': [(712, 979)],  # Kinase domain
        }
        
        if gene_name.upper() in critical_domains:
            for start, end in critical_domains[gene_name.upper()]:
                if start <= position <= end:
                    score += 1.5
                    factors.append(f"Located in critical functional domain (position {position})")
                    break
        
        # 4. Cap the maximum score
        score = min(score, 10.0)
        
        # 5. Determine risk level and classification
        if score >= 7.0:
            risk_level = "VERY HIGH"
            classification = "Likely Pathogenic"
            color = "red"
        elif score >= 5.0:
            risk_level = "HIGH"
            classification = "Possibly Pathogenic"
            color = "orange"
        elif score >= 3.0:
            risk_level = "MODERATE"
            classification = "Uncertain Significance"
            color = "yellow"
        elif score >= 1.0:
            risk_level = "LOW"
            classification = "Likely Benign"
            color = "lightgreen"
        else:
            risk_level = "VERY LOW"
            classification = "Benign"
            color = "green"
        
        # 6. Generate recommendations
        recommendations = self._generate_recommendations(
            score, classification, original_aa, mutated_aa, cosmic_found
        )
        
        return {
            'score': round(score, 2),
            'risk_level': risk_level,
            'classification': classification,
            'color': color,
            'factors': factors,
            'recommendations': recommendations
        }
    
    def _generate_recommendations(self, score, classification, original_aa, 
                                 mutated_aa, cosmic_found):
        """Generate clinical recommendations based on prediction"""
        recommendations = []
        
        if score >= 7.0:
            recommendations.append("Urgent: Recommend immediate functional validation studies")
            recommendations.append("Clinical genetic counseling strongly advised")
            if cosmic_found:
                recommendations.append("Review patient and family cancer history")
            recommendations.append("Consider additional molecular testing")
        
        elif score >= 5.0:
            recommendations.append("Recommend functional validation in relevant cell models")
            recommendations.append("Genetic counseling advised")
            recommendations.append("Monitor patient closely for clinical manifestations")
            if cosmic_found:
                recommendations.append("Consider tumor sequencing if applicable")
        
        elif score >= 3.0:
            recommendations.append("Uncertain significance - require additional evidence")
            recommendations.append("Consider segregation analysis in family members")
            recommendations.append("Recommend periodic re-evaluation as databases are updated")
            recommendations.append("Functional studies may help clarify significance")
        
        elif score >= 1.0:
            recommendations.append("Likely benign - minimal clinical concern")
            recommendations.append("Routine monitoring may be sufficient")
            recommendations.append("Re-evaluate if additional clinical symptoms develop")
        
        else:
            recommendations.append("Benign variant - no clinical action required")
            recommendations.append("Standard population screening guidelines apply")
        
        # Add mutation-specific recommendations
        if mutated_aa == '*':
            recommendations.append("Nonsense-mediated decay may occur - check protein expression")
        elif original_aa == '*':
            recommendations.append("Extended protein may have altered localization or stability")
        
        return recommendations


def predict_pathogenicity(original_aa, mutated_aa, position, gene_name='', 
                         cosmic_found=False):
    """
    Convenience function for pathogenicity prediction
    
    Args:
        original_aa: Original amino acid
        mutated_aa: Mutated amino acid
        position: Position in protein
        gene_name: Gene name
        cosmic_found: Whether found in COSMIC
        
    Returns:
        dict: Prediction results
    """
    predictor = PathogenicityPredictor()
    return predictor.predict(original_aa, mutated_aa, position, 
                            gene_name, cosmic_found)


# Example usage
if __name__ == "__main__":
    # Test the predictor
    predictor = PathogenicityPredictor()
    
    # Test Case 1: Known pathogenic TP53 mutation
    result1 = predictor.predict(
        original_aa='R',
        mutated_aa='H',
        position=175,
        gene_name='TP53',
        cosmic_found=True
    )
    
    print("Test Case 1: TP53 R175H (known pathogenic)")
    print(f"Score: {result1['score']}/10")
    print(f"Risk: {result1['risk_level']}")
    print(f"Classification: {result1['classification']}")
    print(f"Factors: {result1['factors']}")
    print()
    
    # Test Case 2: Synonymous mutation
    result2 = predictor.predict(
        original_aa='L',
        mutated_aa='L',
        position=100,
        gene_name='BRCA1',
        cosmic_found=False
    )
    
    print("Test Case 2: BRCA1 L100L (synonymous)")
    print(f"Score: {result2['score']}/10")
    print(f"Risk: {result2['risk_level']}")
    print(f"Classification: {result2['classification']}")
    print()
    
    # Test Case 3: Nonsense mutation
    result3 = predictor.predict(
        original_aa='Q',
        mutated_aa='*',
        position=50,
        gene_name='TP53',
        cosmic_found=False
    )
    
    print("Test Case 3: TP53 Q50* (nonsense)")
    print(f"Score: {result3['score']}/10")
    print(f"Risk: {result3['risk_level']}")
    print(f"Classification: {result3['classification']}")
