"""
Validation Script for Mutation Impact Analyzer
Tests the fixed analyzer with proper stop codon recognition
"""

from analyzer import (
    validate_sequence,
    translate_codon,
    classify_mutation,
    calculate_pathogenicity_score,
    is_stop_codon
)

def test_stop_codon_recognition():
    """Test that stop codons are properly recognized."""
    print("\n" + "="*80)
    print("TEST 1: Stop Codon Recognition")
    print("="*80)
    
    stop_codons = ['TAA', 'TAG', 'TGA']
    passed = 0
    failed = 0
    
    for codon in stop_codons:
        result = is_stop_codon(codon)
        if result:
            print(f"PASS: {codon} recognized as stop codon")
            passed += 1
        else:
            print(f"FAIL: {codon} NOT recognized as stop codon")
            failed += 1
    
    # Test non-stop codons
    non_stop = ['ATG', 'CGC', 'TTT']
    for codon in non_stop:
        result = is_stop_codon(codon)
        if not result:
            print(f"PASS: {codon} correctly NOT identified as stop codon")
            passed += 1
        else:
            print(f"FAIL: {codon} incorrectly identified as stop codon")
            failed += 1
    
    return passed, failed


def test_nonsense_mutation_classification():
    """Test that nonsense mutations are correctly classified."""
    print("\n" + "="*80)
    print("TEST 2: Nonsense Mutation Classification")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test cases: (original_codon, mutated_codon, expected_type)
    test_cases = [
        ('CAG', 'TAG', 'Nonsense', 'Q to STOP'),
        ('TGG', 'TGA', 'Nonsense', 'W to STOP'),
        ('TAC', 'TAA', 'Nonsense', 'Y to STOP'),
    ]
    
    for original, mutated, expected_type, description in test_cases:
        result = classify_mutation(original, mutated)
        
        if result['type'] == expected_type:
            print(f"PASS: {original} -> {mutated} ({description})")
            print(f"  Classified as: {result['type']}")
            print(f"  Description: {result['description']}")
            passed += 1
        else:
            print(f"FAIL: {original} -> {mutated} ({description})")
            print(f"  Expected: {expected_type}")
            print(f"  Got: {result['type']}")
            failed += 1
    
    return passed, failed


def test_stop_loss_mutation_classification():
    """Test that stop-loss mutations are correctly classified."""
    print("\n" + "="*80)
    print("TEST 3: Stop-loss Mutation Classification")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test cases: (original_codon, mutated_codon, expected_type)
    test_cases = [
        ('TAA', 'CAA', 'Stop-loss', 'STOP to Q'),
        ('TAG', 'CAG', 'Stop-loss', 'STOP to Q'),
        ('TGA', 'TCA', 'Stop-loss', 'STOP to S'),
    ]
    
    for original, mutated, expected_type, description in test_cases:
        result = classify_mutation(original, mutated)
        
        if result['type'] == expected_type:
            print(f"PASS: {original} -> {mutated} ({description})")
            print(f"  Classified as: {result['type']}")
            print(f"  Description: {result['description']}")
            passed += 1
        else:
            print(f"FAIL: {original} -> {mutated} ({description})")
            print(f"  Expected: {expected_type}")
            print(f"  Got: {result['type']}")
            failed += 1
    
    return passed, failed


def test_missense_mutation_classification():
    """Test that missense mutations are correctly classified."""
    print("\n" + "="*80)
    print("TEST 4: Missense Mutation Classification")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test cases: (original_codon, mutated_codon, expected_type)
    test_cases = [
        ('CGC', 'CAC', 'Missense', 'R to H - TP53 R175H'),
        ('GGT', 'GAT', 'Missense', 'G to D - KRAS G12D'),
        ('GTG', 'ATG', 'Missense', 'V to M - BRAF V600M'),
    ]
    
    for original, mutated, expected_type, description in test_cases:
        result = classify_mutation(original, mutated)
        
        if result['type'] == expected_type:
            print(f"PASS: {original} -> {mutated} ({description})")
            print(f"  Classified as: {result['type']}")
            passed += 1
        else:
            print(f"FAIL: {original} -> {mutated} ({description})")
            print(f"  Expected: {expected_type}")
            print(f"  Got: {result['type']}")
            failed += 1
    
    return passed, failed


def test_silent_mutation_classification():
    """Test that silent mutations are correctly classified."""
    print("\n" + "="*80)
    print("TEST 5: Silent Mutation Classification")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test cases: (original_codon, mutated_codon, expected_type)
    test_cases = [
        ('CTG', 'CTA', 'Silent', 'L to L'),
        ('CGT', 'CGC', 'Silent', 'R to R'),
        ('TCT', 'TCC', 'Silent', 'S to S'),
    ]
    
    for original, mutated, expected_type, description in test_cases:
        result = classify_mutation(original, mutated)
        
        if result['type'] == expected_type:
            print(f"PASS: {original} -> {mutated} ({description})")
            print(f"  Classified as: {result['type']}")
            passed += 1
        else:
            print(f"FAIL: {original} -> {mutated} ({description})")
            print(f"  Expected: {expected_type}")
            print(f"  Got: {result['type']}")
            failed += 1
    
    return passed, failed


def test_pathogenicity_scoring():
    """Test pathogenicity scoring for different mutation types."""
    print("\n" + "="*80)
    print("TEST 6: Pathogenicity Scoring")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test nonsense mutation in critical gene
    mutation_info = classify_mutation('CAG', 'TAG')
    pathogenicity = calculate_pathogenicity_score(mutation_info, 'TP53')
    
    print(f"\nNonsense mutation in TP53:")
    print(f"  Score: {pathogenicity['score']}/5")
    print(f"  Risk Level: {pathogenicity['risk_level']}")
    print(f"  Classification: {pathogenicity['classification']}")
    
    if pathogenicity['score'] == 5:
        print("  PASS: Correct score (5/5)")
        passed += 1
    else:
        print(f"  FAIL: Expected 5, got {pathogenicity['score']}")
        failed += 1
    
    # Test missense mutation in critical gene
    mutation_info = classify_mutation('CGC', 'CAC')
    pathogenicity = calculate_pathogenicity_score(mutation_info, 'TP53')
    
    print(f"\nMissense mutation in TP53 (R175H):")
    print(f"  Score: {pathogenicity['score']}/5")
    print(f"  Risk Level: {pathogenicity['risk_level']}")
    print(f"  Classification: {pathogenicity['classification']}")
    
    if pathogenicity['score'] == 4:
        print("  PASS: Correct score (4/5)")
        passed += 1
    else:
        print(f"  FAIL: Expected 4, got {pathogenicity['score']}")
        failed += 1
    
    # Test silent mutation
    mutation_info = classify_mutation('CTG', 'CTA')
    pathogenicity = calculate_pathogenicity_score(mutation_info, 'TP53')
    
    print(f"\nSilent mutation:")
    print(f"  Score: {pathogenicity['score']}/5")
    print(f"  Risk Level: {pathogenicity['risk_level']}")
    print(f"  Classification: {pathogenicity['classification']}")
    
    if pathogenicity['score'] == 0:
        print("  PASS: Correct score (0/5)")
        passed += 1
    else:
        print(f"  FAIL: Expected 0, got {pathogenicity['score']}")
        failed += 1
    
    # Test stop-loss mutation
    mutation_info = classify_mutation('TAA', 'CAA')
    pathogenicity = calculate_pathogenicity_score(mutation_info, 'TP53')
    
    print(f"\nStop-loss mutation:")
    print(f"  Score: {pathogenicity['score']}/5")
    print(f"  Risk Level: {pathogenicity['risk_level']}")
    print(f"  Classification: {pathogenicity['classification']}")
    
    if pathogenicity['score'] == 4:
        print("  PASS: Correct score (4/5)")
        passed += 1
    else:
        print(f"  FAIL: Expected 4, got {pathogenicity['score']}")
        failed += 1
    
    return passed, failed


def test_sequence_validation():
    """Test sequence validation including stop codons."""
    print("\n" + "="*80)
    print("TEST 7: Sequence Validation")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Valid sequence with stop codon
    sequence1 = "ATGCGTCGATAA"  # Has stop codon TAA at end
    is_valid, error = validate_sequence(sequence1)
    
    if is_valid:
        print(f"PASS: Sequence with stop codon accepted")
        print(f"  Sequence: {sequence1}")
        passed += 1
    else:
        print(f"FAIL: Sequence with stop codon rejected")
        print(f"  Error: {error}")
        failed += 1
    
    # Invalid sequence - not divisible by 3
    sequence2 = "ATGCGTCG"  # 8 nucleotides
    is_valid, error = validate_sequence(sequence2)
    
    if not is_valid and "not divisible by 3" in error:
        print(f"PASS: Incomplete codon correctly rejected")
        print(f"  Sequence: {sequence2}")
        passed += 1
    else:
        print(f"FAIL: Incomplete codon not detected")
        failed += 1
    
    # Invalid sequence - contains invalid characters
    sequence3 = "ATGCXGTAA"
    is_valid, error = validate_sequence(sequence3)
    
    if not is_valid and "invalid" in error.lower():
        print(f"PASS: Invalid characters correctly rejected")
        print(f"  Sequence: {sequence3}")
        passed += 1
    else:
        print(f"FAIL: Invalid characters not detected")
        failed += 1
    
    return passed, failed


def test_translation():
    """Test codon translation including stop codons."""
    print("\n" + "="*80)
    print("TEST 8: Codon Translation")
    print("="*80)
    
    passed = 0
    failed = 0
    
    # Test stop codons
    stop_codons = [('TAA', '*'), ('TAG', '*'), ('TGA', '*')]
    
    for codon, expected_aa in stop_codons:
        try:
            result = translate_codon(codon)
            if result == expected_aa:
                print(f"PASS: {codon} -> {result} (STOP)")
                passed += 1
            else:
                print(f"FAIL: {codon} -> {result} (expected {expected_aa})")
                failed += 1
        except Exception as e:
            print(f"FAIL: {codon} raised exception: {e}")
            failed += 1
    
    # Test regular codons
    regular_codons = [('ATG', 'M'), ('CGC', 'R'), ('GGT', 'G')]
    
    for codon, expected_aa in regular_codons:
        try:
            result = translate_codon(codon)
            if result == expected_aa:
                print(f"PASS: {codon} -> {result}")
                passed += 1
            else:
                print(f"FAIL: {codon} -> {result} (expected {expected_aa})")
                failed += 1
        except Exception as e:
            print(f"FAIL: {codon} raised exception: {e}")
            failed += 1
    
    return passed, failed


def run_all_tests():
    """Run all validation tests."""
    print("\n" + "="*80)
    print("MUTATION IMPACT ANALYZER - VALIDATION SUITE")
    print("="*80)
    print("Testing fixed version with proper stop codon recognition")
    
    total_passed = 0
    total_failed = 0
    
    # Run all tests
    tests = [
        test_stop_codon_recognition,
        test_translation,
        test_sequence_validation,
        test_silent_mutation_classification,
        test_missense_mutation_classification,
        test_nonsense_mutation_classification,
        test_stop_loss_mutation_classification,
        test_pathogenicity_scoring,
    ]
    
    for test_func in tests:
        try:
            passed, failed = test_func()
            total_passed += passed
            total_failed += failed
        except Exception as e:
            print(f"\nERROR in {test_func.__name__}: {e}")
            import traceback
            traceback.print_exc()
            total_failed += 1
    
    # Print summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    print(f"Total tests run: {total_passed + total_failed}")
    print(f"Passed: {total_passed}")
    print(f"Failed: {total_failed}")
    
    if total_failed == 0:
        print("\nSTATUS: ALL TESTS PASSED!")
        success_rate = 100.0
    else:
        success_rate = (total_passed / (total_passed + total_failed)) * 100
        print(f"\nSTATUS: {success_rate:.1f}% tests passed")
        print("Some issues found - review failed tests above")
    
    print("="*80)
    
    # Specific check for the original bug
    print("\n" + "="*80)
    print("CRITICAL BUG CHECK: Nonsense Mutation Detection")
    print("="*80)
    
    print("\nTesting CAG -> TAG (Q to STOP):")
    result = classify_mutation('CAG', 'TAG')
    print(f"  Classification: {result['type']}")
    print(f"  Description: {result['description']}")
    
    if result['type'] == 'Nonsense':
        print("  STATUS: BUG FIXED - Nonsense mutations now detected correctly!")
    else:
        print("  STATUS: BUG STILL PRESENT - Incorrectly classified as Missense")
    
    print("="*80)
    
    return total_passed, total_failed


if __name__ == "__main__":
    run_all_tests()