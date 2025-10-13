# Mutation Impact Analyzer

A beginner-friendly Python tool that analyzes DNA mutations and predicts whether they might cause health problems. Perfect for learning bioinformatics!

**Author:** U. Hema  
**Contact:** hema98661@gmail.com  
**License:** MIT (free to use and modify)

---

## Table of Contents

1. [What Is This Project?](#what-is-this-project)
2. [What Will You Learn?](#what-will-you-learn)
3. [How to Install](#how-to-install)
4. [How to Use](#how-to-use)
5. [Understanding the Results](#understanding-the-results)
6. [Project Files Explained](#project-files-explained)
7. [Mutation Types - Detailed Guide](#mutation-types---detailed-guide)
8. [How the Scoring Works](#how-the-scoring-works)
9. [Troubleshooting](#troubleshooting)
10. [Important Disclaimers](#important-disclaimers)

---

## What Is This Project?

This tool helps you understand **what happens when DNA changes**. 

### Real-World Example:
Imagine you have a DNA sequence like this:
```
ATGCGCCATTAG
```

If one letter (nucleotide) changes from `G` to `A`:
```
ATGCACCATTAG
```

This tool will tell you:
- **What changed**: G became A at position 5
- **What it means**: The protein changed from Arginine to Histidine
- **Is it dangerous?**: Risk score of 4/5 (HIGH risk)
- **What to do**: Recommend genetic counseling

### Why Is This Useful?

- **For Students**: Learn how mutations work in real life
- **For Research**: Understand genetic diseases
- **For Biology Class**: Great for assignments and presentations
- **For Portfolio**: Show coding + biology skills

---

## What Will You Learn?

By using and studying this project, you'll understand:

### 1. Biology Concepts
- **DNA Structure**: How A, T, C, G letters code for life
- **Codons**: Groups of 3 letters that make amino acids
- **Translation**: How DNA becomes proteins
- **Mutations**: Different types and their effects

### 2. Programming Skills
- **Python Basics**: Variables, functions, loops
- **File Handling**: Reading FASTA files (biology data format)
- **Error Handling**: Validating user input
- **Data Visualization**: Creating interactive charts

### 3. Bioinformatics
- **Sequence Analysis**: Processing DNA data
- **Database Queries**: Connecting to medical databases
- **Risk Assessment**: Predicting mutation impact
- **Report Generation**: Creating professional outputs

---

## How to Install

### Step 1: Get Python

Make sure you have Python 3.8 or newer installed:
```bash
python --version
```

If not installed, download from [python.org](https://www.python.org/downloads/)

### Step 2: Download This Project

**Option A - Using Git:**
```bash
git clone https://github.com/Hema291002/mutation-impact-analyzer.git
cd mutation-impact-analyzer
```

**Option B - Download ZIP:**
1. Click the green "Code" button on GitHub
2. Select "Download ZIP"
3. Extract the folder
4. Open terminal/command prompt in that folder

### Step 3: Install Required Packages

```bash
pip install -r requirements.txt
```

This installs:
- **biopython**: For DNA sequence analysis
- **pandas**: For organizing data
- **plotly**: For creating interactive graphs
- **numpy**: For mathematical calculations
- **requests**: For accessing online databases

### Step 4: Verify Installation

```bash
python test_install.py
```

If you see "All libraries installed successfully!" - you're ready!

---

## How to Use

### Basic Usage (Beginner)

1. **Start the program:**
```bash
python analyzer.py
```

2. **Choose input method:**
```
Choose input method:
1. Load from FASTA file
2. Paste sequence directly

Enter choice (1 or 2): 2
```

3. **Enter your gene name:**
```
Enter gene name (e.g., TP53, BRCA1): TP53
```

4. **Enter DNA sequence:**
```
Enter DNA sequence (only A, T, C, G):
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA
```

5. **Choose mutation position:**
```
Enter position to mutate (1-60): 5
```
The program shows you what's at that position.

6. **Enter new nucleotide:**
```
Enter new nucleotide (A/T/C/G): C
```

7. **View results!** The program will:
   - Show detailed analysis on screen
   - Save a text report in `output/` folder
   - Create an interactive HTML dashboard

### Using FASTA Files (Advanced)

FASTA is a standard format for storing DNA sequences.

**Example FASTA file** (`sequences/my_gene.fasta`):
```
>TP53_human_partial
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA
GACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATG
```

To use it:
1. Put your FASTA file in the `sequences/` folder
2. Run `python analyzer.py`
3. Choose option 1
4. Enter your filename: `my_gene.fasta`

---

## Understanding the Results

### Example Output Explained

When you run the analyzer, you get two types of output:

### 1. Text Report (in output/ folder)

```
======================================================================
MUTATION ANALYSIS REPORT
======================================================================
Generated: 2025-10-12 14:30:00

BASIC INFORMATION
----------------------------------------------------------------------
Gene: TP53
Position: 524
Original Nucleotide: G
Mutated Nucleotide: A
```

**What this means:**
- **Gene**: Which gene you're studying (TP53 is a famous cancer gene)
- **Position**: Exact location of the change (524th letter in DNA)
- **Original/Mutated**: What changed (G became A)

```
MUTATION CLASSIFICATION
----------------------------------------------------------------------
Type: Missense
Original Codon: CGC (Arginine)
Mutated Codon: CAC (Histidine)
Amino Acid Change: R175H
Description: Changes Arginine to Histidine
Impact: MODERATE
```

**What this means:**
- **Type**: Missense = one amino acid becomes another
- **Codon**: Group of 3 DNA letters (CGC → CAC)
- **Amino Acid**: Building blocks of proteins (R → H)
- **R175H**: Standard notation (Arginine at position 175 becomes Histidine)

```
PATHOGENICITY ASSESSMENT
----------------------------------------------------------------------
Risk Score: 4/5
Risk Level: HIGH
Classification: LIKELY PATHOGENIC
Note: TP53 is a critical disease/cancer gene
```

**What this means:**
- **Risk Score**: 0=safe, 5=dangerous (4/5 is very concerning)
- **Pathogenic**: Likely to cause disease
- **Critical gene**: TP53 mutations are found in 50% of cancers

```
CLINICAL RECOMMENDATIONS
----------------------------------------------------------------------
1. Genetic counseling strongly recommended
2. Clinical confirmatory testing advised
3. Family history evaluation warranted
```

**What this means:**
- Real-world actions if this was a real patient
- Shows the practical importance of the mutation

### 2. Interactive Dashboard (HTML file)

The HTML file shows:

**Visual Elements:**
- **Sequence Viewer**: See exactly where mutation occurred
- **Codon Comparison**: Before/after side-by-side
- **Risk Gauge**: Visual meter showing danger level
- **Mutation Type Card**: Color-coded classification
- **Disease Associations**: Which diseases are linked
- **Pathway Diagrams**: Biological processes affected

**How to view:**
1. Go to `output/` folder
2. Find file named like `TP53_R175H_dashboard.html`
3. Double-click to open in web browser
4. Interact with charts (hover, zoom, click)

---

## Project Files Explained

### Main Files (The Important Ones)

#### 1. `analyzer.py` - The Brain of the Project

This is the main file that does everything. It contains:

**Key Functions:**

```python
def validate_sequence(sequence):
```
- **Purpose**: Checks if your DNA sequence is valid
- **Checks for**: Only A,T,C,G letters; length divisible by 3
- **Returns**: True/False and error message

```python
def translate_codon(codon):
```
- **Purpose**: Converts 3 DNA letters into 1 amino acid
- **Example**: "ATG" → "M" (Methionine)
- **Uses**: The genetic code table

```python
def classify_mutation(original_codon, mutated_codon):
```
- **Purpose**: Determines mutation type
- **Returns**: Type (Missense/Nonsense/Silent/Stop-loss)
- **Example**: CGC → CAC = Missense

```python
def calculate_pathogenicity_score(mutation_info, gene_name):
```
- **Purpose**: Calculates how dangerous the mutation is
- **Considers**: Mutation type + gene importance
- **Returns**: Score 0-5 and recommendations

```python
def analyze_mutation(sequence, gene_name, position, new_nucleotide):
```
- **Purpose**: Main function that ties everything together
- **Does**: Validates → Classifies → Scores → Reports
- **Returns**: Complete analysis dictionary

#### 2. `pathogenicity.py` - The Risk Calculator

This file determines if mutations are dangerous.

**How it works:**

```python
MUTATION_SEVERITY = {
    'Silent/Synonymous': base_score = 0  # No change in protein
    'Missense': base_score = 2           # One amino acid changes
    'Nonsense': base_score = 4           # Creates STOP signal
    'Stop-loss': base_score = 3          # Removes STOP signal
}
```

**Critical Genes List:**
```python
critical_genes = {
    'TP53',    # Tumor suppressor (>50% of cancers)
    'BRCA1',   # Breast cancer gene
    'BRCA2',   # Breast cancer gene
    'KRAS',    # Pancreatic cancer (90%)
    ...
}
```

If mutation is in a critical gene, score increases!

#### 3. `visualize_results.py` - The Artist

Creates beautiful, interactive visualizations using Plotly.

**What it creates:**
- Sequence viewers with color-coded changes
- Codon comparison bars
- Risk gauge meters (like a speedometer)
- Disease association charts
- Pathway network diagrams

**Key function:**
```python
def create_comprehensive_dashboard():
```
Combines all visualizations into one HTML file.

#### 4. `requirements.txt` - The Shopping List

Lists all Python packages needed:
```
biopython>=1.79
requests>=2.31.0
pandas>=1.5.0
plotly>=5.17.0
numpy>=1.24.0
```

#### 5. `test_install.py` - The Checker

Simple script to verify installation:
```python
try:
    from Bio import SeqIO
    import pandas as pd
    import requests
    import plotly.graph_objects as go
    print("All libraries installed successfully!")
except ImportError as e:
    print(f"Error: {e}")
```

### Folder Structure

```
mutation-impact-analyzer/
│
├── analyzer.py                 # Main program - START HERE
├── pathogenicity.py            # Risk scoring logic
├── visualize_results.py        # Creates visualizations
├── requirements.txt            # Required packages
├── test_install.py             # Installation tester
├── README.md                   # This guide
├── LICENSE                     # MIT license
│
├── sequences/                  # INPUT: Put FASTA files here
│   ├── example_TP53.fasta     # Sample file provided
│   └── .gitkeep               # Keeps folder in Git
│
├── output/                     # OUTPUT: Results appear here
│   ├── *_report.txt           # Text summaries
│   ├── *_dashboard.html       # Interactive visualizations
│   └── .gitkeep               # Keeps folder in Git
│
└── data/                       # OPTIONAL: Extra data files
    └── .gitkeep
```

---

## Mutation Types - Detailed Guide

### 1. Silent (Synonymous) Mutations

**What happens:** DNA changes but protein stays the same

**Why:** Multiple codons can code for the same amino acid

**Example:**
```
Original: CTG → Leucine
Mutated:  CTA → Leucine (still Leucine!)
```

**Risk Level:** LOW (0/5)
- Protein function unchanged
- Usually harmless
- No clinical action needed

**Real-world analogy:** Like changing "color" to "colour" - different spelling, same meaning

### 2. Missense Mutations

**What happens:** One amino acid becomes a different amino acid

**Why:** The codon now codes for a different amino acid

**Example:**
```
Original: CGC → Arginine (R)
Mutated:  CAC → Histidine (H)
Notation: R175H
```

**Risk Level:** MODERATE to HIGH (2-4/5)
- Depends on:
  - How different the amino acids are
  - Where in the protein it occurs
  - How important the gene is

**Types of missense:**

**Conservative:**
- Similar amino acids (both positive, both negative, etc.)
- Example: Leucine → Isoleucine (both hydrophobic)
- Lower risk

**Non-conservative:**
- Very different amino acids
- Example: Arginine (positive) → Histidine (polar)
- Higher risk

**Real-world analogy:** Like replacing a word with another word - might change meaning slightly or completely

### 3. Nonsense Mutations

**What happens:** Creates a premature STOP signal

**Why:** Codon changes to TAA, TAG, or TGA (stop codons)

**Example:**
```
Original: CAG → Glutamine (Q)
Mutated:  TAG → STOP (*)
Result: Protein is cut short!
```

**Risk Level:** VERY HIGH (5/5)
- Protein is incomplete
- Usually non-functional
- Can cause severe diseases

**Consequence:**
```
Normal protein:  MKTF...AQRG (200 amino acids)
Mutated protein: MKTF...AQ (only 50 amino acids!)
```

**Real-world analogy:** Like ending a sentence early - "The quick brown fox jumps over the" - incomplete message!

### 4. Stop-loss Mutations

**What happens:** Removes the normal STOP signal

**Why:** Stop codon (TAA/TAG/TGA) changes to a regular codon

**Example:**
```
Original: TAA → STOP (*)
Mutated:  CAA → Glutamine (Q)
Result: Protein is too long!
```

**Risk Level:** HIGH (4/5)
- Protein extends beyond normal
- Extra amino acids added
- Can disrupt function

**Consequence:**
```
Normal protein:  MKTF...AQRG-STOP (200 amino acids)
Mutated protein: MKTF...AQRGQKL... (250+ amino acids!)
```

**Real-world analogy:** Like a sentence that doesn't end - keeps going into the next paragraph where it shouldn't!

---

## How the Scoring Works

### Risk Score Calculation (0-5 scale)

The program uses a smart algorithm:

#### Base Scores:
```
Silent mutation:    0 points
Missense mutation:  2 points
Stop-loss mutation: 3 points
Nonsense mutation:  4 points
```

#### Modifiers:

**1. Gene Importance (+0 to +1 points)**

Critical genes get higher scores:
```python
critical_genes = ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'BRAF', ...]
```

Example:
- Missense in unknown gene: 2 points
- Missense in TP53: 3 or 4 points (upgraded!)

**2. Amino Acid Properties (+0 to +1 points)**

For missense mutations:

**Conservative change (similar amino acids):**
```
Leucine → Isoleucine: +0 points
Both are hydrophobic, similar size
```

**Non-conservative change (different amino acids):**
```
Arginine → Histidine: +1 point
Positive charge → Polar, big change!
```

#### Final Classification:

```
Score 0:     BENIGN              (Green)
Score 1-2:   LIKELY BENIGN       (Light Green)
Score 3:     UNCERTAIN           (Yellow)
Score 4:     LIKELY PATHOGENIC   (Orange)
Score 5:     PATHOGENIC          (Red)
```

### Example Calculations:

**Example 1: Silent mutation in any gene**
```
Base score:        0
Gene modifier:     0
Property modifier: 0
TOTAL:            0/5 → BENIGN
```

**Example 2: Missense in regular gene**
```
Base score:        2
Gene modifier:     0
Property modifier: 0 (conservative)
TOTAL:            2/5 → LIKELY BENIGN
```

**Example 3: Missense in TP53**
```
Base score:        2
Gene modifier:     +2 (critical gene)
Property modifier: 0
TOTAL:            4/5 → LIKELY PATHOGENIC
```

**Example 4: Nonsense in any gene**
```
Base score:        4
Gene modifier:     +1 (any gene)
Property modifier: N/A
TOTAL:            5/5 → PATHOGENIC
```

---

## Troubleshooting

### Problem 1: "Module not found" or "ImportError"

**Error message:**
```
ModuleNotFoundError: No module named 'Bio'
```

**Solution:**
```bash
pip install -r requirements.txt
```

**If still not working:**
```bash
pip install biopython pandas plotly numpy requests
```

**Still issues? Try:**
```bash
python -m pip install --upgrade pip
pip install biopython pandas plotly numpy requests
```

### Problem 2: "Invalid sequence" error

**Error message:**
```
Error: Sequence contains invalid characters: {'X', 'N'}
```

**Solution:**
- Only use A, T, C, G letters
- Remove any spaces or special characters
- Check for typos

**Example - Wrong:**
```
ATGCXGTAA  ← X is not valid
```

**Example - Correct:**
```
ATGCGGTAA
```

### Problem 3: "Not divisible by 3" error

**Error message:**
```
Sequence length (61) is not divisible by 3. Incomplete codons detected.
```

**Solution:**
DNA must be in groups of 3 (codons). Either:
- Add missing nucleotides to complete the last codon
- Remove extra nucleotides

**Example - Wrong (61 letters):**
```
ATGCGTCGATAG  ← 61 is not divisible by 3
```

**Example - Correct (60 letters):**
```
ATGCGTCGATA  ← Removed one letter, now 60 (20 codons)
```

**Or add letters:**
```
ATGCGTCGATAGC  ← Added two letters, now 63 (21 codons)
```

### Problem 4: "File not found" error

**Error message:**
```
FileNotFoundError: FASTA file not found: sequences/my_file.fasta
```

**Solution:**
1. Make sure file is in `sequences/` folder
2. Check spelling of filename (case-sensitive!)
3. Include `.fasta` extension

**Correct path:**
```
sequences/my_gene.fasta  ← File must be here
```

### Problem 5: Position out of range

**Error message:**
```
Position 100 is out of range (1-60)
```

**Solution:**
- Check your sequence length
- Choose a position within 1 to length
- Remember: positions start at 1, not 0!

**Example:**
```
Sequence: ATGCGT (6 letters)
Valid positions: 1, 2, 3, 4, 5, 6
Invalid: 0, 7, 8, ...
```

### Problem 6: Same nucleotide error

**Error message:**
```
Position 5 already contains nucleotide C
```

**Solution:**
- You're trying to "mutate" to the same letter
- Choose a different nucleotide

**Example - Wrong:**
```
Position 5: C
Your input: C  ← No change!
```

**Example - Correct:**
```
Position 5: C
Your input: A  ← Now it's a real mutation
```

### Problem 7: Plotly visualizations not showing

**Solution:**
```bash
pip install --upgrade plotly
```

Or use specific version:
```bash
pip install plotly==5.17.0
```

### Problem 8: Python version issues

**Check your version:**
```bash
python --version
```

**Need Python 3.8 or higher!**

**If too old:**
- Download latest Python from python.org
- Or use: `python3` instead of `python`

---

## Important Disclaimers

### Educational Purpose Only

This project is designed for:
- Learning bioinformatics concepts
- Understanding how mutations work
- Practicing Python programming
- School/college assignments
- Building your portfolio

### NOT for Clinical Use

This tool is NOT intended for:
- Medical diagnosis
- Treatment decisions
- Clinical genetic counseling
- Patient care

### Why?

Real clinical analysis requires:
- Validated laboratory methods
- Multiple confirmation tests
- Professional genetic counselors
- FDA-approved software
- Patient medical history
- Family pedigree analysis

### Database Limitations

Some features use simulated data:
- ClinVar queries (commented in code)
- COSMIC data (limited dataset)
- OMIM information (basic examples)

**For real research:**
- Use official APIs with proper authentication
- Access full databases through institutions
- Follow database usage policies

### Always Verify

If using for serious work:
- Cross-check with multiple tools
- Consult with advisors/professors
- Cite your sources properly
- Validate results independently

### Ethical Considerations

When working with genetic data:
- Respect privacy (don't use real patient data)
- Follow institutional guidelines
- Understand ethical implications
- Get proper permissions for research

---

## How to Cite This Project

If you use this tool in your work:

**For Papers/Reports:**
```
Hema, U. (2025). Mutation Impact Analyzer: A Python tool for 
genetic mutation analysis. GitHub repository. 
https://github.com/Hema291002/mutation-impact-analyzer
```

**For Presentations:**
```
Tool: Mutation Impact Analyzer
Author: U. Hema
GitHub: github.com/Hema291002/mutation-impact-analyzer
```

---

## Learning Resources

Want to learn more? Check these out:

### Online Courses (Free)
- **Coursera**: "Bioinformatics Specialization"
- **edX**: "Introduction to Genomic Technologies"
- **Khan Academy**: "DNA and Genetics"

### Textbooks
- "Introduction to Bioinformatics" by Arthur Lesk
- "Python for Biologists" by Martin Jones
- "Molecular Biology of the Cell" by Alberts et al.

### Websites
- **NCBI**: National Center for Biotechnology Information
- **ClinVar**: Clinical genetic variants
- **COSMIC**: Catalogue Of Somatic Mutations In Cancer
- **UniProt**: Protein database

### Practice Datasets
- Use example files in `sequences/` folder
- Download from NCBI Gene database
- Try mutations in different genes

---

## Future Improvements

Ideas for expanding this project:

### Easy Additions:
- [ ] More gene databases
- [ ] Additional visualization types
- [ ] Export results to PDF
- [ ] Batch processing multiple mutations

### Intermediate:
- [ ] Web interface using Flask
- [ ] Multiple sequence alignment
- [ ] Protein structure prediction
- [ ] Machine learning predictions

### Advanced:
- [ ] Real-time database connections
- [ ] Integration with protein databases
- [ ] 3D protein structure visualization
- [ ] Population frequency analysis

---

## Contributing

Want to improve this project?

### How to contribute:
1. Fork the repository
2. Create a new branch: `git checkout -b feature-name`
3. Make your changes
4. Test thoroughly
5. Submit a pull request

### Ideas welcome for:
- Bug fixes
- New features
- Documentation improvements
- Code optimization
- Additional examples

---

## License

MIT License - Free and open source!

**This means you can:**
- Use it for personal projects
- Use it for school assignments
- Modify the code
- Share with friends
- Learn from it

**Just remember to:**
- Keep the license notice
- Give credit to original author
- Not claim it as entirely your own work

---

## Contact & Support

**Author:** U. Hema  
**Email:** hema98661@gmail.com  
**GitHub:** github.com/Hema291002/mutation-impact-analyzer

### How to get help:

1. **Check this README first** - Most answers are here!
2. **Look at existing issues** - Someone might have had the same problem
3. **Open a new issue** - Describe your problem clearly
4. **Email me** - For specific questions

### When reporting bugs:
- Include error message
- Describe what you were trying to do
- Share your sequence (if not confidential)
- Mention your Python version
- Include operating system (Windows/Mac/Linux)

---

**Happy mutation analyzing! Remember: Every expert was once a beginner.**

**Learn. Build. Share. Repeat.**