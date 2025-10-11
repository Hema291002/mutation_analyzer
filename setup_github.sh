#!/bin/bash
# Simple GitHub Setup - Creates 5 Essential Files
# Save this as: setup_github.sh
# Run with: bash setup_github.sh

echo "🧬 Creating essential files for GitHub upload..."
echo ""

# Get user information
read -p "Enter your name (for LICENSE): " USER_NAME
read -p "Enter your email: " USER_EMAIL
read -p "Enter your GitHub username: " GITHUB_USER

echo ""
echo "Creating files..."
echo ""

# ============================================
# 1. CREATE .gitignore
# ============================================
echo "📄 Creating .gitignore..."

cat > .gitignore << 'EOF'
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Virtual Environment
venv/
ENV/
env/
.venv

# IDE
.vscode/
.idea/
*.swp
*.swo
*~
.DS_Store

# Jupyter Notebook
.ipynb_checkpoints

# Output files (keep folder, ignore contents)
output/*.html
output/*.txt
output/*.png
output/*.pdf

# Keep output folder but ignore contents
!output/.gitkeep

# Data files (large sequences)
data/*.fasta
sequences/*.fasta

# Keep example files
!sequences/example_TP53.fasta
!data/README.md

# Testing
.pytest_cache/
.coverage
htmlcov/

# Environment variables
.env
EOF

echo "✅ .gitignore created"

# ============================================
# 2. CREATE LICENSE
# ============================================
echo "📄 Creating LICENSE..."

YEAR=$(date +%Y)

cat > LICENSE << EOF
MIT License

Copyright (c) $YEAR $USER_NAME

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF

echo "✅ LICENSE created"

# ============================================
# 3. CREATE README.md
# ============================================
echo "📄 Creating README.md..."

cat > README.md << EOF
# 🧬 Mutation Impact Analyzer

![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![BioPython](https://img.shields.io/badge/BioPython-1.79%2B-orange.svg)

**A comprehensive bioinformatics tool for analyzing genetic mutations, predicting disease associations, and mapping biological pathways with interactive visualizations.**

## ✨ Features

- 🔬 **Mutation Classification** - Automatically identify Missense, Nonsense, Silent, and Stop-loss mutations
- 🗄️ **Database Integration** - Real-time queries to ClinVar, COSMIC, and OMIM
- 🧬 **Pathway Mapping** - Connects mutations to KEGG and Reactome biological pathways
- 📊 **Interactive Visualizations** - Beautiful HTML dashboards powered by Plotly
- 🎯 **Pathogenicity Prediction** - Intelligent risk scoring system (0-5 scale)
- 📁 **FASTA Support** - Industry-standard sequence format

## 🚀 Installation

\`\`\`bash
# Clone the repository
git clone https://github.com/$GITHUB_USER/mutation-impact-analyzer.git
cd mutation-impact-analyzer

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\\Scripts\\activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
python test_install.py
\`\`\`

## 🎯 Quick Start

\`\`\`bash
python analyzer.py
\`\`\`

Follow the interactive prompts:
1. Choose input method (FASTA file or paste sequence)
2. Enter gene name (e.g., TP53, BRCA1, KRAS)
3. Specify nucleotide position to mutate
4. Enter new nucleotide (A/T/C/G)
5. View results and generate visualizations

## 📖 Example Analysis

\`\`\`
Gene: TP53
Position: 524 (codon 175)
Original: CGC (Arginine)
Mutated: CAC (Histidine)
Type: Missense (R175H)
Pathogenicity: LIKELY PATHOGENIC
Risk: HIGH (5/5)
\`\`\`

## 🧬 Supported Genes

### Cancer Genes (COSMIC)
- **TP53** - Tumor suppressor (>50% of cancers)
- **BRCA1/BRCA2** - DNA repair genes
- **KRAS** - GTPase signaling (90% pancreatic cancer)
- **EGFR** - Growth factor receptor
- **BRAF** - Kinase (50% melanoma)
- **PIK3CA** - PI3K pathway

### Disease Genes (OMIM)
- **TP53** - Li-Fraumeni syndrome
- **BRCA1/2** - Hereditary breast/ovarian cancer
- **CFTR** - Cystic fibrosis
- **HBB** - Sickle cell anemia
- **DMD** - Duchenne muscular dystrophy

## 📁 Project Structure

\`\`\`
mutation-impact-analyzer/
├── analyzer.py                 # Main analysis engine
├── visualize_results.py        # Plotly visualization generator
├── pathogenicity.py           # Risk assessment module
├── requirements.txt           # Python dependencies
├── test_install.py            # Installation verification
├── README.md                  # This file
├── LICENSE                    # MIT License
│
├── sequences/                 # Input FASTA files
│   └── example_TP53.fasta
│
├── output/                    # Generated reports
│   ├── *_dashboard.html      # Interactive visualizations
│   └── *_report.txt          # Text summaries
│
└── data/                      # Additional data
\`\`\`

## 🔬 How It Works

1. **Sequence Processing** - Load from FASTA or direct input, validate composition
2. **Classification** - Translate DNA to protein, compare amino acids
3. **Database Queries** - Search ClinVar, COSMIC, OMIM for clinical significance
4. **Risk Assessment** - Score based on mutation type and gene importance
5. **Visualization** - Generate interactive Plotly charts and HTML dashboards

## 📊 Output

The tool generates two types of output:

### Text Reports (\`output/*_report.txt\`)
- Mutation summary
- Molecular impact
- Clinical significance
- Recommendations

### Interactive HTML Dashboard (\`output/*_dashboard.html\`)
- Summary cards with key metrics
- Interactive sequence viewer
- Codon change visualization
- Pathogenicity gauge (0-5 scale)
- Pathway network diagrams
- Disease association charts

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (\`git checkout -b feature/AmazingFeature\`)
3. Commit your changes (\`git commit -m 'Add some AmazingFeature'\`)
4. Push to the branch (\`git push origin feature/AmazingFeature\`)
5. Open a Pull Request

## ⚠️ Disclaimer

**This tool is for research and educational purposes only.**

- ❌ Not intended for clinical diagnosis
- ❌ Requires validation by qualified professionals
- ✅ Always consult genetic counselors for clinical decisions

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This project uses:
- **[BioPython](https://biopython.org/)** - Sequence analysis toolkit
- **[Plotly](https://plotly.com/)** - Interactive visualization
- **[NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** - Clinical variants
- **[COSMIC](https://cancer.sanger.ac.uk/cosmic)** - Cancer mutations
- **[OMIM](https://www.omim.org/)** - Genetic diseases
- **[KEGG](https://www.genome.jp/kegg/)** - Pathway database

## 📞 Contact

- **Issues**: [GitHub Issues](https://github.com/$GITHUB_USER/mutation-impact-analyzer/issues)
- **Email**: $USER_EMAIL

## 🌟 Support This Project

If you find this tool useful:
- ⭐ Star this repository
- 🐛 Report bugs and issues
- 💡 Suggest new features
- 🤝 Contribute code

---

**Made with 🧬 and ❤️ for the bioinformatics community**
EOF

echo "✅ README.md created"

# ============================================
# 4. RENAME TEST FILE
# ============================================
echo "📄 Checking test file..."

if [ -f "test_installl.py" ]; then
    mv test_installl.py test_install.py
    echo "✅ Renamed test_installl.py → test_install.py"
elif [ -f "test_install.py" ]; then
    echo "✅ test_install.py already exists"
else
    echo "⚠️  Warning: No test file found"
fi

# ============================================
# 5. CREATE FOLDER STRUCTURE
# ============================================
echo "📁 Creating folder structure..."

mkdir -p sequences
mkdir -p output
mkdir -p data

touch sequences/.gitkeep
touch output/.gitkeep
touch data/.gitkeep

echo "✅ Folders created with .gitkeep files"

# ============================================
# 6. CREATE EXAMPLE SEQUENCE (BONUS)
# ============================================
echo "🧬 Creating example TP53 sequence..."

cat > sequences/example_TP53.fasta << 'EOF'
>NM_000546.6 Homo sapiens tumor protein p53 (TP53) - partial coding sequence
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA
GACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATG
GATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCA
GATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCT
ACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAG
AAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAG
TCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACC
TGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATG
GCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAG
CGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAAT
TTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTAT
GAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGT
TCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCC
AGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGA
GACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCC
CCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAG
AAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATG
TTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGG
GGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCAT
AAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGA
EOF

echo "✅ example_TP53.fasta created"

# ============================================
# SUMMARY
# ============================================
echo ""
echo "=========================================="
echo "✅ ALL FILES CREATED SUCCESSFULLY!"
echo "=========================================="
echo ""
echo "Files created:"
echo "  ✅ .gitignore"
echo "  ✅ LICENSE"
echo "  ✅ README.md"
echo "  ✅ sequences/.gitkeep"
echo "  ✅ output/.gitkeep"
echo "  ✅ data/.gitkeep"
echo "  ✅ sequences/example_TP53.fasta"
echo "  ✅ test_install.py (renamed if needed)"
echo ""
echo "Next steps:"
echo ""
echo "1. Review the files (especially README.md)"
echo ""
echo "2. Initialize Git repository:"
echo "   git init"
echo "   git add ."
echo "   git commit -m 'Initial commit: Mutation Impact Analyzer v1.0'"
echo ""
echo "3. Create GitHub repository:"
echo "   Go to: https://github.com/new"
echo "   Repository name: mutation-impact-analyzer"
echo "   DO NOT initialize with README"
echo ""
echo "4. Push to GitHub:"
echo "   git remote add origin https://github.com/$GITHUB_USER/mutation-impact-analyzer.git"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "🎉 Ready to upload to GitHub!"
echo "" 
