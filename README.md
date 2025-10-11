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

```bash
# Clone the repository
git clone https://github.com/Hema291002/mutation-impact-analyzer.git
cd mutation-impact-analyzer

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
python test_install.py
```

## 🎯 Quick Start

```bash
python analyzer.py
```

Follow the interactive prompts:
1. Choose input method (FASTA file or paste sequence)
2. Enter gene name (e.g., TP53, BRCA1, KRAS)
3. Specify nucleotide position to mutate
4. Enter new nucleotide (A/T/C/G)
5. View results and generate visualizations

## 📖 Example Analysis

```
Gene: TP53
Position: 524 (codon 175)
Original: CGC (Arginine)
Mutated: CAC (Histidine)
Type: Missense (R175H)
Pathogenicity: LIKELY PATHOGENIC
Risk: HIGH (5/5)
```

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

```
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
```

## 🔬 How It Works

1. **Sequence Processing** - Load from FASTA or direct input, validate composition
2. **Classification** - Translate DNA to protein, compare amino acids
3. **Database Queries** - Search ClinVar, COSMIC, OMIM for clinical significance
4. **Risk Assessment** - Score based on mutation type and gene importance
5. **Visualization** - Generate interactive Plotly charts and HTML dashboards

## 📊 Output

The tool generates two types of output:

### Text Reports (`output/*_report.txt`)
- Mutation summary
- Molecular impact
- Clinical significance
- Recommendations

### Interactive HTML Dashboard (`output/*_dashboard.html`)
- Summary cards with key metrics
- Interactive sequence viewer
- Codon change visualization
- Pathogenicity gauge (0-5 scale)
- Pathway network diagrams
- Disease association charts

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
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

- **Issues**: [GitHub Issues](https://github.com/Hema291002/mutation-impact-analyzer/issues)
- **Email**: hema98661@gmail.com

## 🌟 Support This Project

If you find this tool useful:
- ⭐ Star this repository
- 🐛 Report bugs and issues
- 💡 Suggest new features
- 🤝 Contribute code

---

**Made with 🧬 and ❤️ for the bioinformatics community**
