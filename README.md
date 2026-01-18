# Clinical Variant Annotation Pipeline

A tool that annotates genetic variants with clinical databases and prioritizes them by pathogenicity.

## What It Does

When patients get DNA sequencing, the results contain thousands of genetic variants. Most are harmless, but some cause disease. This tool automates finding the important ones by:

- Checking each variant against clinical databases
- Scoring variants by how likely they are to cause disease
- Ranking and filtering results

## Annotations

- **ClinVar**: Known disease associations from NCBI
- **gnomAD**: Population frequency (rare variants are more suspicious)
- **CADD**: Computational pathogenicity prediction
- **VEP**: Gene name and variant consequence (missense, frameshift, etc.)

## Priority Scoring

Variants are scored based on:
- ClinVar status (known pathogenic = higher score)
- Population frequency (rare = higher score)
- CADD score (predicted damaging = higher score)
- Consequence type (frameshift/nonsense = higher score)

Tiers: Critical (80+), High (50-79), Medium (30-49), Low (<30)

## Setup
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

Run the Streamlit app:
```bash
streamlit run app.py
```

Then:
1. Upload a VCF file (or use sample file)
2. Click "Run Annotation Pipeline"
3. Filter and explore results
4. Download CSV

## Project Structure
```
clinical-variant-annotation/
├── app.py                 # Streamlit dashboard
├── src/
│   ├── parser.py          # VCF file parsing
│   ├── annotator.py       # Database annotations
│   └── prioritizer.py     # Scoring algorithm
├── data/
│   └── test_variants.vcf  # Sample file
└── requirements.txt
```

## Technologies

- Python
- Streamlit
- Pandas
- Plotly
- REST APIs (ClinVar, gnomAD, VEP)