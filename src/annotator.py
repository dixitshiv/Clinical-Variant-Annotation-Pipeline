import requests
import pandas as pd

def annotate_with_clinvar(df):
    """Add ClinVar annotations using NCBI E-utilities API."""
    
    clinical_significance = []
    conditions = []
    
    for _, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        # search ClinVar for this variant
        search_term = f"{chrom}[chr] AND {pos}[chrpos37]"
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "clinvar",
            "term": search_term,
            "retmode": "json"
        }
        
        try:
            response = requests.get(search_url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                id_list = data.get('esearchresult', {}).get('idlist', [])
                
                if id_list:
                    clin_sig = "found_in_clinvar"
                    condition = "see_clinvar"
                else:
                    clin_sig = "not_in_clinvar"
                    condition = "none"
            else:
                clin_sig = "api_error"
                condition = "api_error"
        except Exception as e:
            clin_sig = "error"
            condition = "error"
        
        clinical_significance.append(clin_sig)
        conditions.append(condition)
        print(f"ClinVar: {row['variant_id']} -> {clin_sig}")
    
    df['clinvar_significance'] = clinical_significance
    df['clinvar_condition'] = conditions
    
    return df

def annotate_with_gnomad(df):
    """Add population frequency from gnomAD."""
    
    frequencies = []
    
    for _, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        url = f"https://gnomad.broadinstitute.org/api"
        query = """
        query($variantId: String!) {
            variant(variantId: $variantId, dataset: gnomad_r2_1) {
                genome {
                    ac
                    an
                    af
                }
            }
        }
        """
        
        variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        
        try:
            response = requests.post(
                url,
                json={"query": query, "variables": {"variantId": variant_id}},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                variant_data = data.get('data', {}).get('variant')
                if variant_data and variant_data.get('genome'):
                    af = variant_data['genome'].get('af', 0)
                    freq = af if af else 0
                else:
                    freq = 0
            else:
                freq = -1  # api error
        except:
            freq = -1
        
        frequencies.append(freq)
        print(f"gnomAD: {row['variant_id']} -> AF={freq}")
    
    df['gnomad_af'] = frequencies
    
    return df

def annotate_with_cadd(df):
    """Add CADD pathogenicity scores."""
    
    cadd_scores = []
    
    for _, row in df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        url = f"https://cadd.gs.washington.edu/api/v1.0/{chrom}:{pos}_{ref}_{alt}"
        
        try:
            response = requests.get(url, timeout=15)
            
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    # last column is PHRED score
                    fields = lines[1].split('\t')
                    score = float(fields[-1])
                else:
                    score = -1
            else:
                score = -1
        except:
            score = -1
        
        cadd_scores.append(score)
        print(f"CADD: {row['variant_id']} -> score={score}")
    
    df['cadd_score'] = cadd_scores
    
    return df

def annotate_with_vep(df):
    """Add VEP annotations using Ensembl REST API."""
    
    consequences = []
    gene_symbols = []
    
    for _, row in df.iterrows():
        chrom = row['chrom'].replace('chr', '')
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        
        url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}:{pos}/{alt}"
        headers = {"Content-Type": "application/json"}
        
        try:
            response = requests.get(url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    transcript = data[0].get('transcript_consequences', [{}])[0]
                    consequence = transcript.get('consequence_terms', ['unknown'])[0]
                    gene = transcript.get('gene_symbol', 'unknown')
                else:
                    consequence = 'unknown'
                    gene = 'unknown'
            else:
                consequence = 'api_error'
                gene = 'api_error'
        except Exception as e:
            consequence = 'error'
            gene = 'error'
        
        consequences.append(consequence)
        gene_symbols.append(gene)
        print(f"VEP: {row['variant_id']} -> {consequence}, {gene}")
    
    df['consequence'] = consequences
    df['gene_symbol'] = gene_symbols
    
    return df

def annotate_variants(df):
    """Run all annotations on variant dataframe."""
    df = annotate_with_clinvar(df)
    df = annotate_with_gnomad(df)
    df = annotate_with_cadd(df)
    df = annotate_with_vep(df)
    return df