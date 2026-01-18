import pandas as pd

def calculate_priority_score(row):
    """Calculate priority score for a variant. Higher = more important."""
    
    score = 0
    
    # ClinVar significance
    clinvar = row.get('clinvar_significance', '')
    if clinvar == 'found_in_clinvar':
        score += 30
    
    # Population frequency (rare = higher priority)
    af = row.get('gnomad_af', 0)
    if af == 0:
        score += 30  # not seen in population
    elif af < 0.001:
        score += 25  # very rare
    elif af < 0.01:
        score += 15  # rare
    elif af < 0.05:
        score += 5   # uncommon
    # common variants get 0
    
    # CADD score
    cadd = row.get('cadd_score', -1)
    if cadd >= 30:
        score += 30  # very damaging
    elif cadd >= 20:
        score += 20  # damaging
    elif cadd >= 10:
        score += 10  # possibly damaging
    
    # Variant consequence
    consequence = row.get('consequence', '')
    high_impact = ['frameshift_variant', 'stop_gained', 'splice_donor_variant', 
                   'splice_acceptor_variant', 'stop_lost', 'start_lost']
    moderate_impact = ['missense_variant', 'inframe_deletion', 'inframe_insertion']
    
    if consequence in high_impact:
        score += 30
    elif consequence in moderate_impact:
        score += 20
    elif 'upstream' in consequence or 'downstream' in consequence:
        score += 5
    
    return score


def assign_priority_tier(score):
    """Assign priority tier based on score."""
    if score >= 80:
        return 'Critical'
    elif score >= 50:
        return 'High'
    elif score >= 30:
        return 'Medium'
    else:
        return 'Low'


def prioritize_variants(df):
    """Add priority scores and tiers to variant dataframe."""
    
    df['priority_score'] = df.apply(calculate_priority_score, axis=1)
    df['priority_tier'] = df['priority_score'].apply(assign_priority_tier)
    
    # sort by priority score descending
    df = df.sort_values('priority_score', ascending=False).reset_index(drop=True)
    
    return df


if __name__ == '__main__':
    # test with annotated data
    df = pd.read_csv('data/annotated_variants.csv')
    df = prioritize_variants(df)
    print(df[['variant_id', 'gene_symbol', 'consequence', 'priority_score', 'priority_tier']])