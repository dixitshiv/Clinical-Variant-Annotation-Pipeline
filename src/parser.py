import pandas as pd

def create_variant_id(chrom, pos, ref, alt):
    """Create unique variant identifier."""
    return f"{chrom}-{pos}-{ref}-{alt}"

def parse_vcf(filepath):
    """Read VCF file and return dataframe with variant info."""
    
    variants = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # skip header lines
            if line.startswith('##'):
                continue
            # column names line
            if line.startswith('#CHROM'):
                continue
            
            fields = line.strip().split('\t')
            
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            
            # determine variant type
            if len(ref) == 1 and len(alt) == 1:
                var_type = 'SNV'
            elif len(ref) > len(alt):
                var_type = 'deletion'
            else:
                var_type = 'insertion'
            
            var_id = create_variant_id(chrom, pos, ref, alt)
            
            variants.append({
                'variant_id': var_id,
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filt,
                'variant_type': var_type
            })
    
    return pd.DataFrame(variants)


if __name__ == '__main__':
    df = parse_vcf('data/sample.vcf')
    print(df)