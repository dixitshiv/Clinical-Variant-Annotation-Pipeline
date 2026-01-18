from src.parser import parse_vcf
from src.annotator import annotate_variants

df = parse_vcf('data/test_variants.vcf')
print("Parsed variants:")
print(df)

print("\n--- Running full annotation pipeline ---\n")
df = annotate_variants(df)

print("\n--- Final Results ---\n")
print(df.to_string())

# save to CSV
df.to_csv('data/annotated_variants.csv', index=False)
print("\nSaved to data/annotated_variants.csv")