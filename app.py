import streamlit as st
import pandas as pd
import plotly.express as px
from src.parser import parse_vcf
from src.annotator import annotate_variants
from src.prioritizer import prioritize_variants
import tempfile
import os

st.set_page_config(page_title="Variant Annotator", layout="wide")

st.title("🧬 Clinical Variant Annotation Pipeline")
st.markdown("Upload a VCF file to annotate variants with clinical databases and prioritize by pathogenicity.")

# file upload
uploaded_file = st.file_uploader("Upload VCF file", type=['vcf'])

if uploaded_file is not None:
    with tempfile.NamedTemporaryFile(delete=False, suffix='.vcf') as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name
    
    try:
        with st.spinner("Parsing VCF file..."):
            df = parse_vcf(tmp_path)
        st.success(f"Parsed {len(df)} variants")
        
        with st.expander("View parsed variants"):
            st.dataframe(df)
        
        if st.button("Run Annotation Pipeline"):
            with st.spinner("Running annotations... This may take a minute."):
                df = annotate_variants(df)
                df = prioritize_variants(df)
                st.session_state['annotated_df'] = df
            st.success("Annotation complete!")
    
    finally:
        os.unlink(tmp_path)

# show results if available
if 'annotated_df' in st.session_state:
    df = st.session_state['annotated_df']
    
    st.header("Results")
    
    # sidebar filters
    st.sidebar.header("Filters")
    
    # priority filter
    priority_options = df['priority_tier'].unique().tolist()
    selected_priorities = st.sidebar.multiselect(
        "Priority Tier",
        priority_options,
        default=priority_options
    )
    
    # gene filter
    genes = df['gene_symbol'].unique().tolist()
    selected_genes = st.sidebar.multiselect(
        "Gene",
        genes,
        default=genes
    )
    
    # clinvar filter
    clinvar_options = df['clinvar_significance'].unique().tolist()
    selected_clinvar = st.sidebar.multiselect(
        "ClinVar Status",
        clinvar_options,
        default=clinvar_options
    )
    
    # frequency filter
    max_af = st.sidebar.slider(
        "Max Population Frequency",
        0.0, 1.0, 1.0, 0.01
    )
    
    # apply filters
    filtered_df = df[
        (df['priority_tier'].isin(selected_priorities)) &
        (df['gene_symbol'].isin(selected_genes)) &
        (df['clinvar_significance'].isin(selected_clinvar)) &
        (df['gnomad_af'] <= max_af)
    ]
    
    # summary stats
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Showing", f"{len(filtered_df)} / {len(df)}")
    col2.metric("Critical", len(filtered_df[filtered_df['priority_tier'] == 'Critical']))
    col3.metric("High", len(filtered_df[filtered_df['priority_tier'] == 'High']))
    col4.metric("In ClinVar", len(filtered_df[filtered_df['clinvar_significance'] == 'found_in_clinvar']))
    
    # charts
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.subheader("Priority Distribution")
        tier_counts = filtered_df['priority_tier'].value_counts()
        fig = px.pie(values=tier_counts.values, names=tier_counts.index,
                     color=tier_counts.index,
                     color_discrete_map={'Critical': 'red', 'High': 'orange', 
                                        'Medium': 'gold', 'Low': 'green'})
        st.plotly_chart(fig, use_container_width=True)
    
    with col_right:
        st.subheader("Variant Types")
        type_counts = filtered_df['variant_type'].value_counts()
        fig2 = px.bar(x=type_counts.index, y=type_counts.values,
                      labels={'x': 'Type', 'y': 'Count'})
        st.plotly_chart(fig2, use_container_width=True)
    
    # results table
    st.subheader("Annotated Variants")
    st.dataframe(filtered_df, use_container_width=True)
    
    # download button
    csv = filtered_df.to_csv(index=False)
    st.download_button("Download Filtered CSV", csv, "annotated_variants.csv", "text/csv")