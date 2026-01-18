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

# sidebar info
st.sidebar.header("About")
st.sidebar.markdown("""
This tool annotates genetic variants with:
- **ClinVar**: Known disease associations
- **gnomAD**: Population frequency
- **CADD**: Pathogenicity prediction
- **VEP**: Gene and consequence
""")

# file upload
uploaded_file = st.file_uploader("Upload VCF file", type=['vcf'])

# also allow using sample file
use_sample = st.checkbox("Use sample file instead")

if use_sample:
    file_path = 'data/test_variants.vcf'
    if os.path.exists(file_path):
        with st.spinner("Parsing sample VCF file..."):
            df = parse_vcf(file_path)
        st.success(f"Parsed {len(df)} variants from sample file")
        
        with st.expander("View parsed variants"):
            st.dataframe(df)
        
        if st.button("Run Annotation Pipeline", key="sample"):
            progress = st.progress(0)
            status = st.empty()
            
            status.text("Running ClinVar annotation...")
            progress.progress(25)
            
            with st.spinner("Running annotations... This may take a minute."):
                df = annotate_variants(df)
                progress.progress(75)
                df = prioritize_variants(df)
                progress.progress(100)
            
            st.session_state['annotated_df'] = df
            status.text("")
            st.success("Annotation complete!")
    else:
        st.error("Sample file not found")

elif uploaded_file is not None:
    with tempfile.NamedTemporaryFile(delete=False, suffix='.vcf') as tmp:
        tmp.write(uploaded_file.getvalue())
        tmp_path = tmp.name
    
    try:
        with st.spinner("Parsing VCF file..."):
            df = parse_vcf(tmp_path)
        st.success(f"Parsed {len(df)} variants")
        
        with st.expander("View parsed variants"):
            st.dataframe(df)
        
        if st.button("Run Annotation Pipeline", key="upload"):
            progress = st.progress(0)
            status = st.empty()
            
            status.text("Running annotations...")
            progress.progress(25)
            
            with st.spinner("Running annotations... This may take a minute."):
                df = annotate_variants(df)
                progress.progress(75)
                df = prioritize_variants(df)
                progress.progress(100)
            
            st.session_state['annotated_df'] = df
            status.text("")
            st.success("Annotation complete!")
    
    finally:
        os.unlink(tmp_path)

# show results if available
if 'annotated_df' in st.session_state:
    df = st.session_state['annotated_df']
    
    st.header("Results")
    
    # filters in sidebar
    st.sidebar.header("Filters")
    
    priority_options = df['priority_tier'].unique().tolist()
    selected_priorities = st.sidebar.multiselect(
        "Priority Tier",
        priority_options,
        default=priority_options
    )
    
    genes = df['gene_symbol'].unique().tolist()
    selected_genes = st.sidebar.multiselect(
        "Gene",
        genes,
        default=genes
    )
    
    clinvar_options = df['clinvar_significance'].unique().tolist()
    selected_clinvar = st.sidebar.multiselect(
        "ClinVar Status",
        clinvar_options,
        default=clinvar_options
    )
    
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
    st.subheader("Summary")
    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Total", len(df))
    col2.metric("Showing", len(filtered_df))
    col3.metric("Critical", len(filtered_df[filtered_df['priority_tier'] == 'Critical']))
    col4.metric("High", len(filtered_df[filtered_df['priority_tier'] == 'High']))
    col5.metric("In ClinVar", len(filtered_df[filtered_df['clinvar_significance'] == 'found_in_clinvar']))
    
    # charts
    st.subheader("Visualizations")
    col_left, col_mid, col_right = st.columns(3)
    
    with col_left:
        st.markdown("**Priority Distribution**")
        tier_counts = filtered_df['priority_tier'].value_counts()
        fig = px.pie(values=tier_counts.values, names=tier_counts.index,
                     color=tier_counts.index,
                     color_discrete_map={'Critical': '#e74c3c', 'High': '#e67e22', 
                                        'Medium': '#f1c40f', 'Low': '#27ae60'})
        fig.update_layout(margin=dict(t=0, b=0, l=0, r=0))
        st.plotly_chart(fig, use_container_width=True)
    
    with col_mid:
        st.markdown("**Variant Types**")
        type_counts = filtered_df['variant_type'].value_counts()
        fig2 = px.bar(x=type_counts.index, y=type_counts.values,
                      labels={'x': 'Type', 'y': 'Count'})
        fig2.update_layout(margin=dict(t=0, b=0, l=0, r=0))
        st.plotly_chart(fig2, use_container_width=True)
    
    with col_right:
        st.markdown("**Consequence Types**")
        cons_counts = filtered_df['consequence'].value_counts()
        fig3 = px.bar(x=cons_counts.values, y=cons_counts.index, orientation='h',
                      labels={'x': 'Count', 'y': 'Consequence'})
        fig3.update_layout(margin=dict(t=0, b=0, l=0, r=0))
        st.plotly_chart(fig3, use_container_width=True)
    
    # results table
    st.subheader("Annotated Variants")
    
    # highlight critical/high rows
    def highlight_priority(row):
        if row['priority_tier'] == 'Critical':
            return ['background-color: #ffcccc'] * len(row)
        elif row['priority_tier'] == 'High':
            return ['background-color: #f1c232'] * len(row)
        return [''] * len(row)
    
    styled_df = filtered_df.style.apply(highlight_priority, axis=1)
    st.dataframe(styled_df, use_container_width=True)
    
    # download buttons
    st.subheader("Export")
    col_dl1, col_dl2 = st.columns(2)
    
    with col_dl1:
        csv = filtered_df.to_csv(index=False)
        st.download_button("Download Filtered CSV", csv, "filtered_variants.csv", "text/csv")
    
    with col_dl2:
        full_csv = df.to_csv(index=False)
        st.download_button("Download All CSV", full_csv, "all_variants.csv", "text/csv")