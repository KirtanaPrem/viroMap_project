# ViroMap – Streamlit Bioinformatics App
# Final polished version

import streamlit as st
import pandas as pd
import requests
from io import StringIO
from Bio import Entrez  # ✅ ensure biopython is in requirements.txt

# -----------------------
# Page Config & Theme
# -----------------------
st.set_page_config(page_title="ViroMap", layout="wide")

st.markdown("""
<style>
body {background:white;}
h1, h2 {color:#AA336A;}
.stTabs [data-baseweb="tab"]{
    background:#f8f8f8;border:2px solid #ccc;
    border-radius:8px;padding:10px 15px;margin-right:6px;
    font-weight:600;color:#444;
}
.stTabs [aria-selected="true"]{
    background:#AA336A !important;color:white !important;
}
</style>
""", unsafe_allow_html=True)

# -----------------------
# Header
# -----------------------
st.markdown("<h1 style='text-align:center'>🧬 ViroMap – Integrated Viral Bioinformatics Portal</h1>", unsafe_allow_html=True)

# -----------------------
# Data Files Map
# -----------------------
file_map = {
    "🧬 Codon Bias": "data/sars-cov-2_codon_bias.csv",
    "💧 LLPS":       "data/sars-cov-2_llps.csv",
    "🧫 Mimicry":    "data/sars-cov-2_mimicry.csv"
}

# -----------------------
# Tabs
# -----------------------
tabs = st.tabs([
    "🏠 Overview",
    "🧬 Codon Bias",
    "💧 LLPS",
    "🧫 Mimicry",
    "🧠 GNN Predictions",
    "📤 Upload & Explore",
    "🔍 Fetch Spike FASTA"
])

# -----------------------
# Tab 1 – Overview
# -----------------------
with tabs[0]:
    st.subheader("🏠 Overview")
    st.markdown("""
**ViroMap** integrates multi-layered viral bioinformatics:

| Layer | Description |
|-------|-------------|
| Codon Bias | Viral adaptation to host |
| LLPS | Protein phase separation prediction |
| Mimicry | Viral–host mimicry |
| GNN | Drug prediction for spike |
| Fetch FASTA | Download spike sequence by strain name |
""")

# -----------------------
# Tab 2 – Codon Bias
# -----------------------
with tabs[1]:
    st.subheader("🧬 Codon Bias Analysis")
    df = pd.read_csv(file_map["🧬 Codon Bias"])
    query = st.text_input("🔍 Search Codon Bias", key="codon")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "codon_bias.csv")

# -----------------------
# Tab 3 – LLPS
# -----------------------
with tabs[2]:
    st.subheader("💧 LLPS Prediction")
    df = pd.read_csv(file_map["💧 LLPS"])
    query = st.text_input("🔍 Search LLPS", key="llps")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "llps.csv")

# -----------------------
# Tab 4 – Molecular Mimicry
# -----------------------
with tabs[3]:
    st.subheader("🧫 Mimicry Predictions")
    df = pd.read_csv(file_map["🧫 Mimicry"])
    query = st.text_input("🔍 Search Mimicry", key="mimicry")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "mimicry.csv")

# -----------------------
# Tab 5 – GNN Predictions
# -----------------------
with tabs[4]:
    st.subheader("🧠 GNN Drug Predictions for Spike")
    try:
        gnn_df = pd.read_csv("data/real_gnn_predictions.csv")
        gnn_df["GNN_pKd"] = pd.to_numeric(gnn_df["GNN_pKd"], errors="coerce")

        search = st.text_input("🔍 Search Drug", key="gnn")
        filtered = gnn_df[gnn_df["Drug"].str.contains(search, case=False)] if search.strip() else gnn_df

        filtered = filtered.dropna(subset=["GNN_pKd"])
        if not filtered.empty:
            def highlight(v):
                if pd.isna(v): return ""
                return f"color:{'green' if v>7.2 else 'orange' if v>6.9 else 'red'}"
            st.dataframe(filtered.style.applymap(highlight, subset=["GNN_pKd"]), use_container_width=True)

            if len(filtered) > 1:
                st.subheader("📊 Binding Affinity")
                st.bar_chart(filtered.set_index("Drug")["GNN_pKd"])

            st.download_button("📥 Download Filtered", filtered.to_csv(index=False), "filtered_gnn.csv")
        else:
            st.warning("No matches.")
    except FileNotFoundError:
        st.error("`real_gnn_predictions.csv` not found in /data/.")

# -----------------------
# Tab 6 – Upload & Explore
# -----------------------
with tabs[5]:
    st.subheader("📤 Upload and Explore Custom CSV")
    up = st.file_uploader("Upload CSV", type="csv")
    if up:
        df = pd.read_csv(up)
        q = st.text_input("🔍 Search Upload", key="upload")
        if q.strip():
            df = df[df.apply(lambda r: r.astype(str).str.contains(q, case=False).any(), axis=1)]
        st.dataframe(df, use_container_width=True)
        st.download_button("📥 Download Filtered", df.to_csv(index=False), "upload_filtered.csv")
    else:
        st.info("Upload a CSV to explore.")

# -----------------------
# Tab 7 – Fetch Spike FASTA
# -----------------------
Entrez.email = "your-email@example.com"  # <-- Replace with your actual email

def fetch_spike_fasta(strain):
    try:
        handle = Entrez.esearch(db="nucleotide", term=f"{strain} AND spike", retmax=1)
        record = Entrez.read(handle)
        if not record["IdList"]:
            return None, "No sequence found."
        seq_id = record["IdList"][0]
        fasta = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text").read()
        return fasta, None
    except Exception as e:
        return None, str(e)

with tabs[6]:
    st.subheader("🔍 Fetch Spike FASTA from NCBI")
    strain_input = st.text_input("Enter strain name (e.g., SARS-CoV-2 Wuhan):")
    if strain_input:
        fasta_data, err = fetch_spike_fasta(strain_input)
        if err:
            st.error(err)
        else:
            st.code(fasta_data, language="fasta")
            st.download_button("📥 Download FASTA",
                               fasta_data,
                               f"{strain_input.replace(' ','_')}_spike.fasta")
