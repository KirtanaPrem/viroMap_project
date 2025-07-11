import streamlit as st
import pandas as pd
from Bio import Entrez
from functools import lru_cache

# -------------------
# App Setup and Style
# -------------------
st.set_page_config(page_title="ViroMap Unified", layout="wide")

st.markdown("""
<style>
body {background:white;}
h1 {color:#AA336A;}
.stTabs [data-baseweb="tab"] {
    background:#f8f8f8; border:2px solid #ccc;
    padding:10px 15px; margin-right:6px;
    border-radius:8px; font-weight:600; color:#444;
}
.stTabs [aria-selected="true"] {
    background:#AA336A !important; color:white !important;
}
</style>
""", unsafe_allow_html=True)

# -------------------
# App Title and Search
# -------------------
st.markdown("<h1 style='text-align:center'>🧬 ViroMap – Unified Search</h1>", unsafe_allow_html=True)
strain = st.text_input("🔍 Type virus name (e.g. SARS-CoV-2 Omicron):").strip()

# -------------------
# Email for NCBI fetch
# -------------------
Entrez.email = "your-email@example.com"  # Replace with your real email

# -------------------
# Functions for fetching data
# -------------------
@lru_cache(maxsize=10)
def fetch_fasta(strain_name):
    try:
        handle = Entrez.esearch(db="nucleotide", term=f"{strain_name} AND spike", retmax=1)
        record = Entrez.read(handle)
        if not record["IdList"]:
            return None
        seq_id = record["IdList"][0]
        fasta = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text").read()
        return fasta
    except Exception as e:
        return None

def get_demo_data(name):
    # This is where your demo CSVs go
    path = f"data/{name}_demo.csv"
    try:
        return pd.read_csv(path)
    except:
        return pd.DataFrame()

# -------------------
# Tabs for the sections
# -------------------
tab_names = ["🧬 Codon Bias", "💧 LLPS", "🧫 Mimicry", "🧠 GNN Predictions"]
tabs = st.tabs(tab_names)

# -------------------
# If user typed a virus
# -------------------
if strain:
    fasta = fetch_fasta(strain)
    if fasta:
        # You can display FASTA if you want
        st.success("Spike protein FASTA found!")
        st.code(fasta, language="fasta")

        with tabs[0]:
            st.subheader("🧬 Codon Bias")
            df = get_demo_data("codon")
            st.dataframe(df, use_container_width=True)

        with tabs[1]:
            st.subheader("💧 LLPS Prediction")
            df = get_demo_data("llps")
            st.dataframe(df, use_container_width=True)

        with tabs[2]:
            st.subheader("🧫 Mimicry Summary")
            df = get_demo_data("mimicry")
            st.dataframe(df, use_container_width=True)

        with tabs[3]:
            st.subheader("🧠 GNN Drug Predictions")
            df = get_demo_data("gnn")
            st.dataframe(df, use_container_width=True)
    else:
        st.error("❌ Could not fetch spike protein for that strain. Try another.")
else:
    for t in tabs:
        with t:
            st.info("🧠 Type a virus/strain name above to view results.")
