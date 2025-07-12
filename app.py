import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
import requests
import matplotlib.pyplot as plt

# Streamlit config
st.set_page_config(layout="wide", page_title="ViroMap", page_icon="🧬")

st.markdown(
    """
    <style>
        body {
            background-color: white;
        }
        .stTabs [data-baseweb="tab"] {
            background-color: #e3f2fd;
            padding: 10px;
            border-radius: 8px;
            margin: 2px;
        }
        .stTabs [aria-selected="true"] {
            background-color: #90caf9;
            font-weight: bold;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 ViroMap: Unified Viral Prediction Platform")
st.markdown("Search any virus or strain to fetch FASTA and view predictions.")

# Search bar
query = st.text_input("🔍 Enter virus or strain name (e.g., SARS-CoV-2, HIV-1, Dengue):")

# Initialize placeholders
fasta = ""
record_id = ""

# Fetch sequence from NCBI
def fetch_fasta(query):
    Entrez.email = "your@email.com"
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        if not id_list:
            return "", "No sequence found."
        fetch_handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="fasta", retmode="text")
        seq_record = SeqIO.read(fetch_handle, "fasta")
        return str(seq_record.seq), seq_record.id
    except Exception as e:
        return "", f"Error: {e}"

# Codon Bias Dummy
def calculate_codon_bias(seq):
    from collections import Counter
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    codon_freq = Counter(codons)
    return pd.DataFrame(codon_freq.items(), columns=["Codon", "Frequency"]).sort_values(by="Frequency", ascending=False)

# LLPS Mock
def predict_llps(seq):
    return pd.DataFrame({
        "Region": ["N-terminal", "Middle", "C-terminal"],
        "LLPS Propensity": [0.8, 0.4, 0.6]
    })

# Mimicry Demo
def predict_mimicry(seq):
    return pd.DataFrame({
        "Host Protein": ["ACE2", "HLA-A", "TMPRSS2"],
        "Mimicry Score": [0.87, 0.76, 0.65]
    })

# Epitope Prediction Demo
def predict_epitopes(seq):
    return pd.DataFrame({
        "Epitope": ["SYGFQPT", "YGFQPTY", "GFQPTNG"],
        "MHC Binding Affinity": [85, 70, 60]
    })

# GNN Drug Prediction Demo
def gnn_predict_demo():
    return pd.DataFrame({
        "Drug": ["Remdesivir", "Molnupiravir", "Favipiravir"],
        "Predicted Binding Score (pKd)": [7.1, 6.5, 5.9]
    })

# Create tabs
tabs = st.tabs(["FASTA Sequence", "Codon Bias", "LLPS Prediction", "Mimicry", "Epitope Prediction", "GNN Drug Prediction"])

# Process if query exists
if query:
    fasta, record_id = fetch_fasta(query)

    with tabs[0]:
        st.subheader("📄 FASTA Sequence")
        if fasta:
            st.code(fasta, language="fasta")
            st.success(f"Record ID: {record_id}")
        else:
            st.error("FASTA not found.")

    with tabs[1]:
        st.subheader("🧬 Codon Bias Analysis")
        if fasta:
            codon_df = calculate_codon_bias(fasta)
            st.dataframe(codon_df, use_container_width=True)
            st.bar_chart(codon_df.set_index("Codon"))

    with tabs[2]:
        st.subheader("🔬 LLPS Propensity Prediction")
        if fasta:
            llps_df = predict_llps(fasta)
            st.dataframe(llps_df, use_container_width=True)

    with tabs[3]:
        st.subheader("🎭 Molecular Mimicry Prediction (Demo)")
        if fasta:
            mimic_df = predict_mimicry(fasta)
            st.dataframe(mimic_df, use_container_width=True)

    with tabs[4]:
        st.subheader("🎯 Epitope Binding Prediction (NetMHC-like)")
        if fasta:
            epitope_df = predict_epitopes(fasta)
            st.dataframe(epitope_df, use_container_width=True)

    with tabs[5]:
        st.subheader("🧠 GNN-Based Drug Binding Prediction (Demo Only)")
        st.caption("Real GNN temporarily disabled due to cloud install limits.")
        gnn_df = gnn_predict_demo()
        st.dataframe(gnn_df, use_container_width=True)

else:
    st.warning("Please enter a virus name above to begin.")
