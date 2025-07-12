import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO

# Set up page
st.set_page_config(layout="wide", page_title="ViroMap", page_icon="🧬")
st.markdown(
    """
    <style>
        body { background-color: white; }
        .stTabs [data-baseweb="tab"] {
            background-color: #e3f2fd;
            border-radius: 8px;
        }
        .stTabs [aria-selected="true"] {
            background-color: #90caf9;
            font-weight: bold;
        }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🧬 ViroMap: Unified Viral Prediction Dashboard")
query = st.text_input("🔍 Search Virus Name (e.g. SARS-CoV-2, HIV, Influenza):")

# Helper: Fetch sequence
def fetch_fasta(query):
    Entrez.email = "example@email.com"
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        record = Entrez.read(handle)
        if not record["IdList"]:
            return "", ""
        fetch_handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta", retmode="text")
        seq_record = SeqIO.read(fetch_handle, "fasta")
        return str(seq_record.seq), seq_record.id
    except:
        return "", ""

# Prediction Demos
def codon_bias(seq):
    from collections import Counter
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    freq = Counter(codons)
    return pd.DataFrame(freq.items(), columns=["Codon", "Frequency"]).sort_values(by="Frequency", ascending=False)

def llps_demo(seq):
    return pd.DataFrame({
        "Region": ["N-term", "Mid", "C-term"],
        "LLPS Score": [0.82, 0.54, 0.70]
    })

def mimicry_demo(seq):
    return pd.DataFrame({
        "Mimic Protein": ["ACE2", "TMPRSS2", "HLA-B"],
        "Similarity Score": [0.88, 0.67, 0.73]
    })

def epitope_demo(seq):
    return pd.DataFrame({
        "Epitope": ["SYGFQPT", "YGFQPTY", "GFQPTNG"],
        "Binding Affinity": [84, 71, 62]
    })

def gnn_demo():
    return pd.DataFrame({
        "Drug": ["Remdesivir", "Molnupiravir", "Favipiravir"],
        "Predicted pKd": [7.2, 6.8, 6.1]
    })

# Tabs
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Mimicry", "Epitope", "GNN"])

if query:
    fasta, rid = fetch_fasta(query)

    with tabs[0]:
        st.subheader("📄 FASTA Sequence")
        if fasta:
            st.code(fasta, language="fasta")
            st.success(f"Fetched sequence ID: {rid}")
        else:
            st.error("No sequence found.")

    with tabs[1]:
        st.subheader("🧬 Codon Bias")
        if fasta:
            st.dataframe(codon_bias(fasta))

    with tabs[2]:
        st.subheader("🔬 LLPS Prediction")
        if fasta:
            st.dataframe(llps_demo(fasta))

    with tabs[3]:
        st.subheader("🎭 Molecular Mimicry")
        if fasta:
            st.dataframe(mimicry_demo(fasta))

    with tabs[4]:
        st.subheader("🎯 Epitope Prediction")
        if fasta:
            st.dataframe(epitope_demo(fasta))

    with tabs[5]:
        st.subheader("🧠 GNN Drug Repurposing (Demo)")
        st.caption("This is a demo until real GNN is enabled.")
        st.dataframe(gnn_demo())
else:
    st.info("Please enter a virus name to begin.")
