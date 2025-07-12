import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt

# Setup
st.set_page_config(layout="wide", page_title="ViroMap", page_icon="🧬")
st.markdown(
    """
    <style>
        body { background-color: white; }
        section[data-testid="stSidebar"] { background-color: #f5f5f5; }
        .stTabs [data-baseweb="tab"] {
            width: 16%;
            background-color: #e3f2fd;
            padding: 10px;
            border-radius: 8px;
            margin-right: 5px;
            text-align: center;
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
query = st.text_input("🔍 Enter virus name (e.g. COVID-19, HIV, Influenza):")

# Step 1: Search for matching strains from NCBI
def search_strains(query):
    Entrez.email = "example@email.com"
    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)
        record = Entrez.read(handle)
        ids = record["IdList"]
        strains = []
        for id_ in ids:
            fetch = Entrez.efetch(db="nucleotide", id=id_, rettype="gb", retmode="text")
            seq_record = SeqIO.read(fetch, "genbank")
            strains.append((id_, seq_record.description))
        return strains
    except:
        return []

# Step 2: Fetch FASTA
def fetch_fasta_from_id(ncbi_id):
    Entrez.email = "example@email.com"
    fetch = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(fetch, "fasta")
    return str(seq_record.seq), seq_record.id

# Demos
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

# === Workflow ===
selected_strain = None
sequence = ""
record_id = ""

if query:
    strains = search_strains(query)
    if strains:
        options = {desc: id_ for id_, desc in strains}
        choice = st.selectbox("🧬 Select a matching strain:", list(options.keys()))
        selected_strain = options[choice]
        sequence, record_id = fetch_fasta_from_id(selected_strain)
    else:
        st.warning("No matching strains found.")

# Tabs
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Mimicry", "Epitope", "GNN"])

if sequence:
    with tabs[0]:
        st.subheader("📄 FASTA Sequence")
        st.code(sequence, language="fasta")
        st.success(f"Record ID: {record_id}")

    with tabs[1]:
        st.subheader("🧬 Codon Bias")
        df = codon_bias(sequence)
        st.dataframe(df)
        st.bar_chart(df.set_index("Codon"))

    with tabs[2]:
        st.subheader("🔬 LLPS Prediction")
        st.dataframe(llps_demo(sequence))

    with tabs[3]:
        st.subheader("🎭 Molecular Mimicry")
        st.dataframe(mimicry_demo(sequence))

    with tabs[4]:
        st.subheader("🎯 Epitope Prediction")
        st.dataframe(epitope_demo(sequence))

    with tabs[5]:
        st.subheader("🧠 GNN Drug Repurposing (Demo)")
        st.caption("This is a demo. Real GNN model is optional later.")
        st.dataframe(gnn_demo())
else:
    st.info("Start by entering a virus name and selecting a strain.")
