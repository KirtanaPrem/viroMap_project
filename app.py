import streamlit as st
import pandas as pd
from collections import Counter
from functools import lru_cache
from Bio import Entrez

# ============  Page & Style  ============
st.set_page_config(page_title="ViroMap Unified", layout="wide")

st.markdown("""
<style>
body {background:white;}
h1 {color:#AA336A;}
.stTabs [data-baseweb="tab"]{
    background:#f2f2f2; border:2px solid #ccc;
    padding:10px 15px; margin-right:6px;
    border-radius:8px; font-weight:600; color:#333;
}
.stTabs [aria-selected="true"]{
    background:#AA336A !important; color:white !important;
}
input[type="text"]{
    height:36px !important; font-size:14px !important;
}
</style>
""", unsafe_allow_html=True)

# ============  Header  ============
st.markdown("<h1 style='text-align:center'>🧬 ViroMap – Unified Viral Dashboard</h1>", unsafe_allow_html=True)

# ============  Search Bar  ============
strain = st.text_input("🔍  Type virus / strain (e.g. SARS‑CoV‑2 Omicron):").strip()

# ============  NCBI e‑mail ============
Entrez.email = "your-email@example.com"          # <- put your real e‑mail

# ============  Helper Functions ============
@lru_cache(maxsize=10)
def fetch_fasta(strain_name):
    """Get spike FASTA for a strain from NCBI."""
    try:
        h = Entrez.esearch(db="nucleotide", term=f"{strain_name} AND spike", retmax=1)
        rec = Entrez.read(h)
        if not rec["IdList"]:
            return None
        fid = rec["IdList"][0]
        fasta = Entrez.efetch(db="nucleotide", id=fid, rettype="fasta", retmode="text").read()
        return fasta
    except Exception:
        return None

def calculate_codon_usage(fasta_text):
    """Real‑time codon count + frequency from FASTA."""
    seq = "".join(fasta_text.splitlines()[1:]).upper()  # drop header, join
    seq = seq[: len(seq) - len(seq) % 3]                # make length /3
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    counts = Counter(codons)
    total = sum(counts.values())
    data = [
        {"Codon": c, "Count": n, "Frequency": round(n/total, 4)}
        for c, n in sorted(counts.items())
    ]
    return pd.DataFrame(data)

def get_demo_data(name):
    """Load demo CSVs for LLPS, Mimicry, GNN."""
    try:
        return pd.read_csv(f"data/{name}_demo.csv")
    except FileNotFoundError:
        return pd.DataFrame()

# ============  Tabs  ============
tab_names = ["📄 FASTA Sequence", "🧬 Codon Bias", "💧 LLPS", "🧫 Mimicry", "🧠 GNN Predictions"]
tabs = st.tabs(tab_names)

# ============  Logic  ============
if strain:
    fasta = fetch_fasta(strain)
    if fasta:
        # ---- FASTA Tab ----
        with tabs[0]:
            st.subheader("📄 Spike Protein FASTA")
            st.code(fasta, language="fasta")

        # ---- Codon Bias Tab ----
        with tabs[1]:
            st.subheader("🧬 Real‑time Codon Usage")
            codon_df = calculate_codon_usage(fasta)
            st.dataframe(codon_df, use_container_width=True)

        # ---- LLPS Tab ---- (demo csv)
        with tabs[2]:
            st.subheader("💧 LLPS Prediction (demo)")
            llps_df = get_demo_data("llps")
            st.dataframe(llps_df, use_container_width=True)

        # ---- Mimicry Tab ---- (demo csv)
        with tabs[3]:
            st.subheader("🧫 Mimicry Summary (demo)")
            mim_df = get_demo_data("mimicry")
            st.dataframe(mim_df, use_container_width=True)

        # ---- GNN Predictions Tab ---- (demo csv)
        with tabs[4]:
            st.subheader("🧠 GNN Drug Predictions (demo)")
            gnn_df = get_demo_data("gnn")
            st.dataframe(gnn_df, use_container_width=True)
    else:
        for t in tabs:
            with t:
                st.error("❌ Could not fetch spike for that strain. Try another.")
else:
    for t in tabs:
        with t:
            st.info("🧠 Enter a virus/strain name above to load data.")
