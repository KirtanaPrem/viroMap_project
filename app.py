import streamlit as st
import pandas as pd
from collections import Counter
from functools import lru_cache
from Bio import Entrez
import math

# ============  Page Config and Theme  ============
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

# ============  NCBI Setup ============
Entrez.email = "your-email@example.com"  # Replace with your email

# ============  Fetch FASTA ============
@lru_cache(maxsize=10)
def fetch_fasta(strain_name):
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

# ============  Codon Bias Calculator ============
def calculate_codon_usage(fasta_text):
    seq = "".join(fasta_text.splitlines()[1:]).upper()
    seq = seq[: len(seq) - len(seq) % 3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    counts = Counter(codons)
    total = sum(counts.values())
    data = [
        {"Codon": c, "Count": n, "Frequency": round(n/total, 4)}
        for c, n in sorted(counts.items())
    ]
    return pd.DataFrame(data)

# ============  LLPS Predictor ============
def predict_llps_simple(fasta_seq):
    sequence = "".join(fasta_seq.splitlines()[1:]).upper()

    hydrophobicity = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

    scores = []
    for i, aa in enumerate(sequence):
        hydro = hydrophobicity.get(aa, 0.0)
        disorder = 1 if aa in 'PG' else 0
        score = hydro + disorder * 1.5
        scores.append((i+1, aa, round(score, 2)))

    df = pd.DataFrame(scores, columns=["Position", "AA", "LLPS_Score"])
    df["Region"] = pd.cut(df["LLPS_Score"], bins=[-math.inf, 0, 2, math.inf], labels=["Low", "Medium", "High"])
    return df

# ============  Placeholder for Other Tabs ============
def get_demo_data(name):
    try:
        return pd.read_csv(f"data/{name}_demo.csv")
    except FileNotFoundError:
        return pd.DataFrame()

# ============  Tabs Setup ============
tab_names = ["📄 FASTA Sequence", "🧬 Codon Bias", "💧 LLPS", "🧫 Mimicry", "🧠 GNN Predictions"]
tabs = st.tabs(tab_names)

# ============  Main Logic ============
if strain:
    fasta = fetch_fasta(strain)
    if fasta:
        # --- FASTA Tab ---
        with tabs[0]:
            st.subheader("📄 Spike Protein FASTA")
            st.code(fasta, language="fasta")

        # --- Codon Bias Tab ---
        with tabs[1]:
            st.subheader("🧬 Real‑time Codon Usage")
            codon_df = calculate_codon_usage(fasta)
            st.dataframe(codon_df, use_container_width=True)

        # --- LLPS Prediction Tab ---
        with tabs[2]:
            st.subheader("💧 LLPS Prediction")
            llps_df = predict_llps_simple(fasta)
            st.dataframe(llps_df, use_container_width=True)

        # --- Mimicry Tab (demo for now) ---
        with tabs[3]:
            st.subheader("🧫 Mimicry Prediction (demo)")
            mim_df = get_demo_data("mimicry")
            st.dataframe(mim_df, use_container_width=True)

        # --- GNN Tab (demo for now) ---
        with tabs[4]:
            st.subheader("🧠 GNN Drug Prediction (demo)")
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
