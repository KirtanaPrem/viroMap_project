import streamlit as st
import pandas as pd
from collections import Counter
from math import exp, log
from Bio import Entrez

# Safe Entrez email
Entrez.email = "viromap@demo.org"

st.set_page_config(page_title="ViroMap", layout="wide")
st.title("ViroMap")

# --- Input: Virus Keyword ---
keyword = st.text_input("Enter virus keyword (e.g., rabies, influenza, SARS-CoV-2):")

strain_options = []
strain_dict = {}

if keyword:
    try:
        search_handle = Entrez.esearch(db="nucleotide", term=keyword + " AND complete genome", retmax=5)
        search_results = Entrez.read(search_handle)
        id_list = search_results["IdList"]

        if id_list:
            summary_handle = Entrez.esummary(db="nucleotide", id=",".join(id_list))
            summaries = Entrez.read(summary_handle)
            for docsum in summaries:
                uid = docsum["Id"]
                title = docsum["Title"]
                label = f"{uid} | {title[:80]}..."
                strain_options.append(label)
                strain_dict[label] = uid
    except Exception as e:
        st.error(f"Error fetching strains: {e}")

selected_strain = st.selectbox("Select a strain", strain_options) if strain_options else None

def fetch_fasta_by_id(ncbi_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
        return handle.read()
    except:
        return ""

fasta = fetch_fasta_by_id(strain_dict[selected_strain]) if selected_strain else ""

# --- Tabs ---
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Epitope Mimicry", "GNN Prediction"])

# --- Tab 0: FASTA ---
with tabs[0]:
    st.header("FASTA Sequence")
    if fasta and fasta.startswith(">"):
        st.code(fasta, language="fasta")
    elif fasta:
        st.warning("No valid FASTA found.")
    else:
        st.info("Search and select a strain to view its sequence.")

# --- Tab 1: Codon Bias ---
with tabs[1]:
    st.header("Codon Bias Metrics")

    if not fasta or not fasta.startswith(">"):
        st.warning("No valid sequence found. Please try a different strain.")
        st.stop()

    def get_codon_list(seq):
        return [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]

    def rare_codon_freq(codons, rare_set={"ATA", "AGG", "AGA", "CGG", "CTA", "TTA"}):
        return round(sum(1 for c in codons if c in rare_set) / len(codons), 4) if codons else 0

    def codon_pair_bias(codons):
        pairs = [codons[i]+codons[i+1] for i in range(len(codons)-1)]
        expected = 1 / (64 * 64)
        counts = Counter(pairs)
        total = sum(counts.values())
        observed = [c / total for c in counts.values()]
        cpb = sum((o - expected) ** 2 / expected for o in observed)
        return round(cpb, 4)

    def volatility(codons):
        amino_changes = 0
        total = 0
        for codon in codons:
            for i in range(3):
                for base in "ACGT":
                    if base != codon[i]:
                        mutant = codon[:i] + base + codon[i+1:]
                        total += 1
                        if mutant != codon:
                            amino_changes += 1
        return round(amino_changes / total, 4) if total else 0

    def simple_cai(seq):
        ref = {
            "GCT": 0.18, "GCC": 0.55, "GCA": 0.15, "GCG": 0.12,
            "CGT": 0.04, "CGC": 0.11, "CGA": 0.11, "CGG": 0.21,
            "AGA": 0.21, "AGG": 0.32
        }
        codons = get_codon_list(seq)
        weights = [ref.get(codon, 0.01) for codon in codons]
        return round(exp(sum(log(w) for w in weights) / len(weights)), 4) if weights else 0

    def enc(codons):
        total = len(codons)
        freqs = Counter(codons)
        n = sum(freqs.values())
        homo = sum(c * (c - 1) for c in freqs.values())
        return round(61 / (1 + homo / (n * (n - 1))), 2) if n > 1 else 0

    clean_seq = "".join(fasta.splitlines()[1:]).replace(" ", "").replace("\n", "").upper()
    codons = get_codon_list(clean_seq)

    metrics = {
        "Codon Adaptation Index (CAI)": simple_cai(clean_seq),
        "Rare Codon Frequency": rare_codon_freq(codons),
        "Codon Pair Bias (CPB)": codon_pair_bias(codons),
        "Codon Volatility": volatility(codons),
        "Effective Number of Codons (ENC)": enc(codons)
    }

    st.dataframe(pd.DataFrame(metrics.items(), columns=["Metric", "Value"]))

# --- Tab 2: LLPS ---
with tabs[2]:
    st.header("LLPS Prediction")

    if not fasta or not fasta.startswith(">"):
        st.warning("No valid sequence found. Please search and select a strain.")
        st.stop()

    def percent_disorder(seq):
        return round(sum(1 for aa in seq if aa in "PESQKRDG") / len(seq) * 100, 2)

    def detect_prion_like(seq):
        return any(motif in seq for motif in ["QQ", "QN", "NQ", "SS", "YG", "RG"])

    def detect_lcr(seq):
        return sum(1 for i in range(len(seq)-6) if len(set(seq[i:i+6])) <= 2)

    def llps_score(seq):
        if not seq: return 0
        score = percent_disorder(seq)/100 + detect_lcr(seq)/50
        if detect_prion_like(seq): score += 0.3
        return round(min(score, 1.0), 2)

    clean_seq = "".join(fasta.splitlines()[1:]).replace(" ", "").replace("\n", "").replace("T", "U").upper()

    data = {
        "LLPS Propensity Score (0–1)": llps_score(clean_seq),
        "Percent Disorder": percent_disorder(clean_seq),
        "Prion-like Motifs Present": "Yes" if detect_prion_like(clean_seq) else "No",
        "Low Complexity Regions": detect_lcr(clean_seq)
    }

    st.dataframe(pd.DataFrame(data.items(), columns=["Metric", "Value"]))

# --- Tab 3: Epitope Mimicry (Placeholder) ---
with tabs[3]:
    st.header("Epitope Mimicry")
    st.info("This module will identify host-virus mimicry (coming soon).")

# --- Tab 4: GNN Prediction (Placeholder) ---
with tabs[4]:
    st.header("GNN Drug Prediction")
    st.info("This module will use Graph Neural Networks to predict antiviral drugs (coming soon).")
