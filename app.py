pip install streamlit biopython pandas numpy
import streamlit as st
import pandas as pd
import numpy as np
from collections import Counter
from math import exp, log
from Bio import Entrez, SeqIO

# Setup Entrez (use your email to avoid NCBI blocking)
Entrez.email = "your_email@example.com"  # Replace with your email

# Page config
st.set_page_config(page_title="ViroMap", layout="wide")
st.title("🧬 ViroMap: Unified Viral Data Analysis Platform")

# Virus name input
virus_name = st.text_input("🔍 Enter virus name (e.g., SARS-CoV-2, HIV, Zika virus):")

# Function to fetch FASTA from NCBI using Entrez
def fetch_fasta_by_virus_name(name):
    try:
        search = Entrez.esearch(db="nucleotide", term=f"{name}[ORGN] AND complete genome", retmax=1)
        result = Entrez.read(search)
        ids = result["IdList"]
        if ids:
            handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
            fasta_data = handle.read()
            return fasta_data
        else:
            return "No sequence found for this virus."
    except Exception as e:
        return f"Error fetching sequence: {e}"

# Fetch the FASTA dynamically
fasta = fetch_fasta_by_virus_name(virus_name) if virus_name else ""

# Tabs
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Epitope Mimicry", "GNN Prediction"])

# Tab 0: FASTA
with tabs[0]:
    st.header("📄 FASTA Sequence")
    if fasta and fasta.startswith(">"):
        st.code(fasta, language="fasta")
    elif fasta:
        st.warning(fasta)
    else:
        st.info("Search for a virus to display its FASTA sequence.")

# Tab 1: Codon Bias
with tabs[1]:
    st.header("🧬 Codon Bias Metrics")

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

    if fasta.startswith(">"):
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
    elif virus_name:
        st.warning("FASTA not available or too short to compute metrics.")
    else:
        st.info("Enter a virus name above to begin analysis.")

# Tabs 2–4: Placeholders
with tabs[2]:
    st.header("💧 LLPS Prediction (Coming Soon)")
    st.info("This module will predict liquid-liquid phase separation using sequence features.")

with tabs[3]:
    st.header("🎯 Epitope Mimicry (Coming Soon)")
    st.info("This module will detect epitope mimicry with host proteins.")

with tabs[4]:
    st.header("🧠 GNN Drug Prediction (Coming Soon)")
    st.info("This will use Graph Neural Networks to predict antiviral drug interactions.")

# Footer
st.markdown("---")
st.caption("Created with ❤️ by [Your Name] | ViroMap Project | Bioinformatics Streamlit App")
