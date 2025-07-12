import streamlit as st
import pandas as pd
import numpy as np
from collections import Counter
from math import exp, log

# Setup UI
st.set_page_config(page_title="ViroMap", layout="wide")
st.title("🧬 ViroMap: Unified Viral Data Analysis Platform")

# Virus name input
virus_name = st.text_input("🔍 Enter virus name (e.g., SARS-CoV-2, HIV):")

# Dummy dropdown for now
strain = st.selectbox("Select a strain (simulated)", ["NC_045512.2 (Wuhan-Hu-1)", "MT163716.1 (Indian strain)"])

# Simulated FASTA output (normally you'd fetch this using Entrez)
fasta = """>spike_protein
ATGTTTGTTTTTCTTGTTTTA...GTTTGTAGTTAGTGCAACT
""" if virus_name else ""

# Tabs
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Epitope Mimicry", "GNN Prediction"])

# Tab 0: FASTA Display
with tabs[0]:
    st.header("📄 FASTA Sequence")
    if fasta:
        st.code(fasta, language="fasta")
    else:
        st.info("Search for a virus to display its FASTA sequence.")

# Tab 1: Codon Bias
with tabs[1]:
    st.header("🧬 Codon Bias Metrics")

    def get_codon_list(seq):
        return [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

    def rare_codon_freq(codons, rare_set={"ATA", "AGG", "AGA", "CGG", "CTA", "TTA"}):
        return round(sum(1 for c in codons if c in rare_set) / len(codons), 4)

    def codon_pair_bias(codons):
        pairs = [codons[i]+codons[i+1] for i in range(len(codons)-1)]
        expected = 1 / (64*64)
        counts = Counter(pairs)
        total = sum(counts.values())
        observed = [c/total for c in counts.values()]
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

    if fasta:
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
    else:
        st.warning("Please search and select a virus to view codon bias metrics.")

# The rest of your tabs can follow here...
