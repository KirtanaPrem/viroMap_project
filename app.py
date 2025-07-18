import streamlit as st
import pandas as pd
from collections import Counter
from math import exp, log
from Bio import Entrez

# Configuration
Entrez.email = "your_email@example.com"  # Replace with your actual email
st.set_page_config(page_title="ViroMap", layout="wide")
st.title("ViroMap")

# Constants
RARE_CODONS = {"ATA", "AGG", "AGA", "CGG", "CTA", "TTA"}
DISORDERED_AAS = set("PESQKRDG")
PRION_MOTIFS = ["QQ", "QN", "NQ", "SS", "YG", "RG"]
CAI_REFERENCE = {
    "GCT": 0.18, "GCC": 0.55, "GCA": 0.15, "GCG": 0.12,
    "CGT": 0.04, "CGC": 0.11, "CGA": 0.11, "CGG": 0.21,
    "AGA": 0.21, "AGG": 0.32
}

def clean_sequence(fasta_str):
    """Extract and clean sequence from FASTA string."""
    if not fasta_str or not fasta_str.startswith(">"):
        return ""
    return "".join(fasta_str.splitlines()[1:]).replace(" ", "").replace("\n", "").upper()

def fetch_strain_options(keyword):
    """Fetch NCBI strain options based on keyword."""
    try:
        with Entrez.esearch(db="nucleotide", term=f"{keyword} AND complete genome", retmax=5) as handle:
            search_results = Entrez.read(handle)
            id_list = search_results["IdList"]

        if not id_list:
            return [], {}

        with Entrez.esummary(db="nucleotide", id=",".join(id_list)) as handle:
            summaries = Entrez.read(handle)
            
        strain_options = []
        strain_dict = {}
        for docsum in summaries:
            uid = docsum["Id"]
            title = docsum["Title"]
            label = f"{uid} | {title[:80]}..."
            strain_options.append(label)
            strain_dict[label] = uid
        return strain_options, strain_dict
        
    except Exception as e:
        st.error(f"Error fetching strains: {e}")
        return [], {}

def fetch_fasta(ncbi_id):
    """Fetch FASTA sequence by NCBI ID."""
    try:
        with Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text") as handle:
            return handle.read()
    except Exception as e:
        st.error(f"Error fetching FASTA: {e}")
        return ""

# Main UI
keyword = st.text_input("Enter virus keyword (e.g., rabies, influenza, SARS-CoV-2):")

strain_options, strain_dict = fetch_strain_options(keyword) if keyword else ([], {})
selected_strain = st.selectbox("Select a strain", strain_options) if strain_options else None
fasta = fetch_fasta(strain_dict.get(selected_strain, "")) if selected_strain else ""

# Tab setup
tabs = st.tabs(["FASTA", "Codon Bias", "LLPS", "Epitope Mimicry", "GNN Prediction"])

# Tab 0: FASTA
with tabs[0]:
    st.header("FASTA Sequence")
    if fasta and fasta.startswith(">"):
        st.code(fasta, language="fasta")
    elif fasta:
        st.warning("No valid FASTA found.")
    else:
        st.info("Search and select a strain to view its sequence.")

# Tab 1: Codon Bias
with tabs[1]:
    st.header("Codon Bias Metrics")

    clean_seq = clean_sequence(fasta)
    if not clean_seq:
        st.warning("No valid sequence found. Please try a different strain.")
        st.stop()

    def get_codons(seq):
        """Extract codons from DNA sequence."""
        return [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]

    def calculate_rare_codon_freq(codons):
        """Calculate frequency of rare codons."""
        if not codons:
            return 0
        return round(sum(1 for c in codons if c in RARE_CODONS) / len(codons), 4)

    def calculate_codon_pair_bias(codons):
        """Calculate codon pair bias score."""
        if len(codons) < 2:
            return 0
            
        pairs = [codons[i]+codons[i+1] for i in range(len(codons)-1)]
        expected = 1 / (64 * 64)
        counts = Counter(pairs)
        total = sum(counts.values())
        observed = [c / total for c in counts.values()]
        cpb = sum((o - expected) ** 2 / expected for o in observed)
        return round(cpb, 4)

    def calculate_volatility(codons):
        """Calculate codon volatility."""
        if not codons:
            return 0
            
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

    def calculate_cai(seq):
        """Calculate simplified Codon Adaptation Index."""
        codons = get_codons(seq)
        weights = [CAI_REFERENCE.get(codon, 0.01) for codon in codons]
        return round(exp(sum(log(w) for w in weights) / len(weights), 4) if weights else 0

    def calculate_enc(codons):
        """Calculate Effective Number of Codons."""
        if len(codons) < 2:
            return 0
            
        freqs = Counter(codons)
        n = sum(freqs.values())
        homo = sum(c * (c - 1) for c in freqs.values())
        return round(61 / (1 + homo / (n * (n - 1))), 2)

    codons = get_codons(clean_seq)
    
    metrics = {
        "Codon Adaptation Index (CAI)": calculate_cai(clean_seq),
        "Rare Codon Frequency": calculate_rare_codon_freq(codons),
        "Codon Pair Bias (CPB)": calculate_codon_pair_bias(codons),
        "Codon Volatility": calculate_volatility(codons),
        "Effective Number of Codons (ENC)": calculate_enc(codons)
    }

    st.dataframe(pd.DataFrame(list(metrics.items()), columns=["Metric", "Value"]))

# Tab 2: LLPS
with tabs[2]:
    st.header("LLPS Prediction")

    clean_seq = clean_sequence(fasta)
    if not clean_seq:
        st.warning("No valid sequence found. Please try a different strain.")
        st.stop()

    def calculate_disorder(seq):
        """Calculate percentage of disordered amino acids."""
        if not seq:
            return 0
        disordered = sum(1 for aa in seq if aa in DISORDERED_AAS)
        return round(disordered / len(seq) * 100, 2)

    def detect_prion_like(seq):
        """Check for prion-like motifs."""
        return any(motif in seq for motif in PRION_MOTIFS)

    def detect_lcr(seq, window=6, max_unique=2):
        """Count low complexity regions."""
        if len(seq) < window:
            return 0
            
        count = 0
        for i in range(len(seq)-window+1):
            if len(set(seq[i:i+window])) <= max_unique:
                count += 1
        return count

    def calculate_llps_score(seq):
        """Calculate LLPS propensity score."""
        if not seq:
            return 0
            
        score = 0
        disorder = calculate_disorder(seq)
        score += disorder / 100
        score += detect_lcr(seq) / 50
        if detect_prion_like(seq):
            score += 0.3
        return round(min(score, 1.0), 2)

    protein_seq = clean_seq.replace("T", "U")  # Pseudo translation for LLPS analysis
    
    data = {
        "LLPS Propensity Score (0-1)": calculate_llps_score(protein_seq),
        "Percent Disorder": calculate_disorder(protein_seq),
        "Prion-like Motifs Present": "Yes" if detect_prion_like(protein_seq) else "No",
        "Low Complexity Regions": detect_lcr(protein_seq)
    }

    st.dataframe(pd.DataFrame(list(data.items()), columns=["Metric", "Value"]))

# Tab 3: Epitope Mimicry
with tabs[3]:
    st.header("Epitope Mimicry")
    st.info("This module will identify host-virus mimicry (coming soon).")

# Tab 4: GNN Prediction
with tabs[4]:
    st.header("GNN Drug Prediction")
    st.info("This module will use Graph Neural Networks to predict antiviral drugs (coming soon).")
