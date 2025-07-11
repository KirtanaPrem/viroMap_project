import streamlit as st
import pandas as pd
import requests, time
from collections import Counter
from functools import lru_cache
from Bio import Entrez
import math
import xml.etree.ElementTree as ET

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
Entrez.email = "your-email@example.com"  # Replace with your real email

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

# ============  BLAST Mimicry ============
def run_mimicry_blast(spike_fasta):
    base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    payload = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "refseq_protein",
        "QUERY": spike_fasta,
        "ENTREZ_QUERY": "txid9606[Organism]",
        "FORMAT_TYPE": "XML"
    }

    r = requests.post(base_url, data=payload)
    if r.status_code != 200: return None

    lines = r.text.splitlines()
    rid = next((line.split("=")[1] for line in lines if line.startswith("RID")), None)
    if not rid:
        return None

    for _ in range(10):
        time.sleep(5)
        check = requests.get(base_url, params={"CMD": "Get", "RID": rid})
        if "Status=READY" in check.text:
            break

    result = requests.get(base_url, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"})
    tree = ET.fromstring(result.content)

    hits = []
    for hit in tree.findall(".//Hit"):
        hit_id = hit.findtext("Hit_def")
        acc = hit.findtext("Hit_accession")
        length = hit.findtext("Hit_len")
        hsps = hit.findall("Hit_hsps/Hsp")
        if hsps:
            hsp = hsps[0]
            identity = hsp.findtext("Hsp_identity")
            align_len = hsp.findtext("Hsp_align-len")
            evalue = hsp.findtext("Hsp_evalue")
            hits.append({
                "Protein": hit_id,
                "Accession": acc,
                "Length": length,
                "Identity": identity,
                "Align_Len": align_len,
                "E-Value": evalue
            })

    return pd.DataFrame(hits[:10])

# ============  Tabs Setup ============
tab_names = ["📄 FASTA", "🧬 Codon Bias", "💧 LLPS", "🧫 Mimicry", "🧠 GNN"]
tabs = st.tabs(tab_names)

# ============  Main Logic ============
if strain:
    fasta = fetch_fasta(strain)
    if fasta:
        # --- Tab 1: FASTA ---
        with tabs[0]:
            st.subheader("📄 Spike Protein FASTA")
            st.code(fasta, language="fasta")

        # --- Tab 2: Codon Bias ---
        with tabs[1]:
            st.subheader("🧬 Codon Usage Analysis")
            codon_df = calculate_codon_usage(fasta)
            st.dataframe(codon_df, use_container_width=True)

        # --- Tab 3: LLPS Prediction ---
        with tabs[2]:
            st.subheader("💧 LLPS Propensity")
            llps_df = predict_llps_simple(fasta)
            st.dataframe(llps_df, use_container_width=True)

        # --- Tab 4: Mimicry via BLAST ---
        with tabs[3]:
            st.subheader("🧫 Mimicry Prediction (Spike vs Human)")

            with st.spinner("Running BLAST... this takes ~30s ⏳"):
                mim_df = run_mimicry_blast(fasta)

            if mim_df is not None and not mim_df.empty:
                st.success("BLAST complete! Top human mimics:")
                st.dataframe(mim_df, use_container_width=True)
            else:
                st.error("BLAST failed or no hits found.")

        # --- Tab 5: GNN Placeholder ---
        with tabs[4]:
            st.subheader("🧠 GNN Drug Prediction")
            st.warning("GNN prediction coming next...")

    else:
        for t in tabs:
            with t:
                st.error("❌ Could not fetch spike protein. Try a different virus or strain.")
else:
    for t in tabs:
        with t:
            st.info("ℹ️ Enter a virus or strain above to get started.")
