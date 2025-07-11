import streamlit as st
import pandas as pd
import requests, time, math, xml.etree.ElementTree as ET
from collections import Counter
from functools import lru_cache
from Bio import Entrez

# ============ Theme & Layout ============
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

st.markdown("<h1 style='text-align:center'>🧬 ViroMap – Unified Viral Dashboard</h1>", unsafe_allow_html=True)

# ============ Search Input ============
strain = st.text_input("🔍  Type virus / strain (e.g. SARS‑CoV‑2 Wuhan):").strip()
Entrez.email = "your@email.com"  # change to your actual email

@lru_cache(maxsize=10)
def fetch_fasta(strain_name):
    try:
        h = Entrez.esearch(db="nucleotide", term=f"{strain_name} AND spike", retmax=1)
        rec = Entrez.read(h)
        if not rec["IdList"]: return None
        fid = rec["IdList"][0]
        return Entrez.efetch(db="nucleotide", id=fid, rettype="fasta", retmode="text").read()
    except Exception:
        return None

# ============ Codon Bias ============
def calculate_codon_usage(fasta_text):
    seq = "".join(fasta_text.splitlines()[1:]).upper()
    seq = seq[: len(seq) - len(seq) % 3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    cnt = Counter(codons)
    total = sum(cnt.values())
    return pd.DataFrame(
        [{"Codon": c, "Count": n, "Frequency": round(n / total, 4)} for c, n in sorted(cnt.items())]
    )

# ============ LLPS Predictor ============
def predict_llps_simple(fasta_text):
    seq = "".join(fasta_text.splitlines()[1:]).upper()
    hydro = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
             'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
             'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    rows = []
    for i, aa in enumerate(seq):
        h = hydro.get(aa, 0)
        dis = 1 if aa in 'PG' else 0
        score = round(h + 1.5 * dis, 2)
        rows.append([i+1, aa, score])
    df = pd.DataFrame(rows, columns=["Position", "AA", "LLPS_Score"])
    df["Region"] = pd.cut(df["LLPS_Score"], [-math.inf, 0, 2, math.inf], labels=["Low", "Medium", "High"])
    return df

# ============ BLAST Mimicry (safe) ============
def run_mimicry_blast(spike_fasta):
    base = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    put = {
        "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
        "QUERY": spike_fasta, "ENTREZ_QUERY": "txid9606[Organism]",
        "FORMAT_TYPE": "XML"
    }
    try:
        submit = requests.post(base, data=put, timeout=30)
        if submit.status_code != 200:
            return None
        rid = None
        for line in submit.text.splitlines():
            if line.startswith("RID"):
                rid = line.split("=", 1)[1].strip()
                break
        if not rid:
            return None
        for _ in range(12):
            time.sleep(5)
            status = requests.get(base, params={"CMD": "Get", "RID": rid}).text
            if "Status=READY" in status:
                break
        result = requests.get(base, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}).content
        root = ET.fromstring(result)
        hits = []
        for hit in root.findall(".//Hit")[:10]:
            hsp = hit.find("Hit_hsps/Hsp")
            if hsp is None: continue
            hits.append({
                "Protein": hit.findtext("Hit_def")[:60],
                "Accession": hit.findtext("Hit_accession"),
                "Identity": hsp.findtext("Hsp_identity"),
                "Align_Len": hsp.findtext("Hsp_align-len"),
                "E-value": hsp.findtext("Hsp_evalue")
            })
        return pd.DataFrame(hits)
    except Exception:
        return None

# ============ NetMHC Demo ============
def netmhc_demo():
    peptides = ["GVYYPDKVFR", "QPELDSFKEE", "SYGFQPTNGV", "VLSFELLHAP",
                "NLNESLIDLQ", "QLTPTWRVYS", "LTPGDSSSGW", "MESEFRVYSS"]
    return pd.DataFrame({
        "Peptide": peptides,
        "HLA": ["HLA-A*02:01"] * len(peptides),
        "Affinity(nM)": [50, 120, 600, 80, 400, 1500, 30, 200],
        "Rank(%)": [0.2, 0.5, 2.5, 0.3, 1.8, 5.0, 0.1, 1.2],
        "Binder": ["Strong", "Strong", "Weak", "Strong", "Weak", "No", "Strong", "Weak"]
    })

# ============ GNN Demo ============
def gnn_demo():
    return pd.DataFrame({
        "Drug": ["Remdesivir", "Molnupiravir", "Favipiravir"],
        "SMILES": ["CC(C)OC1=NC=C(C(=O)N)N1C", "CC1(CC(=O)NC(=O)N1C)C", "CC1=NC(=O)N(C)C(=N1)N"],
        "Predicted_pKd": [7.2, 6.8, 6.5]
    })

# ============ Tabs ============
tabNames = ["📄 FASTA", "🧬 Codon Bias", "💧 LLPS", "🧫 Mimicry", "🧪 Epitope Mimicry", "🧠 GNN"]
tabs = st.tabs(tabNames)

# ============ Logic ============
if strain:
    fasta = fetch_fasta(strain)
    if fasta:
        with tabs[0]:
            st.subheader("📄 FASTA Sequence")
            st.code(fasta, language="fasta")

        with tabs[1]:
            st.subheader("🧬 Codon Bias")
            st.dataframe(calculate_codon_usage(fasta), use_container_width=True)

        with tabs[2]:
            st.subheader("💧 LLPS Prediction")
            st.dataframe(predict_llps_simple(fasta), use_container_width=True)

        with tabs[3]:
            st.subheader("🧫 Mimicry via BLAST")
            with st.spinner("Running BLAST… Please wait ⏳"):
                mimic_df = run_mimicry_blast(fasta)
            if mimic_df is not None and not mimic_df.empty:
                st.success("Top hits from BLAST vs Human Proteins:")
                st.dataframe(mimic_df, use_container_width=True)
            else:
                st.warning("No strong sequence matches found. Mimicry could still occur at structure or immune level.")

        with tabs[4]:
            st.subheader("🧪 Epitope Mimicry (NetMHC demo)")
            st.dataframe(netmhc_demo(), use_container_width=True)
            st.caption("Lower affinity and rank % = stronger binding. These may mimic human immune responses.")

        with tabs[5]:
            st.subheader("🧠 GNN Drug Repurposing Predictions (demo)")
            st.dataframe(gnn_demo(), use_container_width=True)
            st.caption("Predicted pKd (higher = better binding).")
    else:
        for tab in tabs:
            with tab:
                st.error("❌ FASTA not found for this strain. Try another one.")
else:
    for tab in tabs:
        with tab:
            st.info("ℹ️ Please enter a virus or strain name in the search bar.")
