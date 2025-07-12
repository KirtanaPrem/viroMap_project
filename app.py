import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
import requests, time, xml.etree.ElementTree as ET
from collections import Counter

# ------------ PAGE LOOK ------------
st.set_page_config(page_title="ViroMap", layout="wide", page_icon="🧬")
st.markdown(
    """
    <style>
        body { background:white; }
        .stTabs [data-baseweb="tab"]{
            width: 16%;                 /* equal width tabs */
            text-align: center;
            background:#e3f2fd;
            padding:10px;
            border-radius:8px;
            margin-right:4px;
        }
        .stTabs [aria-selected="true"]{
            background:#90caf9;
            font-weight:bold;
        }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🧬 ViroMap – Unified Viral Prediction Dashboard")

# ------------ SEARCH ------------
query = st.text_input("🔍  Enter virus name (e.g. SARS‑CoV‑2, HIV‑1, Dengue):")

Entrez.email = "example@email.com"        # <-- put your real e‑mail

def ncbi_search_strains(term, retmax=10):
    """Return list of (id, description) for matching nucleotide records."""
    try:
        h = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
        rec = Entrez.read(h)
        ids = rec["IdList"]
        strains = []
        for id_ in ids:
            fetch = Entrez.efetch(db="nucleotide", id=id_, rettype="gb", retmode="text")
            rec_gb = SeqIO.read(fetch, "genbank")
            strains.append((id_, rec_gb.description))
        return strains
    except Exception:
        return []

def fetch_fasta(ncbi_id):
    try:
        f = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
        rec = SeqIO.read(f, "fasta")
        return str(rec.seq), rec.id
    except Exception:
        return "", ""

# ------------ REAL MIMICRY (BLASTP) ------------
def blastp_mimicry(seq, top=10):
    """BLAST spike AA sequence vs human proteins – return top hits DataFrame."""
    base = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    # Prepare BLASTP job
    put_data = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": seq,
        "ENTREZ_QUERY": "txid9606[Organism]",   # human proteins only
        "FORMAT_TYPE": "XML"
    }
    r = requests.post(base, data=put_data, timeout=30)
    if r.status_code != 200:
        return pd.DataFrame()

    # Get RID
    rid = None
    for line in r.text.splitlines():
        if line.startswith("RID"):
            rid = line.split("=", 1)[1].strip()
            break
    if not rid:
        return pd.DataFrame()

    # Poll until ready (max ~1 min)
    for _ in range(12):
        time.sleep(5)
        chk = requests.get(base, params={"CMD": "Get", "RID": rid})
        if "Status=READY" in chk.text:
            break

    xml = requests.get(base, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}).content
    root = ET.fromstring(xml)
    hits = []
    for hit in root.findall(".//Hit")[:top]:
        hsp = hit.find("Hit_hsps/Hsp")
        if hsp is None:
            continue
        hits.append({
            "Human Protein": hit.findtext("Hit_def")[:60],
            "Acc": hit.findtext("Hit_accession"),
            "Identity": hsp.findtext("Hsp_identity"),
            "Align_Len": hsp.findtext("Hsp_align-len"),
            "E‑value": hsp.findtext("Hsp_evalue")
        })
    return pd.DataFrame(hits)

# ------------ QUICK DEMO FUNCTIONS ------------
def codon_bias(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    cnt = Counter(codons)
    total = sum(cnt.values()) or 1
    return pd.DataFrame(
        [{"Codon": c, "Count": n, "Freq": round(n/total, 4)} for c, n in cnt.items()]
    ).sort_values("Freq", ascending=False)

def llps_demo(seq):
    return pd.DataFrame({"Region": ["N‑term", "Mid", "C‑term"], "LLPS score": [0.82, 0.55, 0.71]})

def epitope_demo(seq):
    return pd.DataFrame({"Epitope": ["SYGFQPT", "GFQPTNG"], "MHC Affinity": [84, 62]})

def gnn_demo():
    return pd.DataFrame({"Drug": ["Remdesivir", "Molnupiravir"], "pKd": [7.2, 6.8]})

# ------------ SEARCH & FETCH ------------
sequence = ""; rec_id = ""

if query:
    strains = ncbi_search_strains(query)
    if strains:
        choice = st.selectbox("🧬  Select strain:", {d: i for i, d in strains}.keys())
        sel_id = {d: i for i, d in strains}[choice]
        sequence, rec_id = fetch_fasta(sel_id)
    else:
        st.warning("No matching strains found in NCBI.")

# ------------ TABS ------------
tabs = st.tabs(["FASTA", "Codon", "LLPS", "Mimicry", "Epitope", "GNN"])

if sequence:
    with tabs[0]:
        st.code(sequence, language="fasta")
        st.success(f"NCBI record: {rec_id}")

    with tabs[1]:
        st.dataframe(codon_bias(sequence), height=300)

    with tabs[2]:
        st.dataframe(llps_demo(sequence))

    with tabs[3]:
        st.subheader("Running BLAST vs human proteome…")
        with st.spinner("⏳"):
            mimic_df = blastp_mimicry(sequence)
        if not mimic_df.empty:
            st.dataframe(mimic_df)
        else:
            st.warning("No strong sequence mimicry hits found.")

    with tabs[4]:
        st.dataframe(epitope_demo(sequence))

    with tabs[5]:
        st.caption("Demo GNN results (real model disabled in Cloud).")
        st.dataframe(gnn_demo())
else:
    st.info("Type a virus name and pick a strain to begin.")
