import streamlit as st
import pandas as pd
import requests, time, math, xml.etree.ElementTree as ET
from collections import Counter
from functools import lru_cache
from Bio import Entrez
from DeepPurpose import utils, onedrug

# ------------------- CONFIG & STYLE -------------------
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

# ------------------- NCBI Search -------------------
Entrez.email = "your@email.com"  # Replace with your email

strain = st.text_input("🔍 Type virus / strain (e.g. SARS-CoV-2, HIV-1):").strip()
strain_matches = []
selected_id = None
fasta = None

@lru_cache(maxsize=10)
def get_strain_matches(query):
    try:
        handle = Entrez.esearch(db="nucleotide", term=f'"{query}"[Organism] AND biomol_mrna[PROP]', retmax=10)
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids:
            return []
        summaries = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summary_data = Entrez.read(summaries)
        return [(doc["Title"], doc["Id"]) for doc in summary_data]
    except:
        return []

def fetch_fasta_by_id(ncbi_id):
    try:
        return Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text").read()
    except:
        return None

if strain:
    with st.spinner("🔎 Searching NCBI..."):
        strain_matches = get_strain_matches(strain)

    if strain_matches:
        titles = [f"{title[:70]}..." if len(title) > 70 else title for title, _ in strain_matches]
        selected = st.selectbox("🧬 Select matching strain from NCBI:", titles)
        selected_id = strain_matches[titles.index(selected)][1]
        fasta = fetch_fasta_by_id(selected_id)
    else:
        st.warning("❌ No matching strains found in NCBI.")

# ------------------- Codon Bias -------------------
def calculate_codon_usage(fasta_text):
    seq = "".join(fasta_text.splitlines()[1:]).upper()
    seq = seq[: len(seq) - len(seq) % 3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    cnt = Counter(codons)
    total = sum(cnt.values())
    return pd.DataFrame(
        [{"Codon": c, "Count": n, "Freq": round(n / total, 4)} for c, n in sorted(cnt.items())]
    )

# ------------------- LLPS -------------------
def predict_llps_simple(fasta_text):
    hydro = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
             'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
             'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    seq = "".join(fasta_text.splitlines()[1:]).upper()
    rows=[]
    for i,aa in enumerate(seq):
        h=hydro.get(aa,0)
        dis=1 if aa in 'PG' else 0
        score=round(h+1.5*dis,2)
        rows.append([i+1,aa,score])
    df=pd.DataFrame(rows,columns=["Pos","AA","LLPS_Score"])
    df["Region"]=pd.cut(df["LLPS_Score"],[-math.inf,0,2,math.inf],labels=["Low","Medium","High"])
    return df

# ------------------- Mimicry BLAST -------------------
def run_mimicry_blast(fasta_text):
    base = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    put = {
        "CMD":"Put","PROGRAM":"blastp","DATABASE":"nr",
        "QUERY":fasta_text,"ENTREZ_QUERY":"txid9606[Organism]",
        "FORMAT_TYPE":"XML"
    }
    try:
        r=requests.post(base,data=put,timeout=30)
        if r.status_code!=200: return None
        rid=None
        for line in r.text.splitlines():
            if line.startswith("RID"):
                rid=line.split("=",1)[1].strip(); break
        if not rid: return None
        for _ in range(12):
            time.sleep(5)
            if "Status=READY" in requests.get(base,params={"CMD":"Get","RID":rid}).text:
                break
        xml=requests.get(base,params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"XML"}).content
        root=ET.fromstring(xml)
        hits=[]
        for hit in root.findall(".//Hit")[:10]:
            hsp=hit.find("Hit_hsps/Hsp")
            if hsp is None: continue
            hits.append({
                "Protein":hit.findtext("Hit_def")[:60],
                "Acc":hit.findtext("Hit_accession"),
                "Identity":hsp.findtext("Hsp_identity"),
                "AlignLen":hsp.findtext("Hsp_align-len"),
                "E-value":hsp.findtext("Hsp_evalue")
            })
        return pd.DataFrame(hits)
    except:
        return None

# ------------------- NetMHC Demo -------------------
def netmhc_demo():
    peptides=["GVYYPDKVFR","QPELDSFKEE","SYGFQPTNGV","VLSFELLHAP"]
    return pd.DataFrame({
        "Peptide":peptides,"HLA":["HLA-A*02:01"]*len(peptides),
        "Affinity(nM)":[50,120,600,80],
        "Rank%":[0.2,0.5,2.5,0.3],
        "Binder":["Strong","Strong","Weak","Strong"]
    })

# ------------------- GNN Real Prediction -------------------
@st.cache_resource
def load_model():
    drug_encoding = 'SMILES'
    target_encoding = 'AAC'
    return onedrug.model_initialize(drug_encoding, target_encoding, pretrained=True)

@st.cache_data
def gnn_predict(fasta_seq):
    target_seq = "".join(fasta_seq.splitlines()[1:]).strip()
    drugs = {
        "Remdesivir": "CCC(C)OC(=O)[C@H](C)NP(=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1",
        "Molnupiravir": "CC1=CC=C(C=C1)C(CNC(=O)OC)OC(=O)C(C)C",
        "Favipiravir": "NC(=O)c1nc(F)c[nH]c1=O"
    }
    model = load_model()
    drug_list = list(drugs.values())
    drug_names = list(drugs.keys())
    df = pd.DataFrame([[d, target_seq] for d in drug_list], columns=["SMILES", "Target Sequence"])
    test_data = utils.data_process(df, drug_encoding='SMILES', target_encoding='AAC',
                                   split_method='no_split', mode='pKd')
    scores = model.predict(test_data)
    return pd.DataFrame({
        "Drug": drug_names,
        "pKd": [round(s, 2) for s in scores]
    }).sort_values("pKd", ascending=False)

# ------------------- Tabs -------------------
tab_names=["📄 FASTA","🧬 Codon Bias","💧 LLPS","🧫 Mimicry","🧪 Epitope","🧠 GNN"]
tabs=st.tabs(tab_names)

if fasta:
    tabs[0].code(fasta, language="fasta")
    tabs[1].dataframe(calculate_codon_usage(fasta), use_container_width=True)
    tabs[2].dataframe(predict_llps_simple(fasta), use_container_width=True)

    with tabs[3]:
        st.subheader("🧫 BLAST vs Human (Mimicry)")
        with st.spinner("Running BLAST (~30 s)…"):
            mdf = run_mimicry_blast(fasta)
        if mdf is not None and not mdf.empty:
            st.dataframe(mdf, use_container_width=True)
        else:
            st.warning("No close sequence matches found.")

    tabs[4].dataframe(netmhc_demo(), use_container_width=True)
    tabs[4].caption("Lower Affinity & Rank% = stronger HLA binding (demo).")

    with tabs[5]:
        st.subheader("🧠 GNN-Based Drug Binding Prediction")
        with st.spinner("Predicting using DeepPurpose..."):
            gdf = gnn_predict(fasta)
        st.dataframe(gdf, use_container_width=True)
        st.caption("Predicted pKd: Higher = stronger binding affinity")
else:
    for t in tabs:
        with t:
            st.info("ℹ️ Search a virus and select a matching strain to continue.")
