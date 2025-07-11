import streamlit as st
import pandas as pd

st.set_page_config(page_title="ViroMap", page_icon="🧬", layout="wide")

st.markdown("<h1 style='text-align: center; color: #40E0D0;'>🧬 ViroMap: Virus Feature Integration for Drug Targeting</h1>", unsafe_allow_html=True)

st.sidebar.title("🧠 Navigation Panel")
dataset = st.sidebar.selectbox("Choose a dataset", [
    "Codon Bias", "LLPS", "Molecular Mimicry", "Drug Predictions", "Upload Your Own CSV"
])

file_map = {
    "Codon Bias": "data/sars-cov-2_codon_bias.csv",
    "LLPS": "data/sars-cov-2_llps.csv",
    "Molecular Mimicry": "data/sars-cov-2_mimicry.csv",
    "Drug Predictions": "data/sars-cov-2_drug_predictions.csv"
}

# Upload Option
if dataset == "Upload Your Own CSV":
    uploaded_file = st.sidebar.file_uploader("📂 Upload your CSV file", type="csv")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
    else:
        st.warning("Please upload a CSV file.")
        st.stop()
else:
    df = pd.read_csv(file_map[dataset])

# Search Option
search_query = st.text_input("🔍 Search the table (by keyword)", "")
if search_query:
    df = df[df.apply(lambda row: row.astype(str).str.contains(search_query, case=False).any(), axis=1)]

st.dataframe(df, use_container_width=True)
