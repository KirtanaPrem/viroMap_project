
import streamlit as st
import pandas as pd

st.title("🧬 ViroMap: Virus Feature Integration for Drug Targeting")

st.sidebar.title("Select Dataset")
dataset = st.sidebar.selectbox("Choose a dataset", [
    "Codon Bias", "LLPS", "Molecular Mimicry", "Drug Predictions"
])

file_map = {
    "Codon Bias": "data/sars-cov-2_codon_bias.csv",
    "LLPS": "data/sars-cov-2_llps.csv",
    "Molecular Mimicry": "data/sars-cov-2_mimicry.csv",
    "Drug Predictions": "data/sars-cov-2_drug_predictions.csv"
}

df = pd.read_csv(file_map[dataset])
st.dataframe(df)
