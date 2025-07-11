import streamlit as st
import pandas as pd

st.set_page_config(page_title="ViroMap", layout="wide")

# Tab setup
tabs = st.tabs([
    "🏠 Overview",
    "🧬 Codon Bias",
    "💧 LLPS Prediction",
    "🧫 Molecular Mimicry",
    "🧠 GNN Predictions"
])

# Tab 1 – Overview
with tabs[0]:
    st.title("ViroMap – A Unified Bioinformatics Platform for Viral Drug Discovery")
    st.markdown("""
    Welcome to **ViroMap**, an integrative AI + Bioinformatics dashboard for COVID-19 research.
    
    Navigate the tabs to explore:
    - Codon usage patterns
    - LLPS-prone regions in viral proteins
    - Molecular mimicry with host proteins
    - Real Graph Neural Network (GNN) drug-binding predictions
    """)

# Tab 2 – Codon Bias
with tabs[1]:
    st.header("🧬 Codon Bias Analysis")
    st.info("Upload viral gene sequences to analyze codon usage and host adaptation.")
    # Add your codon bias upload/input logic here

# Tab 3 – LLPS Prediction
with tabs[2]:
    st.header("💧 LLPS (Liquid-Liquid Phase Separation) Analysis")
    st.info("Explore LLPS-prone regions of the spike protein using prediction tools.")
    # Add your LLPS visual/chart code here

# Tab 4 – Molecular Mimicry
with tabs[3]:
    st.header("🧫 Molecular Mimicry Detection")
    st.info("Compare viral and human proteins to detect molecular mimicry.")
    # Add your mimicry result table or upload option here

# Tab 5 – GNN Predictions
with tabs[4]:
    st.header("🧠 Real GNN Predictions for COVID-19 Spike Protein")

    # Load predictions CSV you uploaded into /data/
    try:
        gnn_df = pd.read_csv("data/real_gnn_predictions.csv")
        st.dataframe(gnn_df, use_container_width=True)

        st.subheader("📊 GNN-pKd Score Visualization")
        st.bar_chart(gnn_df.set_index("Drug")["GNN_pKd"])

        st.download_button(
            label="📥 Download Predictions CSV",
            data=gnn_df.to_csv(index=False),
            file_name="real_gnn_predictions.csv",
            mime="text/csv"
        )
    except FileNotFoundError:
        st.error("real_gnn_predictions.csv not found. Please upload it to the /data/ folder in GitHub.")
