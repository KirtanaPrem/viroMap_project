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

    try:
        gnn_df = pd.read_csv("data/real_gnn_predictions.csv")

        st.write("Column names:", gnn_df.columns.tolist())  # TEMPORARY check

        # Update this based on your actual column name
        score_column = "GNN_pKd"  # change if it's different in your CSV

        # Search bar
        search_term = st.text_input("🔍 Search Drug", "")
        filtered_df = gnn_df[gnn_df["Drug"].str.contains(search_term, case=False)]

        if not filtered_df.empty:
            # Color score based on strength
            def highlight_score(val):
                color = 'green' if val > 7.2 else 'orange' if val > 6.9 else 'red'
                return f'color: {color}'

            st.subheader("📋 Prediction Table")
            st.dataframe(filtered_df.style.applymap(highlight_score, subset=[score_column]), use_container_width=True)

            # Bar chart
            st.subheader("📊 Binding Score Chart")
            st.bar_chart(filtered_df.set_index("Drug")[score_column])

            # Download
            st.download_button(
                label="📥 Download CSV",
                data=filtered_df.to_csv(index=False),
                file_name="real_gnn_predictions.csv",
                mime="text/csv"
            )
        else:
            st.warning("No drugs match your search.")

    except FileNotFoundError:
        st.error("real_gnn_predictions.csv not found in `/data/` folder.")


    except FileNotFoundError:
        st.error("Prediction file not found. Please upload `real_gnn_predictions.csv` to `/data/` folder in GitHub.")

