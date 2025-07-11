import streamlit as st
import pandas as pd

st.set_page_config(
    page_title="ViroMap",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 ViroMap: Virus Feature Integration for Drug Targeting")

# Create tabs
tabs = st.tabs(["Drug Predictions", "Codon Bias", "LLPS", "Molecular Mimicry", "GNN Prediction"])

# Define file paths
data_files = {
    "Drug Predictions": "data/sars-cov-2_drug_predictions.csv",
    "Codon Bias": "data/sars-cov-2_codon_bias.csv",
    "LLPS": "data/sars-cov-2_llps.csv",
    "Molecular Mimicry": "data/sars-cov-2_mimicry.csv"
}

# Load and display CSVs in respective tabs
for i, tab in enumerate(tabs[:4]):
    with tab:
        name = list(data_files.keys())[i]
        df = pd.read_csv(data_files[name])
        st.subheader(name)
        st.dataframe(df, use_container_width=True)

        # Add basic chart if "Predicted_Efficacy" column exists
        if "Predicted_Efficacy" in df.columns:
            st.bar_chart(df.set_index("Drug")["Predicted_Efficacy"])

# GNN Prediction tab
with tabs[4]:
    st.subheader("Graph Neural Network Insights")
    st.markdown("""
    This section displays a simulated output from a GNN model predicting virus-host interactions.
    """)

    gnn_data = pd.DataFrame({
        "Viral Protein": ["NSP3", "Spike", "RdRp"],
        "Predicted Target": ["Human ACE2", "TMPRSS2", "IFNAR1"],
        "Confidence Score": [0.91, 0.88, 0.84]
    })

    st.dataframe(gnn_data)
    st.bar_chart(gnn_data.set_index("Viral Protein")["Confidence Score"])
