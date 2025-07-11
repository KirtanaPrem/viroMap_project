import streamlit as st
import pandas as pd

st.set_page_config(
    page_title="ViroMap",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
    <style>
        .main { background-color: #f9f9f9; }
        h1, h2, h3 {
            color: #FF5C8D;
        }
        .stTabs [data-baseweb="tab"] {
            background-color: #ffe6f0;
            border-radius: 8px;
            padding: 10px;
            margin-right: 5px;
            color: black;
            font-weight: bold;
        }
        .stTabs [aria-selected="true"] {
            background-color: #FF5C8D !important;
            color: white !important;
        }
    </style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='text-align: center;'>🧬 ViroMap: Virus Feature Integration Platform</h1>", unsafe_allow_html=True)

# Create colorful tabs
tabs = st.tabs([
    "📁 Codon Bias", 
    "💧 LLPS", 
    "🔬 Molecular Mimicry", 
    "💊 Drug Predictions", 
    "🧠 GNN Predictions", 
    "📤 Upload & Explore"
])

# Load CSV files
file_map = {
    "📁 Codon Bias": "data/sars-cov-2_codon_bias.csv",
    "💧 LLPS": "data/sars-cov-2_llps.csv",
    "🔬 Molecular Mimicry": "data/sars-cov-2_mimicry.csv",
    "💊 Drug Predictions": "data/sars-cov-2_drug_predictions.csv"
}

# Tabs 1–4: Display Datasets
for i, tab in enumerate(tabs[:4]):
    dataset_name = list(file_map.keys())[i]
    with tab:
        st.header(dataset_name)
        df = pd.read_csv(file_map[dataset_name])
        
        search = st.text_input(f"🔍 Search in {dataset_name}", key=i)
        if search:
            df = df[df.apply(lambda row: row.astype(str).str.contains(search, case=False).any(), axis=1)]

        st.dataframe(df, use_container_width=True)
        st.download_button(f"📥 Download {dataset_name}", df.to_csv(index=False), f"{dataset_name}.csv")

# Tab 5: GNN Predictions
with tabs[4]:
    st.header("🧠 GNN Predictions")
    st.info("This is a simulated output from a Graph Neural Network. Actual model integration is in progress.")
    
    gnn_data = pd.DataFrame({
        "Drug": ["Remdesivir", "Favipiravir", "Molnupiravir"],
        "Predicted Target": ["NSP12", "RdRp", "Spike"],
        "GNN Confidence Score": [0.92, 0.88, 0.91]
    })
    st.dataframe(gnn_data, use_container_width=True)
    st.download_button("📥 Download GNN Output", gnn_data.to_csv(index=False), "gnn_predictions.csv")

# Tab 6: Upload & Explore
with tabs[5]:
    st.header("📤 Upload Your Own CSV")
    uploaded = st.file_uploader("Upload a CSV to explore", type="csv")
    if uploaded:
        user_df = pd.read_csv(uploaded)
        st.success("File uploaded successfully!")
        user_search = st.text_input("🔍 Search in uploaded file")
        if user_search:
            user_df = user_df[user_df.apply(lambda row: row.astype(str).str.contains(user_search, case=False).any(), axis=1)]
        st.dataframe(user_df, use_container_width=True)
        st.download_button("📥 Download Filtered Results", user_df.to_csv(index=False), "filtered_upload.csv")
    else:
        st.warning("No file uploaded yet.")
