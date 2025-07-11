import streamlit as st
import pandas as pd

# --------------------------------------
# Page Setup
# --------------------------------------
st.set_page_config(page_title="ViroMap", layout="wide")

# --------------------------------------
# Custom CSS – white bg + unified theme
# --------------------------------------
st.markdown("""
<style>
body {
    background-color: white;
}
h1, h2 {
    color: #AA336A;
}
.stTabs [data-baseweb="tab"] {
    background-color: #f8f8f8;
    border: 2px solid #ccc;
    border-radius: 8px;
    padding: 10px 16px;
    margin-right: 6px;
    font-weight: 600;
    color: #444;
}
.stTabs [aria-selected="true"] {
    background-color: #AA336A !important;
    color: white !important;
}
</style>
""", unsafe_allow_html=True)

# --------------------------------------
# App Header
# --------------------------------------
st.markdown("<h1 style='text-align: center;'>🧬 ViroMap – Integrated Viral Bioinformatics Portal</h1>", unsafe_allow_html=True)

# --------------------------------------
# Tabs
# --------------------------------------
tabs = st.tabs([
    "🏠 Overview",
    "🧬 Codon Bias",
    "💧 LLPS",
    "🧫 Mimicry",
    "🧠 GNN Predictions",
    "📤 Upload & Explore"
])

# --------------------------------------
# File Map (static datasets)
# --------------------------------------
file_map = {
    "🧬 Codon Bias":      "data/sars-cov-2_codon_bias.csv",
    "💧 LLPS":            "data/sars-cov-2_llps.csv",
    "🧫 Mimicry":         "data/sars-cov-2_mimicry.csv"
}

# --------------------------------------
# Tab 1 – Overview
# --------------------------------------
with tabs[0]:
    st.subheader("🏠 Overview")
    st.markdown("""
Welcome to **ViroMap**, a bioinformatics tool for COVID‑19 and viral research.
It combines data layers into one platform:

- 🧬 **Codon Bias** – Viral adaptation to host codons  
- 💧 **LLPS** – Intrinsically disordered regions & phase separation  
- 🧫 **Mimicry** – Viral peptides mimicking human proteins  
- 🧠 **GNN Predictions** – AI-predicted drug binding to spike protein  
- 📤 **Upload** – Bring your own CSV and explore

Use the tabs to explore and download datasets.
""")

# --------------------------------------
# Tab 2 – Codon Bias
# --------------------------------------
with tabs[1]:
    st.subheader("🧬 Codon Bias")
    df = pd.read_csv(file_map["🧬 Codon Bias"])
    query = st.text_input("🔍 Search Codon Bias", key="codon")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "codon_bias.csv")

# --------------------------------------
# Tab 3 – LLPS
# --------------------------------------
with tabs[2]:
    st.subheader("💧 LLPS Prediction")
    df = pd.read_csv(file_map["💧 LLPS"])
    query = st.text_input("🔍 Search LLPS", key="llps")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "llps.csv")

# --------------------------------------
# Tab 4 – Mimicry
# --------------------------------------
with tabs[3]:
    st.subheader("🧫 Molecular Mimicry")
    df = pd.read_csv(file_map["🧫 Mimicry"])
    query = st.text_input("🔍 Search Mimicry", key="mimicry")
    if query.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "mimicry.csv")

# --------------------------------------
# Tab 5 – GNN Predictions
# --------------------------------------
with tabs[4]:
    st.subheader("🧠 Real GNN Predictions (COVID Spike Protein)")

    try:
        gnn_df = pd.read_csv("data/real_gnn_predictions.csv")
        gnn_df["GNN_pKd"] = pd.to_numeric(gnn_df["GNN_pKd"], errors="coerce")
        search = st.text_input("🔍 Search Drug", key="gnn")
        if search.strip():
            filtered = gnn_df[gnn_df["Drug"].str.contains(search, case=False)]
        else:
            filtered = gnn_df.copy()
        filtered = filtered.dropna(subset=["GNN_pKd"])

        if not filtered.empty:
            def highlight(val):
                if pd.isna(val): return ""
                color = "green" if val > 7.2 else "orange" if val > 6.9 else "red"
                return f"color:{color}"
            st.dataframe(filtered.style.applymap(highlight, subset=["GNN_pKd"]), use_container_width=True)

            if len(filtered) > 1:
                st.subheader("📊 GNN Binding Scores")
                st.bar_chart(filtered.set_index("Drug")["GNN_pKd"])
            else:
                st.info("ℹ️ Not enough entries to show chart.")

            st.download_button(
                "📥 Download Filtered GNN CSV",
                filtered.to_csv(index=False),
                file_name="filtered_gnn_predictions.csv",
                mime="text/csv"
            )
        else:
            st.warning("No matching drugs found.")
    except FileNotFoundError:
        st.error("❌ `real_gnn_predictions.csv` not found in `/data/` folder.")

# --------------------------------------
# Tab 6 – Upload CSV
# --------------------------------------
with tabs[5]:
    st.subheader("📤 Upload & Explore Your Data")
    file_up = st.file_uploader("Upload a CSV", type="csv")
    if file_up:
        df = pd.read_csv(file_up)
        st.success("File uploaded successfully!")
        query = st.text_input("🔍 Search Uploaded Data", key="upload")
        if query.strip():
            df = df[df.apply(lambda r: r.astype(str).str.contains(query, case=False).any(), axis=1)]
        st.dataframe(df, use_container_width=True)
        st.download_button("📥 Download Filtered CSV", df.to_csv(index=False), "uploaded_filtered.csv")
    else:
        st.info("Upload a CSV file to explore it here.")
