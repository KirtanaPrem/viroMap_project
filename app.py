import streamlit as st
import pandas as pd

# --------------------------------------------------
#  Basic page config
# --------------------------------------------------
st.set_page_config(
    page_title="ViroMap",
    page_icon="🧬",
    layout="wide"
)

# --------------------------------------------------
#  Simple CSS for colorful tabs
# --------------------------------------------------
st.markdown("""
<style>
h1, h2 {color:#FF5C8D;}
.stTabs [data-baseweb="tab"] {
    background-color:#ffe6f0;
    border-radius:6px;
    padding:8px 14px;
    margin-right:6px;
    font-weight:600;
    color:#333333;
}
.stTabs [aria-selected="true"] {
    background-color:#FF5C8D !important;
    color:white !important;
}
</style>
""", unsafe_allow_html=True)

# --------------------------------------------------
#  App title
# --------------------------------------------------
st.markdown(
    "<h1 style='text-align:center'>🧬 ViroMap – Integrated AI & Bioinformatics Dashboard</h1>",
    unsafe_allow_html=True
)

# --------------------------------------------------
#  Define 6 colorful tabs
# --------------------------------------------------
tabs = st.tabs([
    "🏠 Overview",
    "🧬 Codon Bias",
    "💧 LLPS",
    "🧫 Mimicry",
    "🧠 GNN Predictions",
    "📤 Upload & Explore"
])

# --------------------------------------------------
#  File map for static CSV datasets
# --------------------------------------------------
file_map = {
    "🧬 Codon Bias":          "data/sars-cov-2_codon_bias.csv",
    "💧 LLPS":                "data/sars-cov-2_llps.csv",
    "🧫 Mimicry":             "data/sars-cov-2_mimicry.csv",
    "💊 Drug Predictions":    "data/sars-cov-2_drug_predictions.csv"
}

# --------------------------------------------------
#  TAB 1 – Overview
# --------------------------------------------------
with tabs[0]:
    st.subheader("🏠 Overview")
    st.markdown("""
**ViroMap** integrates multiple bioinformatics layers for viral research:

| Layer | What it shows |
|-------|---------------|
| Codon Bias | Viral adaptation to the host |
| LLPS | Phase‑separation regions in proteins |
| Mimicry | Viral peptides that imitate human proteins |
| **GNN Predictions** | AI‑predicted drug–protein binding strength |

Use the tabs above to explore each dataset, filter with search, or upload your own CSV.
    """)

# --------------------------------------------------
#  TAB 2 – Codon Bias
# --------------------------------------------------
with tabs[1]:
    st.subheader("🧬 Codon Bias")
    df = pd.read_csv(file_map["🧬 Codon Bias"])
    search = st.text_input("🔍 Search Codon Bias", key="codon")
    if search.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(search, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "codon_bias.csv")

# --------------------------------------------------
#  TAB 3 – LLPS
# --------------------------------------------------
with tabs[2]:
    st.subheader("💧 LLPS Prediction")
    df = pd.read_csv(file_map["💧 LLPS"])
    search = st.text_input("🔍 Search LLPS", key="llps")
    if search.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(search, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "llps.csv")

# --------------------------------------------------
#  TAB 4 – Mimicry
# --------------------------------------------------
with tabs[3]:
    st.subheader("🧫 Molecular Mimicry")
    df = pd.read_csv(file_map["🧫 Mimicry"])
    search = st.text_input("🔍 Search Mimicry", key="mimicry")
    if search.strip():
        df = df[df.apply(lambda r: r.astype(str).str.contains(search, case=False).any(), axis=1)]
    st.dataframe(df, use_container_width=True)
    st.download_button("📥 Download CSV", df.to_csv(index=False), "mimicry.csv")

# --------------------------------------------------
#  TAB 5 – REAL GNN Predictions
# --------------------------------------------------
with tabs[4]:
    st.subheader("🧠 Real GNN Predictions (Spike Protein)")

    try:
        gnn_df = pd.read_csv("data/real_gnn_predictions.csv")

        # ---- Search ----
        query = st.text_input("🔍 Search drug name", key="gnn")
        if query.strip():
            gnn_df = gnn_df[gnn_df["Drug"].str.contains(query, case=False)]

        if not gnn_df.empty:
            # ---- Color scores ----
            def color_score(v):
                col = 'green' if v > 7.2 else 'orange' if v > 6.9 else 'red'
                return f'color:{col}'
            st.dataframe(
                gnn_df.style.applymap(color_score, subset=["GNN_pKd"]),
                use_container_width=True
            )

            # ---- Bar chart (only if ≥1 row) ----
            if len(gnn_df) > 0:
                st.subheader("📊 Binding affinity (pKd)")
                st.bar_chart(gnn_df.set_index("Drug")["GNN_pKd"])

            # ---- Download ----
            st.download_button(
                "📥 Download filtered CSV",
                data=gnn_df.to_csv(index=False),
                file_name="gnn_predictions_filtered.csv",
                mime="text/csv"
            )
        else:
            st.warning("No drugs match your search.")

    except FileNotFoundError:
        st.error("`real_gnn_predictions.csv` not found in `/data/` folder.")

# --------------------------------------------------
#  TAB 6 – Upload & Explore
# --------------------------------------------------
with tabs[5]:
    st.subheader("📤 Upload & Explore Your CSV")
    file_up = st.file_uploader("Drop a CSV file here", type="csv")
    if file_up:
        up_df = pd.read_csv(file_up)
        st.success("File uploaded!")
        search_up = st.text_input("🔍 Search uploaded data", key="upload")
        if search_up.strip():
            up_df = up_df[up_df.apply(lambda r: r.astype(str).str.contains(search_up, case=False).any(), axis=1)]
        st.dataframe(up_df, use_container_width=True)
        st.download_button(
            "📥 Download filtered upload",
            data=up_df.to_csv(index=False),
            file_name="filtered_upload.csv",
            mime="text/csv"
        )
    else:
        st.info("Upload a CSV to view and search your own data.")
