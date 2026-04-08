import streamlit as st
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import io

# ✅ SAFE IMPORT (FIXED)
try:
    import py3Dmol
    from stmol import showmol
    PDB_AVAILABLE = True
except:
    PDB_AVAILABLE = False

st.set_page_config(layout="wide")

plt.rcParams.update({
    'figure.autolayout': True,
    'font.size': 8
})

# =========================
# STYLE
# =========================
st.markdown("""
<style>
.stButton>button {
    border-radius: 8px;
    background-color: #2E86C1;
    color: white;
    font-weight: 500;
}
h1, h2, h3 {
    color: #1B4F72;
}
</style>
""", unsafe_allow_html=True)

# =========================
# LOAD DATA
# =========================
@st.cache_data
def load_data():
    return pd.read_csv("nipah_dataset1.csv")

df = load_data()

# =========================
# SIDEBAR
# =========================
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", [
    "Home","Tool","Vaccine Pipeline",
    "Advanced Analytics","Export Results",
    "About","Team"
])

# =========================
# HOME
# =========================
if page == "Home":
    st.title("Nipah Virus Vaccine Design Platform")
    st.markdown("### Computational Immunoinformatics-Based Vaccine Development System")
    st.markdown("---")

    st.subheader("Overview")
    st.info("""
    This platform is designed for in-silico multi-epitope vaccine design targeting the Nipah virus.
    """)

    st.subheader("Why Nipah Virus?")
    st.write("High mortality zoonotic virus. Computational methods accelerate vaccine discovery.")

    st.subheader("Workflow")
    st.write("""
    1. Target Selection  
    2. Pre-Screening  
    3. Epitope Analysis  
    4. Filtering  
    5. Vaccine Construction  
    6. Evaluation  
    7. Final Selection  
    """)

# =========================
# TOOL
# =========================
if page == "Tool":
    st.title("Tool Description")

    st.markdown("---")

    st.subheader("Overview")
    st.write("""
    This application is a Streamlit-based computational platform designed for 
    multi-epitope vaccine design against Nipah virus using immunoinformatics approaches.
    """)

    st.subheader("Methodology")
    st.write("""
    - Protein Selection  
    - Antigenicity-based screening  
    - Epitope classification (B-cell / T-cell)  
    - Filtering based on affinity and conservation  
    - Vaccine construction using linkers  
    """)

    st.subheader("Key Parameters")
    st.write("""
    - Antigenicity  
    - Binding Affinity  
    - Toxicity  
    - Allergenicity  
    - Conservation  
    """)

    st.subheader("Technologies Used")
    st.write("""
    Streamlit, Pandas, Matplotlib, Py3Dmol, Stmol
    """)

    st.success("Output: Ranked vaccine candidates for further validation")

# =========================
# PIPELINE
# =========================
if page == "Vaccine Pipeline":

    st.title("Vaccine Design Pipeline")

    if st.checkbox("Show Dataset"):
        st.dataframe(df)

    # TARGET
    st.header("1. Target Selection")
    protein = st.selectbox("Select Protein", df["Protein"].unique())
    df_protein = df[df["Protein"] == protein]
    st.dataframe(df_protein)

    # GRAPH
    if "Subcellular_Location" in df_protein.columns:
        fig, ax = plt.subplots(figsize=(3.5,2.2))
        df_protein["Subcellular_Location"].value_counts().plot(kind='bar', ax=ax)
        st.pyplot(fig)

        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        st.download_button("Download Graph", buf.getvalue(), "localization.png")

    # 3D STRUCTURE
    st.header("Protein Structure")

    if "PDB_ID" in df.columns:
        pdb_ids = df_protein["PDB_ID"].dropna().unique()

        if len(pdb_ids) > 0:
            selected_pdb = st.selectbox("Select PDB", pdb_ids)

            if st.button("Show Structure"):
                if PDB_AVAILABLE:
                    view = py3Dmol.view(query=f"pdb:{selected_pdb}")
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.zoomTo()
                    showmol(view, height=350)
                else:
                    st.warning("⚠️ 3D structure not supported here")
                    st.markdown(f"🔗 https://www.rcsb.org/structure/{selected_pdb}")

    # PRE-SCREENING
    st.header("2. Pre-Screening")
    thr = st.slider("Antigenicity Threshold", 0.0, 1.0, 0.7)

    if st.button("Run Pre-Screening"):
        filtered = df_protein[
            (df_protein["Antigenicity"] >= thr) &
            (df_protein["Allergenicity"] == "Non-allergen") &
            (df_protein["Toxicity"] == "Non-toxic")
        ]
        st.session_state["filtered"] = filtered
        st.dataframe(filtered)

    if "filtered" in st.session_state:
        fig, ax = plt.subplots(figsize=(3.5,2.2))
        ax.hist(st.session_state["filtered"]["Antigenicity"])
        st.pyplot(fig)

    # FILTER
    st.header("3. Filtering")

    if st.button("Filter"):
        if "filtered" in st.session_state:
            final = st.session_state["filtered"]
            st.session_state["final"] = final
            st.dataframe(final)

    # VACCINE
    st.header("4. Vaccine Construction")

    def construct(combo):
        return "EAAAK" + "AAY".join(combo) + "GPGPG"

    if st.button("Generate Vaccines"):
        if "final" in st.session_state:
            ep = st.session_state["final"]["Epitope"].tolist()
            combos = list(itertools.combinations(ep, 3))
            vaccines = [(c, construct(c)) for c in combos]
            st.session_state["vaccines"] = vaccines

            for _, v in vaccines[:5]:
                st.code(v)

    # EVALUATION
    st.header("5. Evaluation")

    def evaluate(seq):
        return {"Length": len(seq), "Score": len(seq)/100}

    if st.button("Evaluate"):
        if "vaccines" in st.session_state:
            results = [evaluate(v) for _, v in st.session_state["vaccines"]]
            eval_df = pd.DataFrame(results)
            st.session_state["eval"] = eval_df
            st.dataframe(eval_df)

    # FINAL
    st.header("6. Final Selection")

    if st.button("Select Best"):
        if "eval" in st.session_state:
            best = st.session_state["eval"].sort_values(by="Score", ascending=False).iloc[0]
            st.success("Best Vaccine Candidate")
            st.write(best)

# =========================
# ADVANCED
# =========================
if page == "Advanced Analytics":
    st.title("Advanced Analysis")

    if "filtered" in st.session_state:
        data = st.session_state["filtered"]

        fig, ax = plt.subplots(figsize=(4,3))
        corr = data.select_dtypes(include='number').corr()
        cax = ax.matshow(corr, cmap="viridis")
        fig.colorbar(cax)

        st.pyplot(fig)

        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        st.download_button("Download Correlation Matrix", buf.getvalue(), "correlation.png")

# =========================
# EXPORT
# =========================
if page == "Export Results":
    st.title("Export Data")

    if "final" in st.session_state:
        csv = st.session_state["final"].to_csv(index=False).encode()
        st.download_button("Download Final Epitopes", csv, "epitopes.csv")

    if "vaccines" in st.session_state:
        vacc_df = pd.DataFrame([v for _, v in st.session_state["vaccines"]], columns=["Sequence"])
        csv2 = vacc_df.to_csv(index=False).encode()
        st.download_button("Download Vaccines", csv2, "vaccines.csv")

# =========================
# ABOUT
# =========================
if page == "About":
    st.title("About Project")

    st.write("""
    This project demonstrates computational vaccine design using immunoinformatics techniques.
    It accelerates vaccine discovery using in-silico filtering and analysis.
    """)

# =========================
# TEAM
# =========================
if page == "Team":
    st.title("Project Team")

    st.write("""
    Developer: Shravani Dhokate  
    Academic Guide: Dr. Kushagra Kashyap  
    Program: MSc Bioinformatics  
    Domain: Immunoinformatics & Computational Biology  
    """)