import streamlit as st
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import io

try:
    import py3Dmol
    from stmol import showmol
    PDB_AVAILABLE = True
except:
    PDB_AVAILABLE = False

st.set_page_config(layout="wide")


st.markdown("""
<style>
body {
    background-color: #F4F6F7;
}
[data-testid="stAppViewContainer"] {
    background: linear-gradient(135deg, #e3f2fd, #fce4ec);
}
.stButton>button {
    border-radius: 10px;
    background: linear-gradient(90deg, #2E86C1, #48C9B0);
    color: white;
    font-weight: 600;
    border: none;
}
h1, h2, h3 {
    color: #1B4F72;
}
.stDataFrame {
    background-color: white;
    border-radius: 10px;
}
</style>
""", unsafe_allow_html=True)


plt.rcParams.update({
    'figure.autolayout': True,
    'font.size': 8
})

@st.cache_data
def load_data():
    return pd.read_csv("nipah_dataset1.csv")

df = load_data()


st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", [
    "Home","Tool","Vaccine Pipeline",
    "Export Results",
    "Team"
])

if page == "Home":
    st.title("Nipah Virus Vaccine Design Platform")
    st.markdown("### Computational Immunoinformatics-Based Vaccine Development System")
    st.markdown("---")

    st.subheader("Overview")
    st.info("This platform is designed for in-silico multi-epitope vaccine design targeting the Nipah virus.")

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


if page == "Tool":
    st.title("Tool Description")
    st.markdown("---")

    st.subheader("Overview")
    st.write("Streamlit-based platform for multi-epitope vaccine design using immunoinformatics.")

    st.subheader("Methodology")
    st.write("""
    - Protein Selection  
    - Antigenicity Screening  
    - Epitope Classification  
    - Filtering  
    - Vaccine Construction  
    """)

    st.subheader("Key Parameters")
    st.write("""
    - Antigenicity  
    - Binding Affinity  
    - Toxicity  
    - Allergenicity  
    - Conservation  
    """)

    st.success("Output: Ranked vaccine candidates")


if page == "Vaccine Pipeline":

    st.title("Vaccine Design Pipeline")

    if st.checkbox("Show Dataset"):
        st.dataframe(df)

    # TARGET
    st.header("1. Target Selection")
    protein = st.selectbox("Select Protein", df["Protein"].unique())
    df_protein = df[df["Protein"] == protein]
    st.dataframe(df_protein)

    # GRAPH 1 (COLORED BAR)
    if "Subcellular_Location" in df_protein.columns:
        fig, ax = plt.subplots(figsize=(4,2.5))
        df_protein["Subcellular_Location"].value_counts().plot(
            kind='bar',
            ax=ax,
            color=['#5DADE2','#48C9B0','#AF7AC5','#F5B041']
        )

        ax.set_xlabel("Location")
        ax.set_ylabel("Count")
        ax.set_title("Subcellular Localization")

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
                    st.warning("3D not supported")
                    st.markdown(f"https://www.rcsb.org/structure/{selected_pdb}")

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

    # GRAPH 2 (COLORED HIST)
    if "filtered" in st.session_state:
        fig, ax = plt.subplots(figsize=(4,2.5))
        ax.hist(
            st.session_state["filtered"]["Antigenicity"],
            color="#48C9B0",
            edgecolor="black"
        )

        ax.set_xlabel("Antigenicity")
        ax.set_ylabel("Count")
        ax.set_title("Antigenicity Distribution")

        st.pyplot(fig)

        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        st.download_button("Download Histogram", buf.getvalue(), "antigenicity.png")

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


if page == "Export Results":
    st.title("Export Data")

    if "final" in st.session_state:
        csv = st.session_state["final"].to_csv(index=False).encode()
        st.download_button("Download Epitopes", csv, "epitopes.csv")

    if "vaccines" in st.session_state:
        vacc_df = pd.DataFrame([v for _, v in st.session_state["vaccines"]], columns=["Sequence"])
        csv2 = vacc_df.to_csv(index=False).encode()
        st.download_button("Download Vaccines", csv2, "vaccines.csv")


if page == "Team":
    st.title(" Project Team")

    st.write("""
    Developer: Shravani Vinod Dhokate  
    Academic Guide: Dr. Kushagra Kashyap  
    Program: MSc Bioinformatics  
    Domain: Immunoinformatics & Computational Biology  
    """)
