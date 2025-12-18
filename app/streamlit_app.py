import requests
import streamlit as st

API_BASE = "http://127.0.0.1:8000"

st.set_page_config(page_title="ChemiMLOps Demo", layout="centered")
st.title("🧪 ChemiMLOps — ADMET + Similarity Search")
st.caption("Enter a SMILES string → get lipophilicity prediction + top-K similar molecules.")

smiles = st.text_input("SMILES", value="c1ccccc1")
k = st.slider("Top-K similar molecules", min_value=1, max_value=20, value=5)

col1, col2 = st.columns(2)

with col1:
    if st.button("Predict Lipophilicity"):
        try:
            r = requests.post(f"{API_BASE}/predict/lipophilicity", json={"smiles": smiles}, timeout=15)
            if r.status_code != 200:
                st.error(r.json().get("detail", "Prediction failed"))
            else:
                pred = r.json()["lipophilicity_pred"]
                st.success(f"Predicted lipophilicity: **{pred:.4f}**")
        except Exception as e:
            st.error(f"API error: {e}")

with col2:
    if st.button("Find Similar Molecules"):
        try:
            r = requests.get(f"{API_BASE}/similarity/topk", params={"smiles": smiles, "k": k}, timeout=30)
            if r.status_code != 200:
                st.error(r.json().get("detail", "Similarity search failed"))
            else:
                data = r.json()
                st.subheader("Top Similar Molecules")
                st.write(f"Query: `{data['query']}`  |  k={data['k']}")
                st.dataframe(data["results"], use_container_width=True)
        except Exception as e:
            st.error(f"API error: {e}")

st.divider()
st.info("Tip: Keep the FastAPI server running (uvicorn) while using this app.")
