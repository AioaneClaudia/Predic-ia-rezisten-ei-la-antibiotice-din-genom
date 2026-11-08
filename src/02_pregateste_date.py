"""Curata datele pentru proof of concept: ciprofloxacin, E. coli, Norvegia, 2002-2017.

Intrare : data/raw/amr_ecoli.tsv, data/raw/genome_meta_ecoli.tsv
Iesire  : data/processed/cipro_norvegia.tsv (un rand per genom, y = 1 pentru rezistent)
"""
from pathlib import Path

import pandas as pd

RAW = Path("data/raw")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)

amr = pd.read_csv(RAW / "amr_ecoli.tsv", sep="\t", dtype=str)
meta = pd.read_csv(RAW / "genome_meta_ecoli.tsv", sep="\t", dtype=str)

# 1. Doar masuratori de laborator (fara predictiile XGBoost ale BV-BRC), pentru ciprofloxacin
lab = amr[amr["computational_method"].isna() & (amr["antibiotic"] == "ciprofloxacin")]
lab = lab.dropna(subset=["resistant_phenotype"])


# 2. O singura eticheta per genom; daca inregistrarile se contrazic, marcam "Conflict"
def combina(serie):
    return serie.iloc[0] if serie.nunique() == 1 else "Conflict"


etichete = lab.groupby("genome_id")["resistant_phenotype"].agg(combina).reset_index()
etichete = etichete[etichete["resistant_phenotype"].isin(["Resistant", "Susceptible"])]

# 3. Metadate: an, tara, calitate, ST si clustere cgMLST
coloane = ["genome_id", "collection_year", "isolation_country", "genome_quality",
           "mlst", "cgmlst_hc0", "cgmlst_hc2", "cgmlst_hc5", "cgmlst_hc10",
           "cgmlst_hc20", "cgmlst_hc50", "cgmlst_hc100"]
df = etichete.merge(meta[coloane], on="genome_id", how="left")
df["collection_year"] = pd.to_numeric(df["collection_year"], errors="coerce")

# 4. Filtre: Norvegia, 2002-2017, calitate buna
df = df[(df["isolation_country"] == "Norway")
        & df["collection_year"].between(2002, 2017)
        & (df["genome_quality"] == "Good")].copy()
df["y"] = (df["resistant_phenotype"] == "Resistant").astype(int)

df.to_csv(OUT / "cipro_norvegia.tsv", sep="\t", index=False)
print(f"{len(df)} genome, prevalenta rezistenta = {df['y'].mean():.3f}")
print(df.groupby(df["collection_year"] <= 2012)["y"].agg(["sum", "count"]))