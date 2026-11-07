"""Descarca fenotipurile AMR (E. coli) si metadatele genomurilor din BV-BRC prin API HTTPS.

Folosim API-ul, nu FTP-ul, deoarece certificatul serverului FTP nu se potriveste cu numele lui.
Rezultate in data/raw/ (ignorate de git):
  amr_ecoli.tsv          fenotipuri de rezistenta pentru E. coli
  genome_meta_ecoli.tsv  metadate (an, tara etc.) pentru aceleasi genome
"""
import io
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://www.bv-brc.org/api"
ANTIBIOTICE = ["ciprofloxacin", "ampicillin"]
DESTINATIE = Path("data/raw")
DESTINATIE.mkdir(parents=True, exist_ok=True)
SESIUNE = requests.Session()


def ia_tsv(colectie, interogare):
    """Trimite o interogare RQL si intoarce rezultatul ca DataFrame (TSV)."""
    url = f"{API}/{colectie}/?{interogare}"
    raspuns = SESIUNE.get(url, headers={"Accept": "text/tsv"}, timeout=120)
    raspuns.raise_for_status()
    if not raspuns.text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(raspuns.text), sep="\t", dtype=str)


def ia_toate(colectie, interogare, pagina=10000):
    """Parcurge rezultatele pagina cu pagina."""
    bucati, start = [], 0
    while True:
        df = ia_tsv(colectie, f"{interogare}&sort(%2Bid)&limit({pagina},{start})")
        bucati.append(df)
        print(f"  {colectie}: {start + len(df)} randuri")
        if len(df) < pagina:
            break
        start += pagina
    return pd.concat(bucati, ignore_index=True)


# 1. Fenotipuri AMR
CAMPURI_AMR = ["genome_id", "genome_name", "taxon_id", "antibiotic", "resistant_phenotype",
               "measurement", "measurement_value", "measurement_unit",
               "laboratory_typing_method", "computational_method"]
tinta_amr = DESTINATIE / "amr_ecoli.tsv"
if tinta_amr.exists():
    amr = pd.read_csv(tinta_amr, sep="\t", dtype=str)
else:
    bucati = []
    for ab in ANTIBIOTICE:
        print(f"Descarc fenotipuri pentru {ab} ...")
        bucati.append(ia_toate("genome_amr",
                               f"eq(antibiotic,{ab})&select({','.join(CAMPURI_AMR)})"))
    amr = pd.concat(bucati, ignore_index=True)
    amr = amr[amr["genome_name"].str.startswith("Escherichia coli", na=False)]
    amr.to_csv(tinta_amr, sep="\t", index=False)
print(f"{len(amr)} randuri AMR, {amr['genome_id'].nunique()} genome E. coli")

# 2. Metadatele genomurilor (JSON, ca sa pastram toate coloanele disponibile)
tinta_meta = DESTINATIE / "genome_meta_ecoli.tsv"
if not tinta_meta.exists():
    ids = sorted(amr["genome_id"].unique())
    randuri = []
    for i in range(0, len(ids), 100):
        lot = ids[i:i + 100]
        url = f"{API}/genome/?in(genome_id,({','.join(lot)}))&limit({len(lot)})"
        r = SESIUNE.get(url, headers={"Accept": "application/json"}, timeout=120)
        r.raise_for_status()
        randuri.extend(r.json())
        time.sleep(0.2)
        print(f"  genome: {min(i + 100, len(ids))}/{len(ids)} id-uri cerute")
    meta = pd.DataFrame(randuri)
    meta.to_csv(tinta_meta, sep="\t", index=False)
    print(f"{len(meta)} genome primite din {len(ids)} cerute")
    print("Coloane cu 'year' sau 'country':",
          [c for c in meta.columns if "year" in c or "country" in c])