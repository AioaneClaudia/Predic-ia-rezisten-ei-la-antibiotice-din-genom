"""Construieste matricea binara genom x hash din schitele salvate in data/sketches/.

Iesire: data/processed/kmeri_X.npz (matrice rara) si data/processed/kmeri_ids.txt (ordinea genomelor).
"""
from pathlib import Path

import pandas as pd
from scipy import sparse

from kmeri import lista_genome, matrice_binara

DIR_SCHITE = Path("data/sketches")

ids = lista_genome()
ok = [g for g in ids if (DIR_SCHITE / f"{g}.npy").exists()]
print(f"{len(ok)} din {len(ids)} genome au schita")
if len(ok) < len(ids):
    print("lipsesc (primele 10):", [g for g in ids if g not in set(ok)][:10])

X, hash_unice = matrice_binara(ok, DIR_SCHITE)
sparse.save_npz("data/processed/kmeri_X.npz", X)
pd.Series(ok).to_csv("data/processed/kmeri_ids.txt", index=False, header=False)
mb = (X.data.nbytes + X.indices.nbytes + X.indptr.nbytes) / 1e6
print(f"matrice: {X.shape[0]} genome x {X.shape[1]} hash-uri distincte, {X.nnz} valori nenule (~{mb:.0f} MB)")