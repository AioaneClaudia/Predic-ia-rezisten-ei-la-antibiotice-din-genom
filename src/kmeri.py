"""Schite de k-meri (FracMinHash simplificat).

Dintr-un genom pastram doar k-merii canonici al caror hash e sub un prag, adica
aproximativ 1/SCALED dintre toti k-merii. Doua genome inrudite au multe hash-uri comune.
"""
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

K = 21
SCALED = 1000

# tabel de codificare: A,C,G,T -> 0..3 (litere mici la fel); orice altceva -> 255 (invalid)
_COD = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _COD[ord(_b)] = _i
    _COD[ord(_b.lower())] = _i


def _amestecare(x):
    """Functie de hash pe 64 de biti (splitmix64): imprastie uniform k-merii."""
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def schita_contig(secventa, k=K, scaled=SCALED):
    """Hash-urile pastrate pentru un contig (k-meri canonici = min(secventa, complement invers))."""
    arr = np.frombuffer(secventa.encode("ascii", "ignore"), dtype=np.uint8)
    n = len(arr)
    if n < k:
        return np.empty(0, dtype=np.uint64)
    cod = _COD[arr]
    valid = cod < 4
    m = n - k + 1
    # ferestrele care contin baze invalide (N etc.) se ignora
    inv = np.concatenate(([0], np.cumsum(~valid, dtype=np.int64)))
    fereastra_ok = (inv[k:] - inv[:-k]) == 0
    c0 = np.where(valid, cod, 0).astype(np.uint64)
    comp = np.uint64(3) - c0
    fw = np.zeros(m, dtype=np.uint64)
    rc = np.zeros(m, dtype=np.uint64)
    for j in range(k):
        fw = (fw << np.uint64(2)) | c0[j:j + m]               # k-mer pe catena directa
        rc = rc | (comp[j:j + m] << np.uint64(2 * j))          # complementul invers
    canonic = np.minimum(fw, rc)[fereastra_ok]
    h = _amestecare(canonic)
    prag = np.uint64(0xFFFFFFFFFFFFFFFF // scaled)
    return h[h < prag]


def schita_genom(contiguri, k=K, scaled=SCALED):
    """Schita unui genom = hash-urile unice din toate contigurile lui."""
    bucati = [schita_contig(c, k, scaled) for c in contiguri]
    if not bucati:
        return np.empty(0, dtype=np.uint64)
    return np.unique(np.concatenate(bucati))


def matrice_binara(ids, dir_schite):
    """Matrice rara (genom x hash): 1 daca genomul contine hash-ul. Intoarce (X, hash_uri)."""
    schite = [np.load(Path(dir_schite) / f"{g}.npy") for g in ids]
    lungimi = np.array([len(s) for s in schite])
    toate = np.concatenate(schite)
    hash_unice, inv = np.unique(toate, return_inverse=True)
    randuri = np.repeat(np.arange(len(ids)), lungimi)
    X = sparse.csr_matrix((np.ones(len(inv), dtype=np.int8), (randuri, inv.ravel())),
                          shape=(len(ids), len(hash_unice)))
    return X, hash_unice


def lista_genome():
    """Genomele din toate fisierele data/processed/*_norvegia.tsv (fara dubluri, in ordine)."""
    ids = []
    for fisier in sorted(glob.glob("data/processed/*_norvegia.tsv")):
        ids += pd.read_csv(fisier, sep="\t", dtype=str)["genome_id"].tolist()
    return list(dict.fromkeys(ids))