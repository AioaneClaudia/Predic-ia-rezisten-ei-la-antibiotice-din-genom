"""Descarca contigurile fiecarui genom din BV-BRC (API HTTPS), calculeaza schita de k-meri
si o salveaza in data/sketches/. Genomele (FASTA) NU se pastreaza pe disc: doar schitele
(circa 40 KB per genom). Scriptul se poate opri si relua: genomele deja procesate se sar.

Folosire:  python src/03_descarca_si_schiteaza.py --proba 3     (test pe 3 genome)
           python src/03_descarca_si_schiteaza.py               (toate genomele)
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from kmeri import K, SCALED, lista_genome, schita_genom

API = "https://www.bv-brc.org/api"
DEST = Path("data/sketches")
LOG = Path("data/processed/schite_log.csv")


def _cere(gid, valoare):
    url = (f"{API}/genome_sequence/?eq(genome_id,{valoare})"
           f"&select(genome_id,accession,sequence)&limit(20000)")
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=300)
    r.raise_for_status()
    return [x for x in r.json() if str(x.get("genome_id")) == gid]   # pastram doar genomul cerut


def descarca_contiguri(gid, incercari=4):
    """Contigurile unui genom, ca lista de inregistrari JSON (cu campul 'sequence')."""
    for incercare in range(incercari):
        try:
            inreg = _cere(gid, gid)
            if not inreg:                      # ID-ul cu zerouri finale poate fi citit ca numar
                try:
                    inreg = _cere(gid, f"string:{gid}")
                except requests.HTTPError:
                    inreg = []
            return inreg
        except (requests.RequestException, ValueError):
            time.sleep(2 * (incercare + 1))
    raise RuntimeError(f"Nu am putut descarca {gid}")


def proceseaza(gid):
    DEST.mkdir(parents=True, exist_ok=True)
    f = DEST / f"{gid}.npy"
    if f.exists():
        return gid, "exista", (0, 0, 0)
    inreg = descarca_contiguri(gid)
    secvente = [x["sequence"] for x in inreg if x.get("sequence")]
    if not secvente:
        return gid, "fara_secventa", (0, 0, 0)
    schita = schita_genom(secvente, K, SCALED)
    np.save(f, schita)
    return gid, "ok", (len(secvente), sum(len(s) for s in secvente), len(schita))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--proba", type=int, default=0, help="proceseaza doar primele N genome")
    ap.add_argument("--fire", type=int, default=4, help="numar de descarcari in paralel")
    args = ap.parse_args()

    DEST.mkdir(parents=True, exist_ok=True)
    ids = lista_genome()
    if args.proba:
        ids = ids[:args.proba]
    print(f"{len(ids)} genome de procesat ({sum((DEST / f'{g}.npy').exists() for g in ids)} deja gata)")

    randuri, t0 = [], time.time()
    with ThreadPoolExecutor(max_workers=args.fire) as pool:
        viitor = {pool.submit(proceseaza, g): g for g in ids}
        for i, fut in enumerate(as_completed(viitor), 1):
            g = viitor[fut]
            try:
                gid, stare, (nc, lung, nh) = fut.result()
            except Exception as e:                 # un esec nu opreste restul
                gid, stare, nc, lung, nh = g, f"eroare: {e}", 0, 0, 0
            randuri.append({"genome_id": gid, "stare": stare, "contiguri": nc,
                            "lungime": lung, "hash_uri": nh})
            if i % 50 == 0 or i == len(ids):
                print(f"  {i}/{len(ids)} ({time.time() - t0:.0f} s)")

    log = pd.DataFrame(randuri)
    if args.proba:
        print(log)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log.to_csv(LOG, index=False)
    print(log["stare"].value_counts())
    ok = log[log["stare"] == "ok"]
    if len(ok):
        print("lungime medie:", int(ok["lungime"].mean()), "| hash-uri medii:", int(ok["hash_uri"].mean()))