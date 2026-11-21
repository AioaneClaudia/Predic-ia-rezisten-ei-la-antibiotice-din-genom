"""Functii pentru evaluarea modelelor de predictie a rezistentei la antibiotice.

Continut:
  - incarcarea datelor curate
  - impartiri: aleatorie, temporala, pe clustere
  - model doar-linie (baseline), model pe k-meri si distanta pana la cel mai apropiat izolat
  - metrici (AUROC si castigul AUPRC fata de prevalenta)
  - controlul cu distanta potrivita
  - tabelul pentru regresia la nivel de izolat
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

NIVELURI = [0, 2, 5, 10, 20, 50, 100]
MODELE = {"ST": "mlst", **{f"HC{h}": f"cgmlst_hc{h}" for h in NIVELURI}}
ETICHETE_DISTANTA = ["<=HC10", "HC20", "HC50", "HC100", ">HC100"]


# ---------------------------------------------------------------- date
def incarca_date(cale="data/processed/cipro_norvegia.tsv"):
    """Citeste fisierul curat creat de src/02_pregateste_date.py."""
    df = pd.read_csv(cale, sep="\t", dtype=str)
    df["collection_year"] = df["collection_year"].astype(float).astype(int)
    df["y"] = df["y"].astype(int)
    return df.reset_index(drop=True)


# ------------------------------------------------------------ impartiri
def split_temporal(df, an_taiere=2012):
    """Antrenare: ani <= an_taiere. Test: ani > an_taiere."""
    return df[df["collection_year"] <= an_taiere], df[df["collection_year"] > an_taiere]


def split_aleatoriu(df, n_train, n_test, seed):
    """Impartire aleatorie stratificata, cu dimensiuni fixate."""
    return train_test_split(df, train_size=n_train, test_size=n_test,
                            stratify=df["y"], random_state=seed)


def split_pe_clustere(df, nivel, frac_test, seed):
    """Clusterele HC ajung intregi fie in antrenare, fie in test (nicio ruda la nivelul dat)."""
    rng = np.random.default_rng(seed)
    col = f"cgmlst_hc{nivel}"
    marimi = df[col].value_counts()
    clustere = marimi.index.to_numpy().copy()
    rng.shuffle(clustere)
    tinta = frac_test * len(df)
    in_test, n = [], 0
    for c in clustere:
        if n >= tinta:
            break
        in_test.append(c)
        n += marimi[c]
    masca = df[col].isin(in_test)
    return df[~masca], df[masca]


# ----------------------------------------------------- model si distanta
def model_linie(train, test, col, alpha=1.0):
    """Probabilitate de rezistenta = proportia de rezistente din antrenare in acelasi
    cluster/ST. Se trage spre prevalenta globala (alpha); un cluster nevazut primeste prevalenta."""
    p0 = train["y"].mean()
    stat = train.groupby(col)["y"].agg(["sum", "count"])
    s = test[col].map(stat["sum"]).fillna(0)
    n = test[col].map(stat["count"]).fillna(0)
    return ((s + alpha * p0) / (n + alpha)).values


def distanta_minima(train, test):
    """Cel mai fin nivel HC la care izolatul de test are o ruda in antrenare (inf = niciun nivel)."""
    d = pd.Series(np.inf, index=test.index)
    for h in sorted(NIVELURI, reverse=True):
        col = f"cgmlst_hc{h}"
        d[test[col].isin(set(train[col]))] = h
    return d


def bin_distanta(d):
    """Grupeaza distantele in 5 categorii: <=HC10, HC20, HC50, HC100, >HC100."""
    return pd.cut(d, bins=[-1, 10, 20, 50, 100, np.inf], labels=ETICHETE_DISTANTA)


# ---------------------------------------------------- fabrici de modele
# Un "model" este o functie  model_fn(train, test1, test2, ...)  care antreneaza o singura data
# pe `train` si intoarce cate un vector de probabilitati pentru fiecare test dat.
def fabrica_linie(nume_model):
    """Modelul doar-linie ("ST", "HC20", "HC50", "HC100" etc.)."""
    col = MODELE[nume_model]

    def model_fn(train, *teste):
        return tuple(model_linie(train, te, col) for te in teste)
    return model_fn


def fabrica_kmeri(X, id2row, C=0.1, min_df=5):
    """Regresie logistica pe prezenta/absenta k-merilor (matrice rara X, un rand per genom).
    Trasaturile prezente la mai putin de `min_df` genome DIN ANTRENARE se elimina."""
    def model_fn(train, *teste):
        Xtr = X[train["genome_id"].map(id2row).values]
        pastreaza = np.asarray(Xtr.sum(axis=0)).ravel() >= min_df
        Xtr = Xtr[:, pastreaza]
        clf = LogisticRegression(C=C, solver="liblinear", max_iter=200)
        clf.fit(Xtr, train["y"].values)
        return tuple(clf.predict_proba(X[te["genome_id"].map(id2row).values][:, pastreaza])[:, 1]
                     for te in teste)
    return model_fn


# --------------------------------------------------------------- metrici
def metrici(y, scor, rng=None):
    """AUROC si castigul AUPRC (AUPRC minus prevalenta). Daca rng e dat, se face un
    esantion bootstrap al izolatelor de test inainte de calcul."""
    y = np.asarray(y)
    scor = np.asarray(scor)
    if rng is not None:
        idx = rng.integers(0, len(y), len(y))
        y, scor = y[idx], scor[idx]
    if y.min() == y.max():
        return {"AUROC": np.nan, "AUPRC_castig": np.nan}
    return {"AUROC": roc_auc_score(y, scor),
            "AUPRC_castig": average_precision_score(y, scor) - y.mean()}


def evalueaza(train, test):
    """Evalueaza toate modelele doar-linie pe un test; intoarce un tabel."""
    rez = {nume: metrici(test["y"].values, model_linie(train, test, col))
           for nume, col in MODELE.items()}
    return pd.DataFrame(rez).T


# ---------------------------------------------- controlul cu distanta potrivita
def control_distanta(df, an_taiere=2012, frac_train=0.65, nume_model="HC50",
                     seed=0, bootstrap=True, n_min=50, returneaza_seturi=False,
                     model_fn=None, potriveste_prevalenta=False):
    """Un experiment: ACEEASI izolate de antrenare (din perioada veche), doua teste:
      - testul temporal (perioada noua);
      - un control din perioada veche, esantionat astfel incat distributia distantei
        pana la antrenare sa fie aceeasi ca la testul temporal.
    Cu potriveste_prevalenta=True se potrivesc si rezistentele in fiecare categorie de
    distanta, deci controlul are si aceeasi prevalenta ca testul.
    Daca performanta e aceeasi, scaderea in timp se explica prin distanta."""
    rng = np.random.default_rng(seed)
    vechi = df[df["collection_year"] <= an_taiere]
    nou = df[df["collection_year"] > an_taiere]
    tr, pool = train_test_split(vechi, train_size=frac_train,
                                stratify=vechi["y"], random_state=seed)

    b_nou = bin_distanta(distanta_minima(tr, nou)).astype(str)
    b_pool = bin_distanta(distanta_minima(tr, pool)).astype(str)
    if potriveste_prevalenta:
        k_nou = pd.Series(b_nou.values + "|" + nou["y"].astype(str).values)
        k_pool = pd.Series(b_pool.values + "|" + pool["y"].astype(str).values)
    else:
        k_nou, k_pool = pd.Series(b_nou.values), pd.Series(b_pool.values)
    p_nou = k_nou.value_counts(normalize=True)
    categorii = list(p_nou.index)
    disponibil = k_pool.value_counts().reindex(categorii, fill_value=0)

    # cel mai mare control la care se pot potrivi toate categoriile
    n_ctrl = int(min(disponibil[c] / p_nou[c] for c in categorii))
    if n_ctrl < n_min:
        return None
    k_pool.index = pool.index
    alesi = []
    for c in categorii:
        candidati = pool.index[(k_pool == c).values]
        k = min(int(round(p_nou[c] * n_ctrl)), len(candidati))
        alesi.extend(rng.choice(candidati, size=k, replace=False))
    ctrl = pool.loc[alesi]

    if model_fn is None:
        model_fn = fabrica_linie(nume_model)
    s_nou, s_ctrl = model_fn(tr, nou, ctrl)          # un singur antrenament, doua teste
    r_rng = rng if bootstrap else None
    m_nou = metrici(nou["y"].values, s_nou, r_rng)
    m_ctrl = metrici(ctrl["y"].values, s_ctrl, r_rng)
    rezultat = {
        "n_antrenare": len(tr), "n_test_temporal": len(nou), "n_control": len(ctrl),
        "prev_temporal": nou["y"].mean(), "prev_control": ctrl["y"].mean(),
        "AUROC_temporal": m_nou["AUROC"], "AUROC_control": m_ctrl["AUROC"],
        "castig_temporal": m_nou["AUPRC_castig"], "castig_control": m_ctrl["AUPRC_castig"],
        "delta_AUROC": m_nou["AUROC"] - m_ctrl["AUROC"],
        "delta_castig": m_nou["AUPRC_castig"] - m_ctrl["AUPRC_castig"],
    }
    if returneaza_seturi:
        return rezultat, tr, nou, ctrl
    return rezultat


# ---------------------------------------------------- regresia la nivel de izolat
def tabel_regresie(df, an_taiere=2012, nume_model="HC50", model_fn=None):
    """Un rand per izolat de test (impartirea temporala), cu eroarea modelului,
    distanta pana la antrenare si decalajul temporal."""
    tr, te = split_temporal(df, an_taiere)
    if model_fn is None:
        model_fn = fabrica_linie(nume_model)
    p = model_fn(tr, te)[0]
    out = te.copy()
    out["p"] = p
    out["rezidual"] = out["y"] - out["p"]          # cu semn: >0 = modelul subestimeaza rezistenta
    out["eroare"] = out["rezidual"] ** 2           # eroare Brier la nivel de izolat
    out["d_min"] = distanta_minima(tr, te).values
    out["d_bin"] = bin_distanta(out["d_min"]).astype(str)
    out["decalaj"] = out["collection_year"] - an_taiere
    return out