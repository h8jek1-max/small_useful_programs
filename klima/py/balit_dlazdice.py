#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BALENÍ DLAŽDIC – převede databázi ze stavba_db.py (bloky po 450 buňkách) do stejné struktury jako stará databáze
(slucovac.py): jedna velká binární data + rejstřík JSON, dlaždice TILE x TILE buněk, každá dlaždice jeden gzip blok.
Hotové soubory obslouží stávající PHP (get_agric_post.php, jen s jinými názvy souborů) a dekóduje je JavaScript
v prohlížeči (klima_dekoder.js).

Nepotřebuje netCDF data ani dočasné úložiště (tmp_*), jen soubory ve složce databaze/ z stavba_db.py.
Je rychlé a dá se kdykoli spustit znovu. Výstup: databaze/web/

    python balit_dlazdice.py                 zabalí dlaždice, vytvoří rejstřík, ověří a vypíše chyby
    python balit_dlazdice.py --test-js       navíc zapíše web/test_ref.json pro test dekodéru v Node
    python balit_dlazdice.py --rok-od 2011   první rok dat (výchozí 2011)
"""
import sys, os, json, gzip, base64, math, time, types
import numpy as np
import odhad_indikatoru as oi
import stavba_db as S
from odhad_vq2 import FEAT, IX, NF, POPIS, COUNTS, WARP, PRIORITA, SIG_MIN, K1, K2

VYSTUP = S.VYSTUP
WEB = os.path.join(VYSTUP, "web")
ROK_OD = int(S.arg("--rok-od", 2011))
N_DATA = "celysvet_klima.gz.bin"
N_REG = "registrace_klima.json"
N_KODY = "kody_klima.bin"


def cas(s): return S.cas(s)

def cesta(n): return os.path.join(VYSTUP, n)

def sekce_definice(M):
    """Pevné pořadí sekcí v dlaždici: (jméno, bajtů na buňku, lag filtru). Lag se doplní z dat."""
    NY, NT, nper = M["NY"], M["NT"], M["nper"]
    nb_idx = max(1, math.ceil(math.log2(M["K_ANOM"])))
    s = [("p1", nper * 36), ("p2", nper * 36), ("sc", nper * 12), ("idx_lo", NY * 12)]
    for k in range(nb_idx - 8): s.append((f"idx_hi{k}", -(-(NY * 12) // 8)))
    for f in COUNTS:
        s.append((f"sum_{f}_lo", NY)); s.append((f"sum_{f}_hi", -(-NY // 8)))
    for nm in ("tgq", "txxq", "ls", "fa"): s.append((nm, NY))
    s.append(("r_lo", NT)); s.append(("r_hi", NT))
    return s, nb_idx

def vyber_lagy(tabulky, nazvy):
    """Nejčastější lag filtru pro každou část podle původních bloků."""
    cnt = {n: {} for n in nazvy}
    for sklad, bloky in tabulky.items():
        for blok in bloky:
            for (nm, off, ln, lg) in blok:
                if nm in cnt: cnt[nm][lg] = cnt[nm].get(lg, 0) + 1
    return {n: (max(c, key=c.get) if c else 0) for n, c in cnt.items()}


def sekce_z_dat(d, M, nb_idx):
    """Slovník dekódovaných polí (n buněk) -> {jméno sekce: uint8 (n, m)}."""
    n = d["p1"].shape[0]
    out = {"p1": d["p1"], "p2": d["p2"], "sc": d["sc"]}
    idx = d["idx"].astype(np.int32)
    out["idx_lo"] = (idx & 255).astype(np.uint8)
    for k in range(nb_idx - 8): out[f"idx_hi{k}"] = np.packbits(((idx >> (8 + k)) & 1).astype(np.uint8), axis=1)
    for j, f in enumerate(COUNTS):
        v = d["sums"][:, j, :].astype(np.int32)
        out[f"sum_{f}_lo"] = (v & 255).astype(np.uint8)
        out[f"sum_{f}_hi"] = np.packbits(((v >> 8) & 1).astype(np.uint8), axis=1)
    for nm in ("tgq", "txxq", "ls", "fa"): out[nm] = d[nm]
    zz = ((d["r"] << 1) ^ (d["r"] >> 31)).astype(np.uint32)
    out["r_lo"] = (zz & 255).astype(np.uint8); out["r_hi"] = ((zz >> 8) & 255).astype(np.uint8)
    return out


def cast_sekce(nm):
    if nm.startswith("sum_") or nm in ("tgq", "txxq", "ls", "fa"): return "rocni"
    if nm in ("r_lo", "r_hi"): return "tnn"
    return "model"

CASTI = ("rocni", "model", "tnn")

def zabal_casti(sek, exc, sekce, lagy):
    """Vrací {cast: gzip blob}. Části: rocni (přesné roční hodnoty), model (základ+index+výjimky SDII), tnn (zbytek TNn po dekádách)."""
    n = next(iter(sek.values())).shape[0]
    out = {}
    for cast in CASTI:
        casti = [np.array([n, exc.shape[1] if cast == "model" else 0], "<u4").tobytes()]
        for nm, m in sekce:
            if cast_sekce(nm) != cast: continue
            a = S.filtr(np.ascontiguousarray(sek[nm]), lagy[nm])
            assert a.shape == (n, m), (nm, a.shape, n, m)
            casti.append(np.ascontiguousarray(a).tobytes())
        if cast == "model":
            raw = np.ascontiguousarray(np.stack([(exc[0] & 255), (exc[0] >> 8), (exc[1] & 255), (exc[1] >> 8),
                                                 (exc[2] & 255), (exc[2] >> 8)]).astype(np.uint8))
            casti.append(raw.tobytes())
        out[cast] = gzip.compress(b"".join(casti), compresslevel=9, mtime=0)
    return out


def zabal_dlazdici(sek, exc, sekce, lagy):
    """sek: {jméno: (n, m) uint8}, exc: (3, nrec) uint16 (buňka v dlaždici, dekáda, hodnota float16 jako uint16)."""
    n = next(iter(sek.values())).shape[0]
    casti = [np.array([n, exc.shape[1]], "<u4").tobytes()]
    for nm, m in sekce:
        a = S.filtr(np.ascontiguousarray(sek[nm]), lagy[nm])
        assert a.shape == (n, m), (nm, a.shape, n, m)
        casti.append(np.ascontiguousarray(a).tobytes())
    raw = np.ascontiguousarray(np.stack([(exc[0] & 255), (exc[0] >> 8), (exc[1] & 255), (exc[1] >> 8),
                                         (exc[2] & 255), (exc[2] >> 8)]).astype(np.uint8))
    casti.append(raw.tobytes())
    return gzip.compress(b"".join(casti), compresslevel=9, mtime=0)


def rozbal_dlazdici(blobs, reg):
    """Opak zabal_casti: {cast: gzip blob} -> slovník polí pro dekódování (stejné klíče jako S.nacti_blok)."""
    sek = {}; ex = None; n = None
    for cast in CASTI:
        b = gzip.decompress(blobs[cast])
        n, nrec = [int(x) for x in np.frombuffer(b[:8], "<u4")]
        pos = 8
        for s in reg["dekoder"]["casti"][cast]:
            m = s["m"]; a = np.frombuffer(b, np.uint8, n * m, pos).reshape(n, m); pos += n * m
            sek[s["jmeno"]] = S.inv_filtr(a, s["lag"])
        if cast == "model": ex = np.frombuffer(b, np.uint8, 6 * nrec, pos).reshape(6, nrec).astype(np.uint16)
    M = reg["dekoder"]; NY = M["NY"]; nb_idx = M["nb_idx"]
    d = {"p1": sek["p1"], "p2": sek["p2"], "sc": sek["sc"]}
    idx = sek["idx_lo"].astype(np.int32)
    for k in range(nb_idx - 8):
        idx |= np.unpackbits(sek[f"idx_hi{k}"], axis=1)[:, :NY * 12].astype(np.int32) << (8 + k)
    d["idx"] = idx.astype(np.uint16)
    sums = np.empty((n, 8, NY), np.uint16)
    for j, f in enumerate(COUNTS):
        sums[:, j, :] = sek[f"sum_{f}_lo"].astype(np.int32) | (np.unpackbits(sek[f"sum_{f}_hi"], axis=1)[:, :NY].astype(np.int32) << 8)
    d["sums"] = sums
    for nm in ("tgq", "txxq", "ls", "fa"): d[nm] = sek[nm]
    zz = (sek["r_lo"].astype(np.uint32) | (sek["r_hi"].astype(np.uint32) << 8)).astype(np.int64)
    d["r"] = ((zz >> 1) ^ -(zz & 1)).astype(np.int32)
    d["exc"] = np.stack([ex[0] | (ex[1] << 8), ex[2] | (ex[3] << 8), ex[4] | (ex[5] << 8)], axis=0)
    return d


def main():
    t0 = time.time()
    os.makedirs(WEB, exist_ok=True)
    with open(cesta("meta.json"), encoding="utf-8") as f: M = json.load(f)
    with open(cesta("stav.json"), encoding="utf-8") as f: stav = json.load(f)
    with open(cesta("tabulky.json"), encoding="utf-8") as f: tabulky = json.load(f)
    if not stav.get("D_hotovo"): sys.exit("Databáze ještě není dokončená (stavba_db.py).")
    CB = M["blok_bunek"]
    G = types.SimpleNamespace(nper=M["nper"], ny=M["ny"], NY=M["NY"], NT=M["NT"], NC=M["NC"])
    rc = np.frombuffer(open(cesta("bunky.bin"), "rb").read(), "<u2").reshape(-1, 2)
    rows = rc[:, 0].astype(np.int64); cols = rc[:, 1].astype(np.int64)
    lat = np.array(M["lat"]); lon = np.array(M["lon"])
    jih_all = lat[rows] < 0
    T = oi.TILE
    ncols = len(lon); radku = int(np.sum(lat >= oi.LAT_MIN - 1e-6))
    PX, PY = math.ceil(ncols / T), math.ceil(radku / T)
    tid = (rows // T) * PX + (cols // T)
    if np.any(np.diff(tid) < 0): sys.exit("Chyba: pořadí buněk neodpovídá dlaždicím (změnil se oi.TILE?)")
    pocty = np.bincount(tid, minlength=PX * PY)
    starty = np.concatenate([[0], np.cumsum(pocty)])
    p(f"mřížka {ncols} x {len(lat)}, ukládá se {radku} řad, dlaždice {T}x{T} -> {PX} x {PY}; buněk {G.NC}; neprázdných dlaždic {int((pocty > 0).sum())}")

    sekce, nb_idx = sekce_definice({"NY": G.NY, "NT": G.NT, "nper": G.nper, "K_ANOM": M["K_ANOM"]})
    lagy = vyber_lagy(tabulky, [n for n, _ in sekce])
    sekce_reg = [{"jmeno": n, "m": m, "lag": int(lagy[n])} for n, m in sekce]

    cache = {}
    def blok_data(b):
        if b not in cache:
            n = min((b + 1) * CB, G.NC) - b * CB
            cache[b] = S.nacti_blok(G, stav, b, n)
            for k in [k for k in cache if k < b - 2]: del cache[k]
        return cache[b]

    reg_dl = {c: [] for c in CASTI}; off = 0; t1 = time.time(); ntile = PX * PY
    with open(os.path.join(WEB, N_DATA), "wb") as fo:
        for t in range(ntile):
            a, e = int(starty[t]), int(starty[t + 1])
            if e == a:
                for c in CASTI: reg_dl[c].append({"offset": off, "delka": 0})
                continue
            parts = {}; exc_l = []
            for b in range(a // CB, (e - 1) // CB + 1):
                d = blok_data(b); b0 = b * CB
                lo, hi = max(a, b0) - b0, min(e, b0 + d["p1"].shape[0]) - b0
                sub = {k: (v[lo:hi] if k != "exc" else None) for k, v in d.items()}
                for nm, arr in sekce_z_dat({k: v for k, v in sub.items() if k != "exc"}, M, nb_idx).items():
                    parts.setdefault(nm, []).append(arr)
                ex = d["exc"]; gc = ex[0].astype(np.int64) + b0
                m = (gc >= a) & (gc < e)
                exc_l.append(np.stack([(gc[m] - a).astype(np.uint16), ex[1][m], ex[2][m]]))
            sek = {nm: np.concatenate(v, axis=0) for nm, v in parts.items()}
            exc = np.concatenate(exc_l, axis=1)
            for c, blob in zabal_casti(sek, exc, sekce, lagy).items():
                fo.write(blob); reg_dl[c].append({"offset": off, "delka": len(blob)}); off += len(blob)
            if t % 5 == 0 or t == ntile - 1: S.prog(t + 1, ntile, 0, t1, "dlaždice")
    S.prog(ntile, ntile, 0, t1, "dlaždice")

    maska = np.zeros((radku, ncols), bool); maska[rows, cols] = True
    reg = {
        "verze": 3, "popis": "klimatické indikátory po dekádách (bezExt + TNn po dekádách + výjimky SDII)",
        "mrizka": {"ncols": ncols, "radku": radku, "lon0": float(lon[0]), "lat0": float(lat[0]),
                   "dx": float(abs(lon[1] - lon[0])), "dy": float(abs(lat[1] - lat[0])), "tile": T,
                   "dlazdic_x": PX, "dlazdic_y": PY},
        "pocet_kroku": G.NT, "rok_od": ROK_OD, "pocet_bunek": G.NC,
        "maska_pevniny_b64": base64.b64encode(np.packbits(maska.ravel()).tobytes()).decode("ascii"),
        "soubor_kody": N_KODY,
        "dekoder": {"FEAT": FEAT, "COUNTS": COUNTS, "NF": NF, "nper": G.nper, "ny": G.ny, "NY": G.NY, "NT": G.NT,
                    "K_ANOM": M["K_ANOM"], "K1": K1, "K2": K2, "nb_idx": nb_idx, "krok_tnn": M["krok_tnn"],
                    "strop_tnn": M["strop_tnn"], "sdii_prah": M["sdii_prah"], "SIG_MIN": SIG_MIN,
                    "teplota_krok": 0.5, "teplota_od": -75.0,
                    "warp": {k: {"hranice": [float(x) for x in h], "sklony": [float(x) for x in s]} for k, (h, s) in WARP.items()},
                    "priorita": {k: float(v) for k, v in PRIORITA.items()},
                    "casti": {c: [x for x in sekce_reg if cast_sekce(x["jmeno"]) == c] for c in CASTI}},
        "casti": reg_dl,
    }
    with open(os.path.join(WEB, N_REG), "w", encoding="utf-8") as f: json.dump(reg, f, ensure_ascii=False, separators=(",", ":"))
    import shutil
    shutil.copyfile(cesta("kody.bin"), os.path.join(WEB, N_KODY))
    mb = lambda n: os.path.getsize(os.path.join(WEB, n)) / 1e6
    for c in CASTI:
        p(f"  část {c:<7}{sum(x['delka'] for x in reg_dl[c]) / 1e6:9.1f} MB")
    p(f"\nHOTOVO ve složce {WEB}:\n  {N_DATA:<24}{mb(N_DATA):9.1f} MB\n  {N_REG:<24}{mb(N_REG):9.2f} MB\n  {N_KODY:<24}{mb(N_KODY):9.2f} MB (stahuje se jednou)")
    overit(reg, G, jih_all, starty)
    p(f"  čas {cas(time.time() - t0)}")


def p(s=""): print(s)


def nacti_dlazdici(reg, t):
    out = {}
    with open(os.path.join(WEB, N_DATA), "rb") as f:
        for c in CASTI:
            e = reg["casti"][c][t]; f.seek(e["offset"]); out[c] = f.read(e["delka"])
    return out

def overit(reg, G, jih_all, starty):
    p("\nOVĚŘENÍ z hotových dlaždic (jen čtením rejstříku a souboru dat, jako PHP + prohlížeč):")
    KOD = S.nacti_kody()
    z = np.load(cesta("overeni.npz")); buny = [int(c) for c in z["bunka"]]; Tr = z["data"].astype(np.float32)
    T = reg["mrizka"]["tile"]; PX = reg["mrizka"]["dlazdic_x"]
    with open(cesta("bunky.bin"), "rb") as f: rc = np.frombuffer(f.read(), "<u2").reshape(-1, 2)
    E = []; ref = []
    for k, c in enumerate(buny):
        t = (int(rc[c, 0]) // T) * PX + int(rc[c, 1]) // T
        d = rozbal_dlazdici(nacti_dlazdici(reg, t), reg); r = c - int(starty[t])
        dr = {a: v[r:r + 1] for a, v in d.items() if a != "exc"}
        ex = d["exc"]; m = ex[0] == r
        dr["exc"] = np.stack([ex[0][m] * 0, ex[1][m], ex[2][m]])
        R = S.dekodovat(G, KOD, dr, jih_all[c:c + 1])
        E.append(R[0] - Tr[k]); ref.append((c, t, r, R[0]))
    E = np.stack(E)
    p(f"  kontrolních buněk: {len(buny)}")
    p(f"  {'indikátor':<11}{'RMSE':>8}{'MAE':>8}{'bias':>8}{'P95':>8}{'P99':>8}{'max':>9}")
    for j, f in enumerate(FEAT):
        e = E[..., j].ravel(); a = np.abs(e)
        p(f"  {POPIS[f]:<11}{np.sqrt((e * e).mean()):>8.2f}{a.mean():>8.2f}{e.mean():>+8.2f}{np.percentile(a, 95):>8.2f}{np.percentile(a, 99):>8.2f}{a.max():>9.2f}")
    if "--test-js" in sys.argv:
        out = []
        for (c, t, r, R) in ref[:3]:
            out.append({"bunka": c, "dlazdice": t, "radek": r, "jih": bool(jih_all[c]), "R": np.round(R.astype(np.float64), 4).ravel().tolist()})
        with open(os.path.join(WEB, "test_ref.json"), "w") as f: json.dump(out, f)
        p("  zapsáno web/test_ref.json (referenční výsledek z Pythonu pro test dekodéru v JavaScriptu)")


if __name__ == "__main__":
    main()
