#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ZMĚŘÍ, KOLIK by se dalo ušetřit na TNn po dekádách (tnn.dat) lepším kódováním, dřív než se cokoli přepisuje.

    python tnn_entropie.py            (ve složce se stavba_db.py a databaze/; čte vzorek bloků, nic nezapisuje)

Spočítá z uložených zbytků TNn (zigzag hodnoty) entropii: (0) bez kontextu = nejlepší možné bitové balení,
(1) podle zbytku v předchozí dekádě, (2) navíc podle téže dekády v předchozím roce. Entropie je dolní mez toho, co umí
aritmetický/Huffmanův kodér s danou znalostí; srovnává se se skutečným tnn.dat (gzip). Pokud je rozdíl malý, přepisovat nemá smysl.
"""
import sys, os, json, math, types
import numpy as np
import stavba_db as S

def entropie(ctx, z, nz):
    n = len(z); h = np.zeros((int(ctx.max()) + 1, nz), np.int64)
    np.add.at(h, (ctx, z), 1)
    tot = h.sum(1, keepdims=True); p = np.where(h > 0, h / np.maximum(tot, 1), 1)
    return float(-(h * np.log2(p)).sum() / n)

def main():
    with open(S.cesta("meta.json"), encoding="utf-8") as f: M = json.load(f)
    with open(S.cesta("stav.json"), encoding="utf-8") as f: stav = json.load(f)
    G = types.SimpleNamespace(nper=M["nper"], ny=M["ny"], NY=M["NY"], NT=M["NT"], NC=M["NC"])
    CB = M["blok_bunek"]; nb = M["nblok"]; vzorek = list(range(0, nb, max(1, nb // 15)))
    Z = []; 
    for i, b in enumerate(vzorek):
        n = min((b + 1) * CB, G.NC) - b * CB
        r = S.nacti_blok(G, stav, b, n)["r"].astype(np.int64)
        Z.append(((r << 1) ^ (r >> 63)).astype(np.int64))
        sys.stdout.write(f"\r  vzorek bloků: {i + 1}/{len(vzorek)}"); sys.stdout.flush()
    print()
    z = np.concatenate(Z, axis=0)                       # (buňky, NT)
    nz = int(min(z.max() + 1, 64)); z = np.minimum(z, nz - 1)
    cur = z[:, 36:]                                     # od druhého roku, aby existoval předchozí rok
    prev = z[:, 35:-1]; rok = z[:, :-36]
    c = lambda a: np.minimum(a, 7)
    h0 = entropie(np.zeros(cur.size, np.int64), cur.ravel(), nz)
    h1 = entropie(c(prev).ravel(), cur.ravel(), nz)
    h2 = entropie((c(prev) * 8 + c(rok)).ravel(), cur.ravel(), nz)
    nhod = G.NC * G.NT
    skut = os.path.getsize(S.cesta("tnn.dat")) / 1e6
    print(f"  hodnot ve vzorku: {cur.size}; nenulových zbytků: {100 * (cur > 0).mean():.1f} %")
    print(f"  skutečný tnn.dat (gzip): {skut:.1f} MB = {skut * 8e6 / nhod:.2f} bitu na hodnotu")
    for nm, h in (("bez kontextu (bitové balení)", h0), ("podle předchozí dekády", h1), ("podle předchozí dekády a téže dekády loni", h2)):
        print(f"  entropie {nm:<42}: {h:.2f} bitu → {h * nhod / 8e6:6.1f} MB")
    print("  (entropie je dolní mez; skutečný kodér ji nedosáhne úplně, režie kontextů a řídkých hodnot ještě něco přidá)")

if __name__ == "__main__":
    main()
