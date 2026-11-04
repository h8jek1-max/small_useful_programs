# slucovac5.py - náhrada slucovac4.py: totéž + chybějící HWSD2 sloupce + odvozené agro-parametry
# Výstup má stejné názvy souborů i strukturu jako dosud (get_soil.php se nemění), jen přibyla pole.
# Potřebuje: HWSD2_SMU.csv, HWSD2_LAYERS.csv, (volitelně) D_ROOT_DEPTH.csv
import csv, json, gzip, os, re, math

CSV_SMU, CSV_LAYERS, CSV_ROOT = "HWSD2_SMU.csv", "HWSD2_LAYERS.csv", "D_ROOT_DEPTH.csv"
VYSTUP_DATA, VYSTUP_INDEX = "pudy_data_komplet.gz.bin", "pudy_index_komplet.json"

def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None

def rd(v, n=2): return None if v is None else round(v, n)

# ---------- pedotransferové funkce ----------
def saxton_rawls(sand, clay, om):
    """Saxton & Rawls (2006). sand, clay v %, om v % hmotnosti. Vrací bod vadnutí, polní kapacitu,
    nasycení (obj. zlomky) a Ksat (mm/h). Platnost: OM do ~8 %, pro rašeliny neplatí."""
    S, C, OM = sand / 100, clay / 100, min(om, 8.0)
    t1500 = -0.024*S + 0.487*C + 0.006*OM + 0.005*S*OM - 0.013*C*OM + 0.068*S*C + 0.031
    wp = t1500 + 0.14*t1500 - 0.02
    t33 = -0.251*S + 0.195*C + 0.011*OM + 0.006*S*OM - 0.027*C*OM + 0.452*S*C + 0.299
    fc = t33 + (1.283*t33**2 - 0.374*t33 - 0.015)
    ts = 0.278*S + 0.034*C + 0.022*OM - 0.018*S*OM - 0.027*C*OM - 0.584*S*C + 0.078
    s33 = ts + (0.636*ts - 0.107)
    sat = fc + s33 - 0.097*S + 0.043
    wp = max(0.01, wp); fc = max(wp + 0.01, fc); sat = max(fc + 0.01, min(sat, 0.7))
    B = (math.log(1500) - math.log(33)) / (math.log(fc) - math.log(wp))
    return wp, fc, sat, 1930 * (sat - fc) ** (3 - 1 / B)

def k_usle(sand, silt, clay, oc):
    """Faktor erodovatelnosti K (EPIC, Williams 1995), SI: t.ha.h/(ha.MJ.mm)."""
    sn = 1 - sand / 100
    f1 = 0.2 + 0.3 * math.exp(-0.0256 * sand * (1 - silt / 100))
    f2 = (silt / (clay + silt)) ** 0.3 if clay + silt > 0 else 1
    f3 = 1 - 0.25 * oc / (oc + math.exp(3.72 - 2.95 * oc))
    f4 = 1 - 0.7 * sn / (sn + math.exp(-5.51 + 22.9 * sn))
    return round(f1 * f2 * f3 * f4 * 0.1317, 3)

def textura(s, si, c):  # USDA trojúhelník
    if si + 1.5*c < 15: return "písčitá"
    if si + 1.5*c >= 15 and si + 2*c < 30: return "hlinitopísčitá"
    if (7 <= c < 20 and s > 52 and si + 2*c >= 30) or (c < 7 and si < 50 and si + 2*c >= 30): return "písčitohlinitá"
    if 7 <= c < 27 and 28 <= si < 50 and s <= 52: return "hlinitá"
    if (si >= 50 and 12 <= c < 27) or (50 <= si < 80 and c < 12): return "prachovitohlinitá"
    if si >= 80 and c < 12: return "prachovitá"
    if 20 <= c < 35 and si < 28 and s > 45: return "jílovitohlinitopísčitá"
    if 27 <= c < 40 and 20 < s <= 45: return "jílovitohlinitá"
    if 27 <= c < 40 and s <= 20: return "prachovitojílovitohlinitá"
    if c >= 35 and s > 45: return "jílovitopísčitá"
    if c >= 40 and si >= 40: return "prachovitojílovitá"
    if c >= 40: return "jílovitá"
    return "hlinitá"

def odvod_vrstvu(l):
    s, c, oc = num(l["sand"]), num(l["clay"]), num(l["oc"])
    if s is None or c is None or oc is None: return {}
    si = num(l["silt"]); si = max(0.0, 100 - s - c) if si is None else si
    top, bot, bulk = num(l["top"]), num(l["bot"]), num(l["bulk"])
    kam = 1 - (num(l["coarse"]) or 0) / 100          # podíl jemnozemě
    wp, fc, sat, ks = saxton_rawls(s, c, oc * 1.724)
    o = {"wp": rd(wp*100, 1), "fc": rd(fc*100, 1), "sat": rd(sat*100, 1), "ksat": rd(ks, 1),
         "tex": textura(s, si, c), "k_usle": k_usle(s, si, c, oc)}
    if top is not None and bot is not None:
        tl = bot - top
        o["awc_mm"] = rd((fc - wp) * tl * 10 * kam, 1)           # mm vody dostupné rostlinám
        if bulk:
            o["soc_tha"] = rd(oc * bulk * tl * kam, 1)           # t C/ha
            n = num(l["total_n"])
            if n: o["n_tha"] = rd(n / 10 * bulk * tl * kam, 2)   # t N/ha
    return o

def odvod_komponentu(k, root_map):
    root = root_map.get(str(k["root_depth"]).strip(), 100.0)
    taw = soc30 = soc100 = 0.0; ks = []
    for key in sorted(k["layers"]):
        l = k["layers"][key]; top, bot = num(l["top"]), num(l["bot"])
        if top is None or bot is None or "awc_mm" not in l: continue
        tl = bot - top
        taw += l["awc_mm"] * max(0, min(bot, root) - top) / tl
        if top < 100: ks.append(l["ksat"])
        if "soc_tha" in l:
            soc30 += l["soc_tha"] * max(0, min(bot, 30) - top) / tl
            soc100 += l["soc_tha"] * max(0, min(bot, 100) - top) / tl
    kmin = min(ks) if ks else None
    hsg = None if kmin is None else "A" if kmin > 144 else "B" if kmin > 36 else "C" if kmin > 3.6 else "D"
    d1 = k["layers"].get("D1", {})
    ph, ec, esp, al = num(d1.get("ph")), num(d1.get("elec_cond")), num(d1.get("esp")), num(d1.get("alum_sat"))
    bulk, coarse, oc = num(d1.get("bulk")), num(d1.get("coarse")), num(d1.get("oc"))
    om = []
    if ph is not None and ph < 5.5: om.append("kyselá půda (pH < 5,5)")
    if ph is not None and ph > 8.0: om.append("silně zásaditá půda (pH > 8)")
    if ec is not None and ec > 2: om.append("zasolení (EC > 2 dS/m)")
    if esp is not None and esp > 15: om.append("sodická půda (ESP > 15 %)")
    if al is not None and al > 20: om.append("riziko toxicity hliníku")
    if bulk is not None and bulk > 1.45: om.append("utužení (obj. hm. > 1,45)")
    if coarse is not None and coarse > 15: om.append("skeletnatost > 15 %")
    if oc is not None and oc < 1.0: om.append("nízký obsah organického uhlíku")
    if ks and taw < 80: om.append("malá zásoba vody v kořenové zóně")
    if hsg == "D": om.append("vysoké riziko povrchového odtoku")
    return {"root_cm": root, "taw_mm": rd(taw, 0), "soc_0_30": rd(soc30, 1), "soc_0_100": rd(soc100, 1),
            "ksat_min": kmin, "hsg": hsg, "cn2": {"A": 67, "B": 78, "C": 85, "D": 89}.get(hsg), "limity": om}

def nacti_root_depth():
    """Kód -> střední hloubka (cm) z textu D_ROOT_DEPTH.csv (čísla v textu; u rozsahu průměr). Zkontrolujte výsledek!"""
    m = {}
    if os.path.exists(CSV_ROOT):
        with open(CSV_ROOT, encoding="utf-8") as f:
            rdr = csv.reader(f); h = [x.upper() for x in next(rdr)]
            ki = h.index("CODE") if "CODE" in h else 0; vi = h.index("VALUE") if "VALUE" in h else 1
            for row in rdr:
                cis = [int(x) for x in re.findall(r"\d+", row[vi])][:2]
                if cis: m[row[ki].strip()] = sum(cis) / len(cis)
    print("Hloubky kořenů z číselníku:", m)
    return m

# ---------- hlavní program ----------
root_map = nacti_root_depth()
print("1. KROK: SMU...")
smu_data = {}
with open(CSV_SMU, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        sid = int(row["HWSD2_SMU_ID"]); smu_data.setdefault(sid, {"smu_id": sid, "komponenty": []})

print("2. KROK: vrstvy + odvozené parametry...")
NOVA = {"cec_soil": "CEC_SOIL", "cec_eff": "CEC_EFF", "teb": "TEB", "bsat": "BSAT", "alum_sat": "ALUM_SAT",
        "esp": "ESP", "caco3": "TCARBON_EQ", "gypsum": "GYPSUM", "ref_bulk": "REF_BULK",
        "tex_usda": "TEXTURE_USDA", "tex_soter": "TEXTURE_SOTER"}
STARE = {"top": "TOPDEP", "bot": "BOTDEP", "clay": "CLAY", "sand": "SAND", "silt": "SILT", "coarse": "COARSE",
         "bulk": "BULK", "oc": "ORG_CARBON", "ph": "PH_WATER", "total_n": "TOTAL_N", "cn_ratio": "CN_RATIO",
         "cec_clay": "CEC_CLAY", "elec_cond": "ELEC_COND"}
komp_mapa = {}
with open(CSV_LAYERS, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        sid = int(row["HWSD2_SMU_ID"])
        if sid not in smu_data: continue
        seq = row.get("SEQUENCE", "1").strip(); lid = row.get("LAYER", "").strip(); key = f"{sid}_{seq}"
        if key not in komp_mapa:
            komp_mapa[key] = {"sequence": seq, "share": row.get("SHARE", "100"), "root_depth": row.get("ROOT_DEPTH", ""),
                "awc_total": row.get("AWC", ""), "drainage": row.get("DRAINAGE", ""), "fao90": row.get("FAO90", ""),
                "wrb2": row.get("WRB2", ""), "wrb4": row.get("WRB4", ""), "wrb_phases": row.get("WRB_PHASES", ""),
                "il": row.get("IL", ""), "swr": row.get("SWR", ""), "add_prop": row.get("ADD_PROP", ""),
                "phase1": row.get("PHASE1", ""), "phase2": row.get("PHASE2", ""), "roots": row.get("ROOTS", ""),
                "layers": {}}
        if lid:
            l = {k: row.get(v, "") for k, v in {**STARE, **NOVA}.items()}
            l.update(odvod_vrstvu(l))
            komp_mapa[key]["layers"][lid] = l

for key, k in komp_mapa.items():
    k["odv"] = odvod_komponentu(k, root_map)
    smu_data[int(key.split("_")[0])]["komponenty"].append(k)

print("3. KROK: souhrn SMU (vážený podílem), komprese, index...")
index_pud, kontejner, offset = {}, bytearray(), 0
for sid, p in smu_data.items():
    p["komponenty"].sort(key=lambda x: int(x["sequence"]))
    w = [(num(k["share"]) or 0, k["odv"]) for k in p["komponenty"] if k["odv"]["taw_mm"]]
    sw = sum(x for x, _ in w)
    if sw > 0:
        p["souhrn"] = {n: rd(sum(x * o[n] for x, o in w if o[n] is not None) / sw, 1) for n in ("taw_mm", "soc_0_30", "soc_0_100")}
    blok = gzip.compress(json.dumps(p, ensure_ascii=False, separators=(",", ":")).encode("utf-8"), compresslevel=6)
    index_pud[sid] = {"offset": offset, "delka": len(blok)}
    kontejner.extend(blok); offset += len(blok)

with open(VYSTUP_DATA, "wb") as f: f.write(kontejner)
with open(VYSTUP_INDEX, "w", encoding="utf-8") as f: json.dump(index_pud, f, ensure_ascii=False)
print(f"HOTOVO: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB), profilů: {len(index_pud)}")
