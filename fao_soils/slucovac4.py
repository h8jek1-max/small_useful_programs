
import csv
import json
import gzip
import os

CSV_SMU = "HWSD2_SMU.csv"
CSV_LAYERS = "HWSD2_LAYERS.csv"
VYSTUP_DATA = "pudy_data_komplet.gz.bin"
VYSTUP_INDEX = "pudy_index_komplet.json"

print("1. KROK: Načítání hlavních metadat ze SMU.csv (řádek po řádku)...")
smu_data = {}
with open(CSV_SMU, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        smu_id = int(row["HWSD2_SMU_ID"])
        if smu_id not in smu_data:
            smu_data[smu_id] = {
                "smu_id": smu_id,
                "komponenty": []
            }

print("2. KROK: Spojování s detaily vrstev a novými agro-parametry z LAYERS.csv...")
komponenty_mapa = {}

with open(CSV_LAYERS, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        smu_id = int(row["HWSD2_SMU_ID"])
        if smu_id not in smu_data:
            continue
            
        seq = row.get("SEQUENCE", "1").strip()
        layer_id = row.get("LAYER", "").strip()
        comp_klic = f"{smu_id}_{seq}"
        
        if comp_klic not in komponenty_mapa:
            komponenty_mapa[comp_klic] = {
                "sequence": seq,
                "share": row.get("SHARE", "100"),
                "root_depth": row.get("ROOT_DEPTH", ""),
                "awc_total": row.get("AWC", ""),
                "drainage": row.get("DRAINAGE", ""),
                "fao90": row.get("FAO90", ""),
                "wrb2": row.get("WRB2", ""),
                "wrb4": row.get("WRB4", ""),
                "wrb_phases": row.get("WRB_PHASES", ""),
                "il": row.get("IL", ""),
                "swr": row.get("SWR", ""),
                "add_prop": row.get("ADD_PROP", ""),
                "layers": {}
            }
            
        if layer_id:
            # KLÍČOVÁ OPRAVA: Přidány parametry CEC_CLAY, TOTAL_N a ELEC_COND
            komponenty_mapa[comp_klic]["layers"][layer_id] = {
                "top": row.get("TOPDEP", ""),
                "bot": row.get("BOTDEP", ""),
                "clay": row.get("CLAY", ""),
                "sand": row.get("SAND", ""),
                "silt": row.get("SILT", ""),
                "coarse": row.get("COARSE", ""),
                "bulk": row.get("BULK", ""),
                "oc": row.get("ORG_CARBON", ""),
                "ph": row.get("PH_WATER", ""),
                "total_n": row.get("TOTAL_N", ""),
                "cn_ratio": row.get("CN_RATIO", ""),
                "cec_clay": row.get("CEC_CLAY", ""),
                "elec_cond": row.get("ELEC_COND", "")
            }

for komp_klic, komp_obsah in komponenty_mapa.items():
    smu_id = int(komp_klic.split("_")[0])
    smu_data[smu_id]["komponenty"].append(komp_obsah)

print("3. KROK: Komprimace kompletních profilů a zápis indexu...")
index_pud = {}
binarni_kontainer = bytearray()
aktualni_offset = 0

for smu_id, profil in smu_data.items():
    profil["komponenty"] = sorted(profil["komponenty"], key=lambda x: int(x["sequence"]))
    
    json_profil = json.dumps(profil, ensure_ascii=False)
    komprimovany_profil = gzip.compress(json_profil.encode("utf-8"), compresslevel=6)
    delka_bloku = len(komprimovany_profil)
    
    index_pud[smu_id] = {
        "offset": aktualni_offset,
        "delka": delka_bloku
    }
    
    binarni_containers = binarni_kontainer.extend(komprimovany_profil)
    aktualni_offset += delka_bloku

with open(VYSTUP_DATA, "wb") as f: f.write(binarni_kontainer)
with open(VYSTUP_INDEX, "w", encoding="utf-8") as f: json.dump(index_pud, f, ensure_ascii=False)

print("==================================================")
print("KOMPLETNÍ INDEXACE ÚSPĚŠNĚ DOKONČENA!")
print(f"Vytvořen soubor: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB)")
print("==================================================")
