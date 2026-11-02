
import os
import csv
import json
import gzip

# Cesty k souborům
CSV_SMU = "HWSD2_SMU.csv"
CSV_LAYERS = "HWSD2_LAYERS.csv"
VYSTUP_DATA = "pudy_data.gz.bin"
VYSTUP_INDEX = "pudy_index.json"

print("1. KROK: Načítání hlavních metadat ze SMU.csv (řádek po řádku)...")
smu_data = {}
with open(CSV_SMU, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        smu_id = int(row["HWSD2_SMU_ID"])
        smu_data[smu_id] = {
            "root_depth": row.get("ROOT_DEPTH", ""),
            "awc": row.get("AWC", ""),
            "drainage": row.get("DRAINAGE", ""),
            "layers": {} # Zde nasbíráme vrstvy D1-D7
        }

print("2. KROK: Spojování s detaily vrstev z LAYERS.csv...")
with open(CSV_LAYERS, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        smu_id = int(row["HWSD2_SMU_ID"])
        if smu_id not in smu_data:
            continue
            
        layer_id = row.get("LAYER", "").strip()
        if not layer_id:
            continue
            
        # Extrahujeme nejdůležitější fyzikální a chemické parametry pro precizní agro-simulace
        smu_data[smu_id]["layers"][layer_id] = {
            "top": row.get("TOPDEP", ""),
            "bot": row.get("BOTDEP", ""),
            "clay": row.get("CLAY", ""),
            "sand": row.get("SAND", ""),
            "silt": row.get("SILT", ""),
            "oc": row.get("ORG_CARBON", ""),
            "ph": row.get("PH_WATER", ""),
            "awc_layer": row.get("AWC", "")
        }

print("3. KROK: Balení profilů do binárního GZIP souboru a tvorba indexu...")
index_pud = {}
binarni_kontainer = bytearray()
aktualni_offset = 0

for smu_id, profil in smu_data.items():
    # Převedeme celý komplexní profil (metadata + všech 7 vrstev) na úsporný JSON text
    json_profil = json.dumps(profil, ensure_ascii=False)
    # Zkomprimujeme text profilu pomocí Gzipu
    komprimovany_profil = gzip.compress(json_profil.encode("utf-8"), compresslevel=6)
    delka_bloku = len(komprimovany_profil)
    
    # Uložíme si pozici do indexu
    index_pud[smu_id] = {
        "offset": aktualni_offset,
        "delka": delka_bloku
    }
    
    binarni_kontainer.extend(komprimovany_profil)
    aktualni_offset += delka_bloku

# Zápis výsledků na disk
with open(VYSTUP_DATA, "wb") as f:
    f.write(binarni_kontainer)

with open(VYSTUP_INDEX, "w", encoding="utf-8") as f:
    json.dump(index_pud, f, ensure_ascii=False)

velikost_mb = os.path.getsize(VYSTUP_DATA) / (1024 * 1024)
print("==================================================")
print("ÚSPĚŠNĚ DOKONČENO!")
print(f"Vytvořen datový soubor: {VYSTUP_DATA} ({velikost_mb:.2f} MB)")
print(f"Vytvořen index: {VYSTUP_INDEX}")
print("==================================================")
