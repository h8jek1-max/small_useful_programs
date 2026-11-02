import csv
import json
import glob
import os

SLOZKA_CSV = "./" # Složka, kde máte rozbalené CSV tabulky FAO
VYSTUP_JSON = "ciselniky.json"

# Seznam číselníkových tabulek, které chceme sloučit
ciselnikove_tabulky = [
    "D_ADD_PROP", "D_AWC", "D_COVERAGE", "D_DRAINAGE", "D_FAO90", 
    "D_IL", "D_KOPPEN", "D_PHASE", "D_ROOTS", "D_ROOT_DEPTH", 
    "D_SWR", "D_TEXTURE", "D_TEXTURE_SOTER", "D_TEXTURE_USDA", 
    "D_WRB2", "D_WRB4", "D_WRB_PHASES"
]

vsechny_ciselniky = {}

for nazev in ciselnikove_tabulky:
    cesta = os.path.join(SLOZKA_CSV, f"{nazev}.csv")
    if not os.path.exists(cesta):
        print(f"Varování: Číselník {nazev}.csv nebyl ve složce nalezen, přeskakuji.")
        continue
        
    print(f"Zpracovávám číselník: {nazev}.csv")
    vsechny_ciselniky[nazev] = {}
    
    with open(cesta, mode="r", encoding="utf-8") as f:
        # Detekujeme, jaké sloupce tabulka má (některé mají CODE, jiné ID nebo SYMBOL)
        reader = csv.reader(f)
        hlavicka = [h.upper() for h in next(reader)]
        
        # Určíme, který sloupec je klíč a který hodnota (VALUE)
        klic_idx = 0
        if "CODE" in hlavicka: klic_idx = hlavicka.index("CODE")
        elif "ID" in hlavicka: klic_idx = hlavicka.index("ID")
        elif "SYMBOL" in hlavicka: klic_idx = hlavicka.index("SYMBOL")
            
        val_idx = hlavicka.index("VALUE") if "VALUE" in hlavicka else 1

        for radek in reader:
            if len(radek) > max(klic_idx, val_idx):
                klic = radek[klic_idx].strip()
                hodnota = radek[val_idx].strip()
                vsechny_ciselniky[nazev][klic] = hodnota

with open(VYSTUP_JSON, "w", encoding="utf-8") as f:
    json.dump(vsechny_ciselniky, f, ensure_ascii=False, indent=2)

print(f"\nHOTOVO! Všechny číselníky byly sloučeny do souboru: {VYSTUP_JSON}")

