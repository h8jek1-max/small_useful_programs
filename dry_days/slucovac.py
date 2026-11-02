import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

SLOZKA_VSTUP = "databaze/"
SOUBOR_DATA = "celysvet_cdd.gz.bin"
SOUBOR_REJSTRIK = "registrace_cdd.json"

TILE_SIZE = 45
NCOLS = 720
NROWS = 360

# Hledáme soubory CDD
soubory = sorted(glob.glob(os.path.join(SLOZKA_VSTUP, "CDD_*.nc")))
if not soubory:
    print(f"Chyba: Žádné CDD soubory .nc ve složce '{SLOZKA_VSTUP}'!")
    exit()

print(f"Nalezeno {len(soubory)} souborů pro sezónní index CDD. Zahajuji zpracování celého severu...")

vsechny_kroky = []
for cesta_souboru in soubory:
    print(f"  -> Načítám: {os.path.basename(cesta_souboru)}")
    with xr.open_dataset(cesta_souboru) as ds:
        data_array = ds["CDD"].squeeze()
        if data_array.dims != ('time', 'lat', 'lon'):
            data_array = data_array.transpose('time', 'lat', 'lon')
        data = np.array(data_array.values)
        
    # Vyčištění hodnot (65535 nebo -999 pro oceány převedeme na jednotnou masku 255) 
    data_vycistena = np.where(np.isnan(data) | (data < 0), -999.0, data)
    kvantovana_data = np.where(data_vycistena == -999.0, 255, np.clip(data_vycistena, 0, 254)).astype(np.uint8)
    
    vsechny_kroky.append(kvantovana_data)

# Spojíme časovou osu (tvar: Sezónní kroky, Lat, Lon)
matice_cdd = np.concatenate(vsechny_kroky, axis=0)
celkovy_pocet_kroku = int(matice_cdd.shape[0])
print(f"\nCelkový počet zindexovaných sezónních kroků: {celkovy_pocet_kroku}")

# KLÍČOVÁ KOREKCE SEVERU: Překlopíme osu Y zrcadlově, aby index 0 byl severní pól (+89.75)
matice_cdd = np.flip(matice_cdd, axis=1)

# Překlopíme na tvar (Lat, Lon, Čas), aby časová osa jednoho pixelu běžela v kuse
matice_cdd = np.transpose(matice_cdd, (1, 2, 0))

rejstrik = []
binarni_data = bytearray()
aktualni_offset = 0

pocet_dlazdic_y = NROWS // TILE_SIZE
pocet_dlazdic_x = NCOLS // TILE_SIZE

for ty in range(pocet_dlazdic_y):
    for tx in range(pocet_dlazdic_x):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        
        dlaždice = matice_cdd[y_od:y_do, x_od:x_do, :]
        
        if np.all(dlaždice == 255):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        surove_bajty = dlaždice.tobytes()
        komprimovane = gzip.compress(surove_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
        binarni_data.extend(komprimovane)
        aktualni_offset += velikost_bloku

with open(SOUBOR_DATA, "wb") as f: f.write(binarni_data)
with open(SOUBOR_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({"pocet_kroku": celkovy_pocet_kroku, "dlaždice": rejstrik}, f, ensure_ascii=False)

print("==================================================")
print("SEZÓNNÍ DATABÁZE CDD ÚSPĚŠNĚ VYTVOŘENA!")
print(f"Soubor: {SOUBOR_DATA} ({os.path.getsize(SOUBOR_DATA)/(1024*1024):.2f} MB)")
print("==================================================")
