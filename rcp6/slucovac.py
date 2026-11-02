import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

SLOZKA_VSTUP = "./"
SOUBOR_DATA = "celysvet_gsl.gz.bin"
SOUBOR_REJSTRIK = "registrace_gsl.json"

TILE_SIZE = 45
NCOLS = 720
NROWS = 360 # Celý svět od pólu k pólu

soubory = sorted(glob.glob(os.path.join(SLOZKA_VSTUP, "GSL_*.nc")))
if not soubory:
    print(f"Chyba: Žádné GSL soubory v '{SLOZKA_VSTUP}'!")
    exit()

vsechny_roky = []
for cesta_souboru in soubory:
    print(f"  -> Načítám: {os.path.basename(cesta_souboru)}")
    with xr.open_dataset(cesta_souboru) as ds:
        data_array = ds["GSL"].squeeze()
        if data_array.dims != ('time', 'lat', 'lon'):
            data_array = data_array.transpose('time', 'lat', 'lon')
        data = np.array(data_array.values)
    
    # Bezpečné ošetření NaN, hodnoty podlahy a stropu
    data_vycistena = np.where(np.isnan(data) | (data < 0), 65535, data)
    # Ukládáme jako Unsigned 16-bit int (0 až 65535), takže 365 dní projde bez jakéhokoli ořezu!
    vsechny_roky.append(data_vycistena.astype(np.uint16))

matice_gsl = np.concatenate(vsechny_roky, axis=0)
celkovy_pocet_let = int(matice_gsl.shape[0])

# OPRAVA SMĚRU: Překlopíme osu Y zrcadlově, aby index 0 byl severní pól (+89.75)
matice_gsl = np.flip(matice_gsl, axis=1)
matice_gsl = np.transpose(matice_gsl, (1, 2, 0))

print(f"Tvar matice s celým severem: {matice_gsl.shape}")

rejstrik = []
binarni_data = bytearray()
aktualni_offset = 0

pocet_dlazdic_y = NROWS // TILE_SIZE
pocet_dlazdic_x = NCOLS // TILE_SIZE

for ty in range(pocet_dlazdic_y):
    for tx in range(pocet_dlazdic_x):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        
        dlaždice = matice_gsl[y_od:y_do, x_od:x_do, :]
        
        if np.all(dlaždice == 65535):
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
    json.dump({"pocet_let": celkovy_pocet_let, "dlaždice": rejstrik}, f, ensure_ascii=False)

print(f"HOTOVO! GSL databáze opravena a uložena pro kompletní severní polokouli.")
