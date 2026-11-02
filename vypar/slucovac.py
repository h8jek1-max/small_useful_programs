import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

# Konfigurace
VYSTUP_DATA = "celysvet_evapo.gz.bin"
VYSTUP_REJSTRIK = "registrace_evapo.json"

TILE_SIZE = 32
NCOLS = 256
NROWS = 128

# Automatické vyhledání souboru evapotranspirace ve složce
nalezeno = glob.glob("evspsbl_*.nc")
if not nalezeno and os.path.exists("databaze"):
    nalezeno = glob.glob("databaze/evspsbl_*.nc")

if not nalezeno:
    print("Chyba: Nebyl nalezen žádný soubor s maskou 'evspsbl_*.nc'!")
    exit()

SOUBOR_NC = nalezeno[0]
print(f"Zahajuji konverzi globální evapotranspirace ze souboru: {SOUBOR_NC}")

with xr.open_dataset(SOUBOR_NC) as ds:
    data_array = ds["evspsbl"].squeeze()
    if data_array.dims != ('time', 'lat', 'lon'):
        data_array = data_array.transpose('time', 'lat', 'lon')
    
    # Načtení surových hodnot (kg/m2/s)
    surove_kg = np.array(data_array.values)

celkovy_pocet_mesicu = int(surove_kg.shape[0])
print(f"Detekováno {celkovy_pocet_mesicu} měsíců. Přepočítávám na mm/měsíc...")

# Přepočtová konstanta: kg/m2/s -> mm/měsíc (průměrný měsíc má 30.4375 dní)
# 30.4375 dní * 86400 sekund = 2630016 sekund
sekund_v_mesici = 30.4375 * 86400.0
mm_za_mesic = surove_kg * sekund_v_mesici

# Ošetření neplatných hodnot a masek oceánu
mm_vycistena = np.where(np.isnan(mm_za_mesic) | (mm_za_mesic < 0.0) | (mm_za_mesic > 500.0), -999.0, mm_za_mesic)
bajtova_data = np.where(mm_vycistena == -999.0, 255, np.round(mm_vycistena)).astype(np.uint8)

# ORIENTAČNÍ KOREKCE SEVERU: Překlopíme osu Y (lat), aby index 0 začínal na severním pólu
bajtova_data = np.flip(bajtova_data, axis=1)

# Překlopíme dimenze na (Lat, Lon, Čas) pro plynulé lineární vyřezávání historie pixelu
matice_evapo = np.transpose(bajtova_data, (1, 2, 0))

print(f"Matice připravena v rozměru: {matice_evapo.shape}. Balím do dlaždic {TILE_SIZE}x{TILE_SIZE}...")

rejstrik = []
binarni_kontainer = bytearray()
aktualni_offset = 0

pocet_dlazdic_y = NROWS // TILE_SIZE # 4
pocet_dlazdic_x = NCOLS // TILE_SIZE # 8

for ty in range(pocet_dlazdic_y):
    for tx in range(pocet_dlazdic_x):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        
        dlazdice = matice_evapo[y_od:y_do, x_od:x_do, :]
        
        if np.all(dlazdice == 255):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        surove_bajty = dlazdice.tobytes()
        komprimovane = gzip.compress(surove_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
        binarni_kontainer.extend(komprimovane)
        aktualni_offset += velikost_bloku

with open(VYSTUP_DATA, "wb") as f: f.write(binarni_kontainer)
with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({"pocet_mesicu": celkovy_pocet_mesicu, "dlaždice": rejstrik}, f, ensure_ascii=False)

print("==================================================")
print("DATABÁZE EVAPOTRANSPIRACE ÚSPĚŠNĚ SESTAVENA!")
print(f"Soubor: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB)")
print("==================================================")
