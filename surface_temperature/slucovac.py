import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

# Konfigurace - skript prohledá složku a najde soubor automaticky
VYSTUP_DATA = "celysvet_teplota.gz.bin"
VYSTUP_REJSTRIK = "registrace_teploty.json"

TILE_SIZE = 32
NCOLS = 256
NROWS = 128

# Automatické nalezení souboru teploty ve složce
nalezené_soubory = glob.glob("tas_*.nc")
if not nalezené_soubory:
    # Zkusíme se podívat i do podsložky databaze, pokud existuje
    if os.path.exists("databaze"):
        nalezené_soubory = glob.glob("databaze/tas_*.nc")

if not nalezené_soubory:
    print("Chyba: V adresáři nebyl nalezen žádný soubor teploty s maskou 'tas_*.nc'!")
    print("Zkontrolujte, zda jste soubor správně rozbalili ze ZIP archivu a zda má příponu .nc")
    exit()

SOUBOR_NC = nalezené_soubory[0]
print(f"Zahajuji konverzi globální teplotní databáze ze souboru: {SOUBOR_NC}")

with xr.open_dataset(SOUBOR_NC) as ds:
    # Vynutíme si standardní pořadí dimenzí (time, lat, lon)
    data_array = ds["tas"].squeeze()
    if data_array.dims != ('time', 'lat', 'lon'):
        data_array = data_array.transpose('time', 'lat', 'lon')
    
    # Načtení surových Kelvinů do numpy matice
    surove_kelviny = np.array(data_array.values)

celkovy_pocet_mesicu = int(surove_kelviny.shape[0])
print(f"Detekováno: {celkovy_pocet_mesicu} měsíců (2015-2099). Přepočítávám na Celsia...")

# Převod Kelvin -> Celsius
celsia = surove_kelviny - 273.15

# ÚPRAVA PRO PŮLSTUPŇOVOU PŘESNOST (0.5 °C):
celsia_vycistena = np.where(np.isnan(celsia) | (celsia < -50.0) | (celsia > 50.0), -999.0, celsia)
# Vynásobíme dvěma a posuneme o 100 nahoru, čímž eliminujeme záporná čísla do 1 bajtu
bajtova_data = np.where(celsia_vycistena == -999.0, 255, np.round(celsia_vycistena * 2.0 + 100.0)).astype(np.uint8)


# ORIENTAČNÍ KOREKCE SEVERU: Překlopíme osu Y (lat), aby severní pól začínal na indexu 0
bajtova_data = np.flip(bajtova_data, axis=1)

# Překlopíme dimenze na (Lat, Lon, Čas) pro sekvenční fseek čtení historie pixelu
matice_teploty = np.transpose(bajtova_data, (1, 2, 0))

print(f"Matice připravena v rozměru: {matice_teploty.shape}. Balím do dlaždic {TILE_SIZE}x{TILE_SIZE}...")

rejstrik = []
binarni_kontainer = bytearray()
aktualni_offset = 0

pocet_dlazdic_y = NROWS // TILE_SIZE # 128 / 32 = 4
pocet_dlazdic_x = NCOLS // TILE_SIZE # 256 / 32 = 8

for ty in range(pocet_dlazdic_y):
    for tx in range(pocet_dlazdic_x):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        
        dlazdice = matice_teploty[y_od:y_do, x_od:x_do, :]
        
        # Pokud je dlaždice kompletně tvořena oceánem (maska 255)
        if np.all(dlazdice == 255):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        surove_bajty = dlazdice.tobytes()
        komprimovane = gzip.compress(surove_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
        binarni_kontainer.extend(komprimovane)
        aktualni_offset += velikost_bloku

# Zápis souborů na disk
with open(VYSTUP_DATA, "wb") as f: 
    f.write(binarni_kontainer)
    
with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({"pocet_mesicu": celkovy_pocet_mesicu, "dlaždice": rejstrik}, f, ensure_ascii=False)

print("==================================================")
print("TEPLOTNÍ DATABÁZE PROJEKCE ÚSPĚŠNĚ SESTAVENA!")
print(f"Soubor: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB)")
print("==================================================")
