import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

VYSTUP_DATA = "celysvet_srazky.gz.bin"
VYSTUP_REJSTRIK = "registrace_srazky.json"

TILE_SIZE = 32
NCOLS = 256
NROWS = 128
VELIKOST_DAVKY_MESICU = 100 # Zpracování po malých časových úsecích pro úsporu RAM

# Automatické vyhledání souboru srážek
nalezeno = glob.glob("pr_*.nc")
if not nalezeno and os.path.exists("databaze"):
    nalezeno = glob.glob("databaze/pr_*.nc")

if not nalezeno:
    print("Chyba: Nebyl nalezen žádný soubor s maskou 'pr_*.nc'!")
    exit()

SOUBOR_NC = list(nalezeno)[0]
print(f"Zahajuji ŠETRNOU DAVKOVOU konverzi srážek: {os.path.basename(SOUBOR_NC)}")

# Inicializace prázdné výsledné globální matice v uint8 (zabírá minimum RAM)
with xr.open_dataset(SOUBOR_NC) as ds:
    celkovy_pocet_mesicu = int(ds.dims['time'])
    # Připravíme si finální pole, do kterého budeme dávkově sypat výsledky
    matice_srazky = np.zeros((NROWS, NCOLS, celkovy_pocet_mesicu), dtype=np.uint8)

sekund_v_mesici = 30.4375 * 86400.0

# --- DÁVKOVÉ ZPRACOVÁNÍ PO ČASOVÉ OSE (Úspora RAM a CPU) ---
for start_m in range(0, celkovy_pocet_mesicu, VELIKOST_DAVKY_MESICU):
    konec_m = min(celkovy_pocet_mesicu, start_m + VELIKOST_DAVKY_MESICU)
    
    # Načteme z NetCDF pouze vyříznutý časový úsek (Slicing)
    with xr.open_dataset(SOUBOR_NC) as ds:
        data_array = ds["pr"].isel(time=slice(start_m, konec_m)).squeeze()
        if data_array.dims != ('time', 'lat', 'lon'):
            data_array = data_array.transpose('time', 'lat', 'lon')
        surove_kg_davka = np.array(data_array.values)

    # Přepočet na mm/měsíc pro danou dávku
    mm_davka = surove_kg_davka * sekund_v_mesici
    
    # Ošetření neplatných hodnot
    maska_nan = np.isnan(mm_davka) | (mm_davka < 0.0)
    mm_cista = np.where(maska_nan, 0.0, mm_davka)

    # BLESKOVÁ KVANTOVÁNÍ POMOCÍ ČISTÉHO NUMPY (Rychlost v řádu milisekund)
    # Výchozí stav: zaokrouhlené hodnoty 1:1 (0 až 240 mm)
    bajtova_davka = np.round(mm_cista)
    
    # Postupné aplikování progresivních masek pro vysoké tropické srážky nad 240 mm
    bajtova_davka = np.where((mm_cista > 240.0) & (mm_cista <= 255.0), 241, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 255.0) & (mm_cista <= 275.0), 242, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 275.0) & (mm_cista <= 300.0), 243, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 300.0) & (mm_cista <= 330.0), 244, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 330.0) & (mm_cista <= 365.0), 245, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 365.0) & (mm_cista <= 405.0), 246, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 405.0) & (mm_cista <= 450.0), 247, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 450.0) & (mm_cista <= 500.0), 248, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 500.0) & (mm_cista <= 550.0), 249, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 550.0) & (mm_cista <= 600.0), 250, bajtova_davka)
    bajtova_davka = np.where(mm_cista > 600.0, 254, bajtova_davka)
    
    # Dosazení systémové masky pro chyby/NaN
    bajtova_davka = np.where(maska_nan, 255, bajtova_davka).astype(np.uint8)

    # Překlopení severního pólu pro Leaflet kompatibilitu
    bajtova_davka = np.flip(bajtova_davka, axis=1)

    # Transpozice dávky na rozměr (Lat, Lon, Čas)
    bajtova_davka_transponovana = np.transpose(bajtova_davka, (1, 2, 0))

    # Uložení rozpracované dávky do finální velké matice
    matice_srazky[:, :, start_m:konec_m] = bajtova_davka_transponovana

    # --- UKAZATEL PRŮBĚHU V % ---
    procento = (konec_m / celkovy_pocet_mesicu) * 100
    print(f"  -> Zpracováno: {konec_m}/{celkovy_pocet_mesicu} měsíců ({procento:.1f} %)")

print("Dávkové zpracování dokončeno. Balím výslednou matici do dlaždic...")

# --- PORCOVÁNÍ DO GZIP DLAŽDIC (Stejný systém jako u výparu) ---
rejstrik = []
binarni_kontainer = bytearray()
aktualni_offset = 0

pocet_dlazdic_y = NROWS // TILE_SIZE
pocet_dlazdic_x = NCOLS // TILE_SIZE

for ty in range(pocet_dlazdic_y):
    for tx in range(pocet_dlazdic_x):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        
        dlazdice = matice_srazky[y_od:y_do, x_od:x_do, :]
        
        if np.all(dlazdice == 255):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        surove_bajty = dlazdice.tobytes()
        komprimovane = gzip.compress(surove_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
        binarni_kontainer.extend(komprimovane)
        aktualni_offset += velikost_bloku

with open(VYSTUP_DATA, "wb") as f: 
    f.write(binarni_kontainer)
    
with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({"pocet_mesicu": celkovy_pocet_mesicu, "dlaždice": rejstrik}, f, ensure_ascii=False)

print("==================================================")
print("ULTRA-RYCHLÁ DATABÁZE SRÁŽEK ÚSPĚŠNĚ SESTAVENA!")
print(f"Soubor: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB)")
print("==================================================")
