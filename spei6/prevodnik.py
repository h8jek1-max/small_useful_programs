import glob
import os
import json
import gzip
import shutil
import numpy as np
import xarray as xr

# Nastavení cest
SLOZKA_VSTUP = "databaze/"
SOUBOR_DATA = "celysvet_sucho_gzip.bin"
SOUBOR_REJSTRIK = "registrace_gzip.json"
SLOZKA_TEMP = "docasne_dlaždice/"

TILE_SIZE = 45 
pocet_dlaždic_y = 720 // TILE_SIZE  # 16
pocet_dlaždic_x = 1440 // TILE_SIZE # 32

if os.path.exists(SLOZKA_TEMP):
    shutil.rmtree(SLOZKA_TEMP)
os.makedirs(SLOZKA_TEMP, exist_ok=True)

soubory = sorted(glob.glob(os.path.join(SLOZKA_VSTUP, "*.nc")))
if not soubory:
    print(f"Chyba: Žádné soubory ve složce '{SLOZKA_VSTUP}'!")
    exit()

print(f"Nalezeno {len(soubory)} souborů sucha. Spouštím DÁVKOVÝ sken KOMPLETNÍHO SVĚTA...")

celkovy_pocet_mesicu = 0

for index_s, cesta_souboru in enumerate(soubory):
    print(f" Zpracovávám [{index_s+1}/{len(soubory)}]: {os.path.basename(cesta_souboru)}")
    
    with xr.open_dataset(cesta_souboru) as ds:
        data = np.array(ds["spe06"].squeeze().values)
        
    if data.ndim == 2:
        data = np.expand_dims(data, axis=0)
        
    mesicu_v_souboru = data.shape[0]
    celkovy_pocet_mesicu += mesicu_v_souboru
    
    # 4-bitová kvantizace (0 až 15)
    data_vycistena = np.where(np.isnan(data) | (data < -4.5) | (data > 4.5), -999.0, data)
    data_oříznutá = np.clip(data_vycistena, -4.0, 4.0)
    kvantovana_roky = np.where(
        data_vycistena == -999.0,
        15,
        ((data_oříznutá + 4.0) * (14.0 / 8.0)).astype(np.uint8)
    ).astype(np.uint8)
    
    # !!! ZMĚNA: ÚPLNĚ JSME SMAZALI SEVERNÍ OŘEZ PLOCHY !!!
    # Ponechali jsme pouze logický ořez Antarktidy (pod -60° j.š.), kde pole ani lesy nerostou
    for m in range(mesicu_v_souboru):
        for lat_idx in range(720):
            lat_val = 89.875 - (lat_idx * 0.25)
            if lat_val < -60.0:
                kvantovana_roky[m, lat_idx, :] = 15

    # Rozmístění do dočasných složek
    for ty in range(pocet_dlaždic_y):
        for tx in range(pocet_dlaždic_x):
            y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
            x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
            
            dlaždice_rok = kvantovana_roky[:, y_od:y_do, x_od:x_do]
            dlaždice_rok = np.transpose(dlaždice_rok, (1, 2, 0))
            
            jmeno_souboru = f"tile_{ty}_{tx}.raw"
            cesta_temp_souboru = os.path.join(SLOZKA_TEMP, jmeno_souboru)
            
            with open(cesta_temp_souboru, "ab") as tf:
                tf.write(dlaždice_rok.tobytes())

print(f"\nCelkový počet měsíců: {celkovy_pocet_mesicu}. Sestavuji GZIP kontainer...")

def zabal_4bit(surove_bajty):
    ploche_pole = np.frombuffer(surove_bajty, dtype=np.uint8)
    if len(ploche_pole) % 2 != 0:
        ploche_pole = np.append(ploche_pole, 15)
    horni = ploche_pole[0::2] << 4
    dolni = ploche_pole[1::2]
    return (horni | dolni).astype(np.uint8).tobytes()

rejstrik = []
binarni_data_final = bytearray()
aktualni_offset = 0

for ty in range(pocet_dlaždic_y):
    for tx in range(pocet_dlaždic_x):
        jmeno_souboru = f"tile_{ty}_{tx}.raw"
        cesta_temp_souboru = os.path.join(SLOZKA_TEMP, jmeno_souboru)
        
        with open(cesta_temp_souboru, "rb") as tf:
            surove_bajty_dlaždice = tf.read()
            
        pole_k_overeni = np.frombuffer(surove_bajty_dlaždice, dtype=np.uint8)
        if np.all(pole_k_overeni == 15):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        kompaktni_bajty = zabal_4bit(surove_bajty_dlaždice)
        komprimovane = gzip.compress(kompaktni_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
        binarni_data_final.extend(komprimovane)
        aktualni_offset += velikost_bloku

with open(SOUBOR_DATA, "wb") as f: f.write(binarni_data_final)
with open(SOUBOR_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({"pocet_mesicu": celkovy_pocet_mesicu, "dlaždice": rejstrik}, f, ensure_ascii=False)

shutil.rmtree(SLOZKA_TEMP)
print(f"HOTOVO! Databáze sucha uložena včetně celého severního pólu: {os.path.getsize(SOUBOR_DATA)/(1024*1024):.2f} MB")
