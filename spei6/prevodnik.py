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

# Vyčištění a příprava dočasné složky
if os.path.exists(SLOZKA_TEMP):
    shutil.rmtree(SLOZKA_TEMP)
os.makedirs(SLOZKA_TEMP, exist_ok=True)

# Hledání .nc souborů
soubory = sorted(glob.glob(os.path.join(SLOZKA_VSTUP, "*.nc")))
if not soubory:
    print(f"Chyba: Žádné soubory ve složce '{SLOZKA_VSTUP}'!")
    exit()

print(f"Nalezeno {len(soubory)} souborů. Zahajuji DÁVKOVÉ zpracování na disk...")

celkovy_pocet_mesicu = 0

# 1. ETAPA: Zpracování souborů po jednom a plnění dočasných souborů na disku
for index_s, cesta_souboru in enumerate(soubory):
    print(f" Zpracovávám [{index_s+1}/{len(soubory)}]: {os.path.basename(cesta_souboru)}")
    
    # Načteme pouze JEDEN soubor
    with xr.open_dataset(cesta_souboru) as ds:
        data = np.array(ds["spe06"].squeeze().values)
        
    if data.ndim == 2:  # Ošetření pokud by měl soubor jen 1 měsíc
        data = np.expand_dims(data, axis=0)
        
    mesicu_v_souboru = data.shape[0]
    celkovy_pocet_mesicu += mesicu_v_souboru
    
    # Kvantizace na 4 bity (0 až 15) přímo pro tento jeden rok
    data_vycistena = np.where(np.isnan(data) | (data < -4.5) | (data > 4.5), -999.0, data)
    data_oříznutá = np.clip(data_vycistena, -4.0, 4.0)
    kvantovana_roky = np.where(
        data_vycistena == -999.0,
        15,
        ((data_oříznutá + 4.0) * (14.0 / 8.0)).astype(np.uint8)
    ).astype(np.uint8)
    
    # Aplikace zeměpisného ořezu (-60° až +60°)
    for m in range(mesicu_v_souboru):
        for lat_idx in range(720):
            lat_val = 89.875 - (lat_idx * 0.25)
            if lat_val > 60.0 or lat_val < -60.0:
                kvantovana_roky[m, lat_idx, :] = 15

    # Rozmístění dat roku do jednotlivých dočasných souborů dlaždic
    for ty in range(pocet_dlaždic_y):
        for tx in range(pocet_dlaždic_x):
            y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
            x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
            
            # Výřez dlaždice pro tento rok: tvar (měsíce, 45, 45)
            dlaždice_rok = kvantovana_roky[:, y_od:y_do, x_od:x_do]
            
            # Překlopíme na (45, 45, měsíce), aby časová osa běžela v kuse
            dlaždice_rok = np.transpose(dlaždice_rok, (1, 2, 0))
            
            # Zapíšeme surové bajty na konec příslušného dočasného souboru
            jmeno_souboru = f"tile_{ty}_{tx}.raw"
            cesta_temp_souboru = os.path.join(SLOZKA_TEMP, jmeno_souboru)
            
            with open(cesta_temp_souboru, "ab") as tf:
                tf.write(dlaždice_rok.tobytes())

print(f"\nCelkový zjištěný počet měsíců: {celkovy_pocet_mesicu}")
print("\n2. ETAPA: Zabalení dlaždic z disku do finálního GZIP archivu...")

rejstrik = []
binarni_data_final = bytearray()
aktualni_offset = 0

# Pomocná funkce pro sbalení 4-bitových hodnot (dva pixely do 1 bajtu)
def zabal_4bit(surove_bajty):
    ploche_pole = np.frombuffer(surove_bajty, dtype=np.uint8)
    if len(ploche_pole) % 2 != 0:
        ploche_pole = np.append(ploche_pole, 15)
    horni = ploche_pole[0::2] << 4
    dolni = ploche_pole[1::2]
    return (horni | dolni).astype(np.uint8).tobytes()

# Procházíme dočasné soubory, balíme je přes GZIP a stavíme rejstřík
for ty in range(pocet_dlaždic_y):
    for tx in range(pocet_dlaždic_x):
        jmeno_souboru = f"tile_{ty}_{tx}.raw"
        cesta_temp_souboru = os.path.join(SLOZKA_TEMP, jmeno_souboru)
        
        # Načteme nashromážděná data dlaždice za všechny roky z disku
        with open(cesta_temp_souboru, "rb") as tf:
            surove_bajty_dlaždice = tf.read()
            
        # Převedeme pole bajtů na numpy abychom ověřili, zda to není samý oceán (hodnota 15)
        pole_k_overeni = np.frombuffer(surove_bajty_dlaždice, dtype=np.uint8)
        
        # Pokud dlaždice obsahuje POUZE oceán/póly, do finálního souboru ji vůbec nedáme
        if np.all(pole_k_overeni == 15):
            rejstrik.append({"offset": aktualni_offset, "delka": 0})
            continue
            
        # Zabalíme dva 4-bitové pixely do jednoho bajtu
        kompaktni_bajty = zabal_4bit(surove_bajty_dlaždice)
        
        # Zkomprimujeme pomocí GZIP
        komprimovane = gzip.compress(kompaktni_bajty, compresslevel=6)
        velikost_bloku = len(komprimovane)
        
        rejstrik.append({
            "offset": aktualni_offset,
            "delka": velikost_bloku
        })
        
        binarni_data_final.extend(komprimovane)
        aktualni_offset += velikost_bloku

# 3. ETAPA: Finální zápis a úklid
print("\n3. ETAPA: Ukládání výsledků a čištění disku...")
with open(SOUBOR_DATA, "wb") as f:
    f.write(binarni_data_final)

meta_struktura = {
    "pocet_mesicu": celkovy_pocet_mesicu,
    "dlaždice": rejstrik
}

with open(SOUBOR_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump(meta_struktura, f, ensure_ascii=False)

# Smazání dočasné složky z disku
shutil.rmtree(SLOZKA_TEMP)

velikost_mb = os.path.getsize(SOUBOR_DATA) / (1024 * 1024)
print(f"==================================================")
print(f"ÚSPĚŠNĚ DOKONČENO BEZ ZAMRZNUTÍ PC!")
print(f"Vytvořen rejstřík: {SOUBOR_REJSTRIK}")
print(f"Finální velikost GZIP databáze na disku: {velikost_mb:.2f} MB")
print(f"==================================================")
