import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

# --- NASTAVENÍ CEST ---
SLOZKA_VSTUP = "databaze/"
SOUBOR_DATA = "celysvet_sucho_landmask.bin"
SOUBOR_REJSTRIK = "registrace_landmask.json"
SOUBOR_MASKY = "land_mask_1440x720.raw" # Statická maska pevniny

TILE_SIZE = 45 
X_POINTS, Y_POINTS = 1440, 720
pocet_dlaždic_y = Y_POINTS // TILE_SIZE  # 16
pocet_dlaždic_x = X_POINTS // TILE_SIZE  # 32
celkovy_pocet_dlaždic = pocet_dlaždic_y * pocet_dlaždic_x

soubory = sorted(glob.glob(os.path.join(SLOZKA_VSTUP, "*.nc")))
if not soubory:
    print(f"Chyba: Žádné soubory ve složce '{SLOZKA_VSTUP}'!")
    exit()

print(f"Nalezeno {len(soubory)} souborů. Krok 1/2: Generování / Načítání Land Masky...")

# --- KROK 1: VYTVOŘENÍ NEBO NAČTENÍ LAND MASKY ---
# Maska: 1 = pevnina (farmy, pole, lesy), 0 = moře/oceán/Antarktida
if not os.path.exists(SOUBOR_MASKY):
    print("Maska nenalezena. Procházím soubory pro detekci pevniny...")
    globalni_maska = np.zeros((Y_POINTS, X_POINTS), dtype=bool)
    
    for idx, cesta in enumerate(soubory):
        print(f"  Skenuji masku [{idx+1}/{len(soubory)}]: {os.path.basename(cesta)}")
        with xr.open_dataset(cesta) as ds:
            # Přejmenujte 'SPEI3' případně na 'spe06' podle potřeby
            data = np.array(ds["SPEI3"].squeeze().values)[:Y_POINTS, :X_POINTS]
            # Kdekoliv jsou data platná (nejsou NaN a nejsou -9999), je to pevnina
            platna_data = ~np.isnan(data) & (data > -9000)
            if platna_data.ndim == 3: # Pokud má soubor více měsíců
                platna_data = np.any(platna_data, axis=0)
            globalni_maska |= platna_data
            
    # Ořez Antarktidy (pod -60° j.š. zemědělství neprobíhá)
    for lat_idx in range(Y_POINTS):
        lat_val = 89.875 - (lat_idx * 0.25)
        if lat_val < -60.0:
            globalni_maska[lat_idx, :] = False
            
    globalni_maska.astype(np.uint8).tofile(SOUBOR_MASKY)
    print(f"-> Maska pevniny vytvořena a uložena do: {SOUBOR_MASKY}")
else:
    print(f"-> Načítám existující masku pevniny z: {SOUBOR_MASKY}")
    globalni_maska = np.fromfile(SOUBOR_MASKY, dtype=np.uint8).reshape((Y_POINTS, X_POINTS)).astype(bool)

# --- KROK 2: SEKVENČNÍ ZPRACOVÁNÍ DLAŽDIC (Šetrné k HDD i CPU) ---
print("\nKrok 2/2: Extrakce dat a zápis do finální databáze...")

# Nejdříve načteme VŠECHNA časová data do paměti postupně, ale ořezaná o moře, 
# abychom plotnový disk netrápili neustálým otevíráním souborů dokola.
vsechny_mesice_list = []

for idx, cesta in enumerate(soubory):
    print(f" Načítám data do RAM [{idx+1}/{len(soubory)}]: {os.path.basename(cesta)}")
    with xr.open_dataset(cesta) as ds:
        data = np.array(ds["SPEI3"].squeeze().values)[:Y_POINTS, :X_POINTS]
    if data.ndim == 2:
        data = np.expand_dims(data, axis=0)
        
    # Kvantizace z -4.0 až 4.0 na 4 bity (hodnoty 0 až 14). Moře/chyba bude 15.
    data_oříznutá = np.clip(data, -4.0, 4.0)
    kvantovana = ((data_oříznutá + 4.0) * (14.0 / 8.0)).astype(np.uint8)
    
    # Kde je moře podle naší masky, natvrdo hodíme 15
    for m in range(kvantovana.shape[0]):
        kvantovana[m][~globalni_maska] = 15
        
    vsechny_mesice_list.append(kvantovana)

# Spojíme všechny měsíce ze všech souborů nad sebou [Cas, Y, X]
kvantovana_data_komplet = np.concatenate(vsechny_mesice_list, axis=0)
celkovy_pocet_mesicu = kvantovana_data_komplet.shape[0]
print(f"Celkem nalezeno {celkovy_pocet_mesicu} měsíců ke zpracování.")

# Pomocná funkce pro sbalení dvou 4bitových hodnot do jednoho bajtu
def zabal_pole_4bit(pole_uint8):
    if len(pole_uint8) % 2 != 0:
        pole_uint8 = np.append(pole_uint8, np.uint8(15))
    horni = pole_uint8[0::2] << 4
    dolni = pole_uint8[1::2]
    return (horni | dolni).astype(np.uint8).tobytes()

rejstrik = []
aktualni_offset = 0
zpracovane_dlaždice = 0

# Otevíráme finální soubor pro sekvenční zápis (HDD zápis poběží hladce jako po másle)
with open(SOUBOR_DATA, "wb") as f_out:
    for ty in range(pocet_dlaždic_y):
        for tx in range(pocet_dlaždic_x):
            y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
            x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
            
            # Výřez masky a dat pro tuto konkrétní dlaždici
            maska_dlaždice = globalni_maska[y_od:y_do, x_od:x_do]
            
            # Pokud je v celé dlaždici maska samá nula (čisté moře)
            if not np.any(maska_dlaždice):
                rejstrik.append({"offset": aktualni_offset, "delka": 0})
                zpracovane_dlaždice += 1
                continue
                
            # Extrahujeme data pouze pro pixely pevniny napříč všemi měsíci
            # Výsledkem je ploché 1D pole obsahující JEN pevninu řazenou v čase
            dlaždice_data = kvantovana_data_komplet[:, y_od:y_do, x_od:x_do]
            
            # Trik: Maskování ve 3D (přesuneme masku pevniny do všech časů)
            maska_3d = np.broadcast_to(maska_dlaždice, dlaždice_data.shape)
            jen_pevnina_data = dlaždice_data[maska_3d] # Tímto kompletně zmizí moře!
            
            # Sbalíme data do 4 bitů a zkomprimujeme pomocí GZIP (úroveň 9)
            sbalena_data = zabal_pole_4bit(jen_pevnina_data)
            komprimovane = gzip.compress(sbalena_data, compresslevel=9)
            velikost_bloku = len(komprimovane)
            
            # Zápis přímo na disk - žádné skákání hlaviček HDD
            f_out.write(komprimovane)
            
            rejstrik.append({"offset": aktualni_offset, "delka": velikost_bloku})
            aktualni_offset += velikost_bloku
            
            # Výpis procentuálního stavu na obrazovku
            zpracovane_dlaždice += 1
            procenta = (zpracovane_dlaždice / celkovy_pocet_dlaždic) * 100
            print(f" Pokrok: {procenta:.1f}% | Zpracována dlaždice {ty}_{tx} (velikost: {velikost_bloku} B)", end="\r")

# Uložení indexu (rejstříku) pro aplikaci
with open(SOUBOR_REJSTRIK, "w", encoding="utf-8") as f_json:
    json.dump({"pocet_mesicu": celkovy_pocet_mesicu, "dlaždice": rejstrik}, f_json, ensure_ascii=False)

velikost_mb = os.path.getsize(SOUBOR_DATA) / (1024 * 1024)
print(f"\n\nHOTOVO! Výsledný soubor má pouhých {velikost_mb:.2f} MB.")
print(f"Vytvořen rejstřík: {SOUBOR_REJSTRIK} a maska pevniny: {SOUBOR_MASKY}")
