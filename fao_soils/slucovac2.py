
import os
import json
import gzip
import struct

# Konfigurace vstupů a výstupů
VSTUP_RASTR = "HWSD2.bil"
VYSTUP_DATA = "svet_mapa_precizni.gz.bin"
VYSTUP_REJSTRIK = "registrace_mapy.json"

# Nativní parametry FAO HWSD2
NCOLS = 43200
NROWS = 21600
XDIM = 0.00833333333333333
YDIM = 0.00833333333333333
ULX = -179.995833333333
ULY = 89.9958333333333
NODATA = 65535

# Definice dlaždic pro 1km rastr (1080 bodů beze zbytku dělí 43200 i 21600)
# 43200 / 1080 = 40 dlaždic horizontálně
# 21600 / 1080 = 20 dlaždic vertikálně
TILE_SIZE = 1080
DLAZDIC_X = NCOLS // TILE_SIZE # 40
DLAZDIC_Y = NROWS // TILE_SIZE # 20

# Výpočet řádkových limitů pro ořez šířky (-60° až +60°)
# 60° s.š. odpovídá řádku cca 3600, -60° j.š. odpovídá řádku cca 18000
ROW_START = int((ULY - 60.0) / YDIM)
ROW_END = int((ULY - (-60.0)) / YDIM)

print("Zahajuji precizní ořez a dlaždicovou kompresi 2GB rastru...")
print(f"Ořez nastaven na řádky {ROW_START} až {ROW_END} (Zóna -60° až +60° šířky).")

rejstrik_mapy = []
binarni_data = bytearray()
aktualni_offset = 0

# Otevřeme obří soubor pro přímé seekování po jednotlivých bajtech (RAM šetrné)
with open(VSTUP_RASTR, "rb") as input_f:
    
    # Procházíme svět po obřích 1km dlaždicích
    for ty in range(DLAZDIC_Y):
        row_offset_tile = ty * TILE_SIZE
        
        for tx in range(DLAZDIC_X):
            col_offset_tile = tx * TILE_SIZE
            
            # Kontrola, zda dlaždice nespadá celá mimo náš povolený ořez -60 až +60
            if (row_offset_tile + TILE_SIZE < ROW_START) or (row_offset_tile > ROW_END):
                # Dlaždice leží celá v zamrzlých pólech -> vynecháme ji
                rejstrik_mapy.append({"offset": aktualni_offset, "delka": 0})
                continue
            
            # Dočasný zásobník pro načtení surových dat jedné dlaždice z disku
            dlaždice_bity = bytearray()
            je_cisty_ocean = True
            
            # Čteme čtverec 1080x1080 bod po bodu řádek po řádku
            for local_row in range(TILE_SIZE):
                global_row = row_offset_tile + local_row
                
                # Výpočet přesné bajtové pozice začátku řádku této dlaždice v 2GB souboru
                # Každý pixel má 2 bajty
                byte_offset = (global_row * NCOLS + col_offset_tile) * 2
                
                input_f.seek(byte_offset)
                # Načteme celý řádek dlaždice naráz (1080 pixelů * 2 bajty = 2160 bajtů)
                radek_bajty = input_f.read(TILE_SIZE * 2)
                
                if len(radek_bajty) == TILE_SIZE * 2:
                    dlaždice_bity.extend(radek_bajty)
                    
                    # Rychlá kontrola, zda řádek obsahuje pevninu (cokoli jiného než 65535)
                    # Rozbalujeme po 2 bajtech (Little-endian '<H')
                    if je_cisty_ocean:
                        for p_idx in range(TILE_SIZE):
                            val = struct.unpack_interleaved = struct.unpack('<H', radek_bajty[p_idx*2 : (p_idx+1)*2])[0]
                            if val != NODATA:
                                je_cisty_ocean = False
                                break
                else:
                    # Pokud by soubor neočekávaně skončil, vyplníme zbytek maskou nodata
                    dlaždice_bity.extend(struct.pack('<H', NODATA) * TILE_SIZE)

            # VYHODNOCENÍ DLAŽDICE
            if je_cisty_ocean:
                # Dlaždice obsahuje 100% čistý oceán -> do výsledného souboru ji vůbec nezapíšeme
                rejstrik_mapy.append({"offset": aktualni_offset, "delka": 0})
                continue
            
            # Pokud dlaždice obsahuje byť kousek pevniny, zkomprimujeme ji pomocí Gzipu
            komprimovany_blok = gzip.compress(dlaždice_bity, compresslevel=5)
            velikost_bloku = len(komprimovany_blok)
            
            # Zápis metadat do JSON rejstříku
            rejstrik_mapy.append({
                "offset": aktualni_offset,
                "delka": velikost_bloku
            })
            
            # Přidáme zkomprimovaná data do hlavního balíku
            binarni_data.extend(komprimovany_blok)
            aktualni_offset += velikost_bloku

print("\nZápis výsledných souborů na pevný disk...")
with open(VYSTUP_DATA, "wb") as f:
    f.write(binarni_data)

# Rejstřík bude obsahovat i základní geometrické parametry dlaždic pro PHP
meta_rejstrik = {
    "tile_size": TILE_SIZE,
    "dlazdic_x": DLAZDIC_X,
    "dlazdic_y": DLAZDIC_Y,
    "mapa": rejstrik_mapy
}

with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump(meta_rejstrik, f, ensure_ascii=False)

velikost_mb = os.path.getsize(VYSTUP_DATA) / (1024 * 1024)
print("==================================================")
print("ÚSPĚŠNĚ DOKONČENO BEZ PROBLÉMŮ S RAM!")
print(f"Vytvořena 1km precizní mapa světa: {VYSTUP_DATA} ({velikost_mb:.2f} MB)")
print(f"Vytvořen rejstřík dlaždic: {VYSTUP_REJSTRIK}")
print("==================================================")
