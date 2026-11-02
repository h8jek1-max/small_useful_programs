import os
import json
import gzip
import struct

VSTUP_RASTR = "HWSD2.bil"
VYSTUP_DATA = "svet_mapa_precizni.gz.bin"
VYSTUP_REJSTRIK = "registrace_mapy.json"

NCOLS = 43200
NROWS = 21600
XDIM = 0.00833333333333333; YDIM = 0.00833333333333333
ULX = -179.995833333333; ULY = 89.9958333333333
NODATA = 65535

TILE_SIZE = 1080
DLAZDIC_X = NCOLS // TILE_SIZE # 40
DLAZDIC_Y = NROWS // TILE_SIZE # 20

# !!! KLÍČOVÁ ZMĚNA: ZAČÍNÁME OD INDEXU 0 (SEVERNÍ PÓL +90°) !!!
ROW_START = 0 
ROW_END = int((ULY - (-60.0)) / YDIM) # Skenujeme dolů k Antarktidě

print("Zahajuji precizní ořez 1km rastru půd pro KOMPLETNÍ SEVERNÍ POLOKOULI...")

rejstrik_mapy = []
binarni_data = bytearray()
aktualni_offset = 0

with open(VSTUP_RASTR, "rb") as input_f:
    for ty in range(DLAZDIC_Y):
        row_offset_tile = ty * TILE_SIZE
        
        for tx in range(DLAZDIC_X):
            col_offset_tile = tx * TILE_SIZE
            
            # Kontrola spodního ořezu Antarktidy
            if row_offset_tile > ROW_END:
                rejstrik_mapy.append({"offset": aktualni_offset, "delka": 0})
                continue
            
            dlaždice_bity = bytearray()
            je_cisty_ocean = True
            
            for local_row in range(TILE_SIZE):
                global_row = row_offset_tile + local_row
                byte_offset = (global_row * NCOLS + col_offset_tile) * 2
                
                input_f.seek(byte_offset)
                radek_bajty = input_f.read(TILE_SIZE * 2)
                
                if len(radek_bajty) == TILE_SIZE * 2:
                    dlaždice_bity.extend(radek_bajty)
                    if je_cisty_ocean:
                        for p_idx in range(TILE_SIZE):
                            val = struct.unpack('<H', radek_bajty[p_idx*2 : (p_idx+1)*2])[0]
                            if val != NODATA:
                                je_cisty_ocean = False
                                break
                else:
                    dlaždice_bity.extend(struct.pack('<H', NODATA) * TILE_SIZE)

            if je_cisty_ocean:
                rejstrik_mapy.append({"offset": aktualni_offset, "delka": 0})
                continue
            
            komprimovany_blok = gzip.compress(dlaždice_bity, compresslevel=5)
            velikost_bloku = len(komprimovany_blok)
            
            rejstrik_mapy.append({"offset": aktualni_offset, "delka": velikost_bloku})
            binarni_data.extend(komprimovany_blok)
            aktualni_offset += velikost_bloku

print("\nZapisuji výslednou 1km severskou mapu půd na disk...")
with open(VYSTUP_DATA, "wb") as f: f.write(binarni_data)

with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({
        "tile_size": TILE_SIZE, "dlazdic_x": DLAZDIC_X, "dlazdic_y": DLAZDIC_Y, "mapa": rejstrik_mapy
    }, f, ensure_ascii=False)

print(f"HOTOVO! 1km severská mapa půd úspěšně uložena: {os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB")
