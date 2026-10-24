import glob
import os
import numpy as np
import xarray as xr

VYSTUP_DATA = "celysvet_evapo.raw.bin"

NCOLS = 256
NROWS = 128
CELKEM_MESICU = 1020

nalezeno = glob.glob("evspsbl_*.nc")
if not list(nalezeno) and os.path.exists("databaze"):
    nalezeno = glob.glob("databaze/evspsbl_*.nc")

if not nalezeno:
    print("Chyba: Nebyl nalezen žádný soubor s maskou 'evspsbl_*.nc'!")
    exit()

SOUBOR_NC = list(nalezeno)[0]
print(f"Zahajuji ZERO-COPY nelineární konverzi celého světa ze souboru: {os.path.basename(SOUBOR_NC)}")

with xr.open_dataset(SOUBOR_NC) as ds:
    data_array = ds["evspsbl"].squeeze()
    if data_array.dims != ('time', 'lat', 'lon'):
        data_array = data_array.transpose('time', 'lat', 'lon')
    surove_kg = np.array(data_array.values)

print(f"Detekováno {surove_kg.shape} měsíců. Přepočítávám na mm/měsíc s progresivní mřížkou...")

sekund_v_mesici = 30.4375 * 86400.0
mm_za_mesic = surove_kg * sekund_v_mesici

# Ošetření neplatných NaN hodnot na masku 255
maska_nan = np.isnan(mm_za_mesic) | (mm_za_mesic < 0.0)
mm_cista = np.where(maska_nan, 0.0, mm_za_mesic)

# Funkce pro nelineární zakódování jednoho čísla do 1 bajtu (0-254)
def zakoduj_hodnotu_evapo(val):
    if val <= 240.0:
        return int(round(val))
    elif val <= 255.0: return 241
    elif val <= 275.0: return 242
    elif val <= 300.0: return 243
    elif val <= 330.0: return 244
    elif val <= 365.0: return 245
    elif val <= 405.0: return 246
    elif val <= 450.0: return 247
    elif val <= 500.0: return 248
    elif val <= 550.0: return 249
    elif val <= 600.0: return 250
    else: return 254  # Extrémní strop nad 600 mm

# Vektorizace funkce pro bleskový chod nad celým numpy polem
zakoduj_vektorově = np.vectorize(zakoduj_hodnotu_evapo)
bajtova_data = zakoduj_vektorově(mm_cista).astype(np.uint8)

# Pokud bylo původně NaN, zapíšeme natvrdo systémovou masku 255
bajtova_data = np.where(maska_nan, 255, bajtova_data).astype(np.uint8)

# GEOGRAFICKÁ KOREKCE SEVERU: Překlopíme osu Y, aby index 0 začínal na severním pólu
bajtova_data = np.flip(bajtova_data, axis=1)

# Překlopíme na tvar (Lat, Lon, Čas) pro lineární historii pixelu
matice_final = np.transpose(bajtova_data, (1, 2, 0))

print("Zapisuji nelineární binární stream na disk...")
with open(VYSTUP_DATA, "wb") as f:
    f.write(matice_final.tobytes())

velikost_souboru = os.path.getsize(VYSTUP_DATA)
ocekavana_velikost = NROWS * NCOLS * CELKEM_MESICU

print("==================================================")
print("ZERO-COPY PROGRESIVNÍ DATABÁZE ÚSPĚŠNĚ VYTVOŘENA!")
print(f"Cílový soubor: {VYSTUP_DATA} ({velikost_souboru / (1024*1024):.2f} MB)")
if velikost_souboru == ocekavana_velikost:
    print("► Kontrola integrity: OK (Matice sedí přesně na bajt!)")
print("==================================================")
