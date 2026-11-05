import glob
import os
import json
import gzip
import numpy as np
import xarray as xr

VYSTUP_DATA = "celysvet_rychlost_srazky.gz.bin"
VYSTUP_REJSTRIK = "registrace_rychlost_srazky.json"

VELIKOST_DAVKY_MESICU = 100  # Zpracování po malých časových úsecích pro úsporu RAM
TILE_SIZE = 30               # 30 dělí beze zbytku 360 (12 dlaždic) i 180 (6 dlaždic)

# Automatické vyhledání souboru srážek
nalezeno = (glob.glob("rx1day_*.nc")
            or glob.glob("databaze/rx1day_*.nc")
            or glob.glob("*.nc")
            or glob.glob("databaze/*.nc"))

if not nalezeno:
    print("Chyba: Nebyl nalezen žádný soubor .nc!")
    raise SystemExit(1)

SOUBOR_NC = sorted(nalezeno)[0]
print(f"Zahajuji ŠETRNOU DAVKOVOU konverzi srážek: {os.path.basename(SOUBOR_NC)}")

with xr.open_dataset(SOUBOR_NC) as ds:
    celkovy_pocet_mesicu = int(ds.sizes['time'])
    NROWS = int(ds.sizes['lat'])
    NCOLS = int(ds.sizes['lon'])
    print(f"Detekované rozlišení: {NCOLS}x{NROWS} bodů, měsíců: {celkovy_pocet_mesicu}")
    print(f"Jednotky rx1day v souboru: {ds['rx1day'].attrs.get('units')}")
    matice_srazky = np.zeros((NROWS, NCOLS, celkovy_pocet_mesicu), dtype=np.uint8)

# --- DÁVKOVÉ ZPRACOVÁNÍ PO ČASOVÉ OSE ---
for start_m in range(0, celkovy_pocet_mesicu, VELIKOST_DAVKY_MESICU):
    konec_m = min(celkovy_pocet_mesicu, start_m + VELIKOST_DAVKY_MESICU)

    with xr.open_dataset(SOUBOR_NC) as ds:
        # Nejvyšší hodnota z 21 členů (member)
        data_array = ds["rx1day"].isel(time=slice(start_m, konec_m)).max(dim="member").squeeze()
        if data_array.dims != ('time', 'lat', 'lon'):
            data_array = data_array.transpose('time', 'lat', 'lon')
        # rx1day je už v mm, žádný přepočet (×86400) se nedělá
        mm_davka = np.array(data_array.values, dtype=np.float32)

    # Neplatné hodnoty (oceány, masky, fill hodnoty)
    maska_nan = np.isnan(mm_davka) | (mm_davka < 0.0) | (mm_davka > 1e10)
    mm_cista = np.where(maska_nan, 0.0, mm_davka)

    # Kvantování: 0–250 = mm, 251–254 = pásma extrémů, 255 = maska
    bajtova_davka = np.minimum(np.round(mm_cista), 250.0)
    bajtova_davka = np.where((mm_cista > 250.0) & (mm_cista <= 300.0), 251, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 300.0) & (mm_cista <= 400.0), 252, bajtova_davka)
    bajtova_davka = np.where((mm_cista > 400.0) & (mm_cista <= 500.0), 253, bajtova_davka)
    bajtova_davka = np.where(mm_cista > 500.0, 254, bajtova_davka)
    bajtova_davka = np.where(maska_nan, 255, bajtova_davka).astype(np.uint8)

    # Překlopení severního pólu pro Leaflet
    bajtova_davka = np.flip(bajtova_davka, axis=1)

    # (čas, lat, lon) -> (lat, lon, čas)
    matice_srazky[:, :, start_m:konec_m] = np.transpose(bajtova_davka, (1, 2, 0))

    procento = (konec_m / celkovy_pocet_mesicu) * 100
    print(f"  -> Zpracováno: {konec_m}/{celkovy_pocet_mesicu} měsíců ({procento:.1f} %)")

print("Dávkové zpracování dokončeno. Balím výslednou matici do dlaždic...")

# --- GZIP DLAŽDICE (čisté oceány se vynechají) ---
rejstrik = []
binarni_kontainer = bytearray()
aktualni_offset = 0

for ty in range(NROWS // TILE_SIZE):
    for tx in range(NCOLS // TILE_SIZE):
        y_od, y_do = ty * TILE_SIZE, (ty + 1) * TILE_SIZE
        x_od, x_do = tx * TILE_SIZE, (tx + 1) * TILE_SIZE
        dlazdice = matice_srazky[y_od:y_do, x_od:x_do, :]

        if np.all(dlazdice == 255):
            rejstrik.append({"offset": 0, "delka": 0})
            continue

        komprimovane = gzip.compress(dlazdice.tobytes(), compresslevel=6)
        rejstrik.append({"offset": aktualni_offset, "delka": len(komprimovane)})
        binarni_kontainer.extend(komprimovane)
        aktualni_offset += len(komprimovane)

with open(VYSTUP_DATA, "wb") as f:
    f.write(binarni_kontainer)

with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f:
    json.dump({
        "pocet_mesicu": celkovy_pocet_mesicu,
        "nrows": NROWS,
        "ncols": NCOLS,
        "tile_size": TILE_SIZE,
        "dlaždice": rejstrik
    }, f, ensure_ascii=False)

# Kontrola: rozložení kódů (po opravě by 254 neměl dominovat)
hodnoty, pocty = np.unique(matice_srazky, return_counts=True)
print("Nejčastější kódy:", sorted(zip(pocty.tolist(), hodnoty.tolist()), reverse=True)[:6])
print("==================================================")
print("DATABÁZE SRÁŽEK ÚSPĚŠNĚ SESTAVENA!")
print(f"Soubor: {VYSTUP_DATA} ({os.path.getsize(VYSTUP_DATA)/(1024*1024):.2f} MB)")
print(f"Počet měsíců: {celkovy_pocet_mesicu}  (v HTML musí být this.mesicu přesně stejné)")
print("==================================================")
