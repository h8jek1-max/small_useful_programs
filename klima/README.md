# Agroclimatic indicators and climatic agricultural potential

A web application (Leaflet map, chart, table export) built on a database of Copernicus (C3S) agroclimatic indicators in
10-day steps for the years 2011–2099, plus a simple estimate of "climatic agricultural potential" based on a table of crop species.
Global land (60° S to 90° N), 0.5° × 0.5° grid, 67,420 cells, finished database about 180 MB.

The database is based on the agro-climatic indicators Copernicus (C3S) [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=overview) experiment RCP6.0 (According to the methodology of the Intergovernmental Panel on Climate Change (IPCC) it is an older medium-high scenario, in which greenhouse gas emissions peak around 2080 and subsequently decrease. The overall radiation effect will thus be stabilized after 2100 to the limit of 6.0 W/m2.) in a grid of 720 × 360 points with a step of 0.5°, in a time step of ten days for 2011 to 2099. It contains 67 420 cells of land.

Created with the help of AI Claude version for free.

> The potential is an indicative estimate based on climatic limits of species. It does not include evapotranspiration, soil, pests
> or management and must not be used as a yield forecast.

## Repository contents

| file | purpose |
|---|---|
| `odhad_indikatoru.py`, `odhad_vq2.py` | netCDF reading and shared functions (value warping, codebooks, annual-total distribution) |
| `odhad_vq3.py` | **optional** estimate of accuracy and size of storage variants (documents why the database is stored the way it is) |
| `stavba_db.py` | **main database build**, batched and resumable (netCDF → `databaze/`) |
| `balit_dlazdice.py` | repacks into 30 × 30 cell tiles, three parts per tile → `databaze/web/` |
| `tnn_entropie.py` | **optional** measurement of whether the decadal TNn could be stored more compactly |
| `klima_agro.html` | the page (map, chart, table, settings bubbles) |
| `klima_dekoder.js`, `klima_agro.js` | in-browser tile decoding and calculations (parameters, statistics, potential); no libraries |
| `get_klima_post.php` | server only sends packed tile parts according to an index (it decodes nothing) |
| `test_dekoder.js` | **optional** test of the JS decoder against Python (needs Node.js) |

## Requirements

- Python 3 with `numpy` and `netCDF4` (`pip install numpy netCDF4`)
- about 4.5 GB free disk space during the build (temporary store `tmp_*`, can be deleted afterwards)
- a server with PHP to run (any web server that can run the PHP script `get_klima_post.php`)
- Node.js is not required, only for the optional test

## Input data

1. **Copernicus (C3S) agroclimatic indicators**, netCDF, 720 × 360 grid with 0.5° spacing, 10-day time step
   (3 files per indicator: 2011–2040, 2041–2070, 2071–2099): `FD`, `TR`, `RR1`, `R20mm`, `SU` (thresholds 20, 25, 30, 35 °C),
   `TG`, `TN`, `TX`, `TNn`, `TXx`, `SDII`. Scenario and model: *(fill in according to the downloaded dataset)*.
   File paths and variable names are set in `odhad_indikatoru.py`.

## Reproduction steps

Order and approximate times (low-end computer, whole world):

```text
1. python stavba_db.py                          # build, resumes after interruption
       phase A  netCDF -> temporary store         ~12 min
       phase B  codebook training                 ~6 min
       phase C  batched encoding                  ~5 h 20 min
       phase D  writes meta.json, bunky.bin etc.
2. python stavba_db.py --overit                 # check from the finished files (~20 s)
3. python balit_dlazdice.py                     # tiles and index                           ~1-2 min
5. python priprav_druhy.py sloucena_ecocrop.csv # -> druhy_klima.json
```

Useful `stavba_db.py` switches: `--max-minut 60` (stops after an hour and continues next time), `--davka 1` (less memory),
`--krok 0.5 --strop 10` (decadal TNn: step in °C and the ceiling above which TNn is not stored exactly), `--sdii-prah 50`
(precipitation intensity from this value is stored exactly), `--smazat-temp` (deletes the temporary store), `--znovu` (start from scratch),
`--znovu-kodovani` (repeat phase C only).

### What the build produces

`databaze/web/` holds what the server needs:

- `celysvet_klima.gz.bin`: tiles; each has three parts (`rocni` ≈ 32 MB in total, `model` ≈ 56 MB, `tnn` ≈ 89 MB)
- `registrace_klima.json`: index (part offsets), land mask, decoder parameters
- `kody_klima.bin`: codebooks (about 1.9 MB, downloaded once)

The page downloads only the parts needed for the chosen view (e.g. a map of exact annual values needs only the `rocni` part).

### Deployment

Put into one server folder: `klima_agro.html`, Leaflet (`leaflet.min.js`, `leaflet.min.css`), `klima_dekoder.js`, `klima_agro.js`,
`druhy_klima.json`, `registrace_klima.json`, `celysvet_klima.gz.bin`, `kody_klima.bin`, `get_klima_post.php`.

## How the data are stored (briefly)

For better compression, values are "warped" by a piecewise linear mapping. For every cell the typical annual course in decades is stored
(a two-stage codebook), for every month and year a number from a codebook of deviations (9 bits) and the monthly spread. Stored exactly:
annual totals of day counts, annual mean temperature (0.5 °C step), annual maximum temperature, and the dates of the last spring and
first autumn frost. Decadal TNn is stored as a residual against the estimate from the codes (zigzag, 0.5 °C step up to +10 °C). SDII
from 50 mm/day is stored exactly as exceptions. On decoding, annual totals are distributed back into decades.

## Accuracy (check on 50 cells from around the world)

| indicator | RMSE | P95 of error |
|---|---|---|
| temperatures TG, TN, TX | 1.2 to 1.6 °C | 2.5 to 3.4 °C |
| TNn (up to +10 °C) | 1.15 °C | 2.7 °C |
| decadal TXx | 2.8 °C | 5.8 °C |
| day counts per decade (FD, TR, SU, RR1, R20mm) | 0.4 to 1.5 days | 1 to 4 days |
| SDII | 3.1 mm/day | 5.8 mm/day |

Annual totals of day counts, annual TG, annual TXx and frost dates are exact. The data contain an exact SDII ceiling of 400 mm/day
(9,131 decades), which probably comes from the source.

## Limitations

- Values per decade (except TNn up to +10 °C) are estimates from the codes and are marked ≈ in the application.
- Total precipitation is not in the database; the application estimates it as number of wet days × intensity (SDII).
- Crop potential: see the warning above. A species' weight is kcal per 100 g, species without calories are not counted.
- `celysvet_klima.gz.bin` is about 180 MB, above GitHub's 100 MB per-file limit. Use Git LFS or attach it to a release (Releases)
  and keep only the scripts in the repository.

 The species table is a compilation made with AI search Google assistance and has not been expertly reviewed.


czech: 

# Agroklimatické ukazatele a klimatický zemědělský potenciál

Webová aplikace (Leaflet mapa, graf, výpis) nad databází agroklimatických indikátorů Copernicus (C3S) v dekádách
pro roky 2011 až 2099 a jednoduchý odhad „klimatického zemědělského potenciálu“ podle tabulky druhů plodin.
Celý svět (pevnina mezi 60° j. š. a 90° s. š.), mřížka 0,5° × 0,5°, 67 420 buněk, hotová databáze má asi 180 MB.

Databáze vychází z agroklimatických indikátorů Copernicus (C3S) [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=overview) experiment RCP6.0 (Podle metodiky Mezivládního panelu pro změnu klimatu (IPCC) jde o starší středně vysoký scénář, ve kterém emise skleníkových plynů vrcholí kolem roku 2080 a následně klesají. Celkové radiační působení se tak po roce 2100 ustálí na hranici 6,0 W/m².) v mřížce 720 × 360 bodů s krokem 0,5°, v časovém kroku po deseti dnech pro roky 2011 až 2099. Obsahuje 67 420 buněk pevniny. 

Vytvořeno s pomocí AI Claude verzí zdarma.

> Potenciál je orientační odhad podle klimatických mezí druhů. Nezahrnuje například evapotranspiraci, půdu, škůdce ani
> hospodaření a nemá se brát jako předpověď úrody.

## Co v repozitáři je

| soubor | k čemu slouží |
|---|---|
| `odhad_indikatoru.py`, `odhad_vq2.py` | čtení netCDF dat a sdílené funkce (pokřivení hodnot, kódové knihy, rozdělení ročních součtů) |
| `odhad_vq3.py` | **volitelný** odhad přesnosti a velikosti variant ukládání (v repozitáři dokumentuje, proč je databáze ukládána tak, jak je) |
| `stavba_db.py` | **hlavní stavba databáze** po dávkách s pokračováním (netCDF → složka `databaze/`) |
| `balit_dlazdice.py` | přebalení do dlaždic 30 × 30 buněk, tři části na dlaždici → `databaze/web/` |
| `tnn_entropie.py` | **volitelné** měření, zda by šlo TNn po dekádách uložit úsporněji |
| `klima_agro.html` | stránka (mapa, graf, výpis, bubliny nastavení) |
| `klima_dekoder.js`, `klima_agro.js` | dekódování dlaždic v prohlížeči a výpočty (parametry, statistiky, potenciál); bez knihoven |
| `get_klima_post.php` | server posílá jen zabalené části dlaždic podle rejstříku (nic nedekóduje) |
| `test_dekoder.js` | **volitelný** test JS dekodéru proti Pythonu (potřebuje Node.js) |

## Požadavky

- Python 3 s balíčky `numpy` a `netCDF4` (`pip install numpy netCDF4`)
- na disku asi 4,5 GB volného místa při stavbě (dočasné úložiště `tmp_*`, po stavbě jde smazat)
- server s PHP pro provoz (jinak stačí soubory servírovat přes jakýkoli webový server, který umí PHP skript `get_klima_post.php`)
- Node.js není potřeba, jen pro volitelný test

## Vstupní data

1. **Agroklimatické indikátory Copernicus (C3S)**, netCDF, v mřížce 720 × 360 s krokem 0,5°, v časovém kroku 10 dnů
   (3 soubory na indikátor: 2011 až 2040, 2041 až 2070, 2071 až 2099): `FD`, `TR`, `RR1`, `R20mm`, `SU` (prahy 20, 25, 30, 35 °C),
   `TG`, `TN`, `TX`, `TNn`, `TXx`, `SDII`.
   Cesty k souborům a názvy proměnných se nastavují v `odhad_indikatoru.py`.
2. **FAO Ecocrop** (teplotní a srážkové meze, délka cyklu): uložené HTML stránky `datasheet_*.html` a `view_*.html`.
3. **Tabulka druhů** `sloucena_tabulka.csv` (USDA zóna, BEDD_C, tolerance sucha, kalorie, délka vegetační doby a další).
   Byla sestavena ručně s pomocí AI vyhledávače, proto hodnoty před použitím ověřte u zdrojů.

## Postup reprodukce

Pořadí kroků a orientační časy (slabý počítač, celý svět):

```text
1. python stavba_db.py                          # stavba, pokračuje po přerušení
       fáze A  netCDF -> dočasné úložiště         ~12 min
       fáze B  výcvik kódových knih               ~6 min
       fáze C  kódování po dávkách                ~5 h 20 min
       fáze D  zápis meta.json, bunky.bin atd.
2. python stavba_db.py --overit                 # kontrola z hotových souborů (~20 s)
3. python balit_dlazdice.py                     # dlaždice a rejstřík                      ~1 až 2 min
```

Užitečné přepínače `stavba_db.py`: `--max-minut 60` (po hodině se zastaví a příště pokračuje), `--davka 1` (méně paměti),
`--krok 0.5 --strop 10` (TNn po dekádách: krok v °C a strop, nad kterým se TNn neukládá přesně), `--sdii-prah 50`
(intenzita srážek od této hodnoty se ukládá přesně), `--smazat-temp` (smaže dočasné úložiště), `--znovu` (začne od nuly),
`--znovu-kodovani` (zopakuje jen fázi C).

### Co stavba vyrobí

`databaze/web/` obsahuje to, co je potřeba na serveru:

- `celysvet_klima.gz.bin`: dlaždice; každá má tři části (`rocni` ≈ 32 MB celkem, `model` ≈ 56 MB, `tnn` ≈ 89 MB)
- `registrace_klima.json`: rejstřík (offsety částí), maska pevniny, parametry dekodéru
- `kody_klima.bin`: kódové knihy (asi 1,9 MB, stahuje se jednou)

Stránka stahuje jen části, které uživatel pro zvolené zobrazení potřebuje (např. mapa přesných ročních hodnot jen část `rocni`).

### Nasazení

Do jedné složky na serveru: `klima_agro.html`, Leaflet (`leaflet.min.js`, `leaflet.min.css`), `klima_dekoder.js`, `klima_agro.js`,
`druhy_klima.json`, `registrace_klima.json`, `celysvet_klima.gz.bin`, `kody_klima.bin`, `get_klima_post.php`.

## Jak jsou data uložena (stručně)

Hodnoty se pro lepší kompresi „pokřiví“ po částech lineárním převodem. Pro každou buňku se uloží typický roční průběh v dekádách
(kódová kniha ve dvou stupních), pro každý měsíc a rok číslo z kódové knihy odchylek (9 bitů) a rozptyl měsíců. Přesně se ukládají
roční součty počtů dní, roční průměrná teplota (krok 0,5 °C), roční nejvyšší teplota a termíny posledního jarního a prvního
podzimního mrazu. TNn po dekádách se ukládá jako zbytek proti odhadu z kódů (zigzag, krok 0,5 °C do +10 °C). SDII od 50 mm/den
se ukládá přesně jako výjimky. Při dekódování se roční součty rozdělují zpět do dekád.

## Přesnost (kontrola na 50 buňkách z celého světa)

| indikátor | RMSE | P95 chyby |
|---|---|---|
| teploty TG, TN, TX | 1,2 až 1,6 °C | 2,5 až 3,4 °C |
| TNn (do +10 °C) | 1,15 °C | 2,7 °C |
| TXx po dekádách | 2,8 °C | 5,8 °C |
| počty dní po dekádách (FD, TR, SU, RR1, R20mm) | 0,4 až 1,5 dne | 1 až 4 dny |
| SDII | 3,1 mm/den | 5,8 mm/den |

Roční součty počtů dní, roční TG, roční TXx a termíny mrazů jsou přesné. V datech je přesný strop SDII 400 mm/den
(9 131 dekád), který pravděpodobně pochází ze zdroje.

## Omezení

- Hodnoty po dekádách (kromě TNn do +10 °C) jsou odhady z kódů, v aplikaci jsou označeny ≈.
- Úhrn srážek není v databázi; aplikace ho odhaduje jako počet deštivých dnů × intenzita (SDII).
- Potenciál plodin: viz upozornění nahoře. Váha druhu je kcal na 100 g, druhy bez kalorií se nezapočítají.
- Soubor `celysvet_klima.gz.bin` má kolem 180 MB, což je nad limitem GitHubu 100 MB na jeden soubor. Použijte Git LFS nebo ho
  přiložte k vydání (Releases) a v repozitáři ponechte jen skripty.

Tabulka druhů je kompilace s AI vyhledávače Google pomocí a neprošla odbornou kontrolou.

