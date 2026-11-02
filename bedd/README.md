
Biologically efficient temperature and precipitation, forecast from 2011 to 2099

 The database comes from the agro-climatic dataset Copernicus (C3S) The spatial grid has a dimension of 720 × 360 points with a step of 0.5°. The time step is a decade (after ten days). The moon is divided into 3 decades (5th, 15th and 25th days of the month). 
 
   BEDD (Biologically Effective Degree Days) 
   The BEDD sum is the heat usable by plants in selected temperature bands. A biologically effective sum of temperatures in ten days. It calculates accumulated heat adjusted for biological crop limits. The model shows 6 different biological scenarios (summations above 0 °C, 2 °C, 4 °C, 6 °C, 8 °C and 10 °C). It compares with the demands of the crop or variety (how much heat it needs to maturity); they are indicated by the breeders. 
   
    RR (Precipitation) 
    Total rainfall in ten days in millimeters (mm). 
    
    The growth time of vegetation
      This data is not part of the Copernicus database, but a rough estimate calculated from their data. This calculation seeks to express the percentage state to the total maximum productivity of vegetation. From the available data (RR and BEDD per decade, plus width and date) it can be calculated: Ra radiation is calculated exactly from width and date. The average temperature is based on the BEDD0 divided by the number of days. PET is calculated by Oudin's formula from Ra and temperature. Snow (daily melting) and Ks according to FAO-56.   

 Method of compression of the database (Original size 4Gb, compressed 140Mb): ranges of values h1/level bdd 0c from the value when the subsequent level of h2 is no longer zero to the value when h2 has a difference with h3 over 20, the others can be inferred without coefficients, according to the rule above the coefficient begins to be calculated, it is calculated to the level where it begins zero or to the level from which the difference begins. Your addition: The rule of “difference 20”: applies only to 10-day decades. The third decade of the month has 8 to 11 days, so the script divides the number of days (takes it from time_bnds). Cropping on 360: H1.1.5 trimmed hot spots (Sahara, Arabia) where a decade will give over 400 °C·days. H1 is a byte with a step of 0.16 °C/day. The level 2 to 10 °C calculates JS from the table. The table has one line per H1 code (251 lines) and 4 width bands. This is the most subtle possible division. Precipitation as 1 byte on a scale with a step of max(0.25 mm; 5% of value), 255 = missing. The scale is in JSON as a table, so it can be tuned (RR_REL) without changing the JS. The database has a grid of 0.5°, decades since 2011 and according to the name of the BEDD file it is a model of GFDL-ESM2M with RCP6.0. 2. 

 Source: Copernicus C3S, model GFDL-ESM2M, scenario RCP6.0, grid 0.5°, decade step. Vegetation, PET and moisture balance are indicative estimates.
        
Database extracted from the site[cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/) .



Biologicky efektivní suma teplot a srážky, předpověď od roku 2011 do 2099

Databáze pochází z agroklimatického datasetu Copernicus (C3S) Prostorová mřížka má rozměr 720 × 360 bodů s krokem 0,5°. Časový krok je dekádní (po deseti dnech). Měsíc je rozdělen na 3 dekády (5., 15. a 25. den v měsíci). 

BEDD (Biologically Effective Degree Days)
Suma BEDD je teplo využitelné rostlinami ve vybraných teplotních pásmech. Biologicky efektivní suma teplot za deset dní. Počítá akumulované teplo upravené o biologické limity plodin. Model ukazuje 6 různých biologických scénářů (sumace nad 0 °C, 2 °C, 4 °C, 6 °C, 8 °C a 10 °C). Porovnává se s nároky plodiny nebo odrůdy (kolik tepla potřebuje do zralosti); udávají je šlechtitelé. 
RR (Precipitation)
Celkový srážkový úhrn za deset dní v milimetrech (mm).  
Doba růstu vegetace
Tyto data nejsou součástí databáze Copernicus, ale hrubý odhad vypočítaný z jejich dat. Tento výpočet se snaží o vyjádření procentního stavu k celkové maximální produktivity vegetace. Z dostupných dat (RR a BEDD za dekádu, plus šířka a datum) se dá spočítat: Záření Ra se spočítá přesně ze šířky a data. Průměrná teplota vychází z BEDD0 děleného počtem dní. PET se počítá Oudinovým vzorcem z Ra a teploty. Sníh (stupňodenní tání) a Ks podle FAO-56. 

Způsob komprimace databáze (Původní velikost 4Gb, komprimovaná 140Mb): 
 rozsahy hodnot h1/hladina bdd 0c od hodnoty kdy následná úroveň h2 již není nula do hodnoty kdy h2 má rozdíl s h3 přes 20, ostatní se dají vyvodit bez koeficientů, podle pravidla výše se začne vypočítávat koeficient, počítá se k hladině kde začíná nula nebo k hladině od které začíná rozdíl snižování nad 20. Vaše doplnění: Pravidlo „rozdíl 20“: platí jen pro 10denní dekády. Třetí dekáda měsíce má 8 až 11 dní, takže skript dělí počtem dní (bere ho z time_bnds). Ořez na 360: H1/1,5 ořezávalo horká místa (Sahara, Arábie), kde dekáda dá přes 400 °C·dní. H1 je bajt s krokem 0,16 °C/den. Hladiny 2 až 10 °C dopočítá JS z tabulky. Tabulka má jeden řádek na každý kód H1 (251 řádků) a 4 šířková pásma. To je nejjemnější možné dělení. Srážky jako 1 bajt na stupnici s krokem max(0,25 mm; 5 % hodnoty), 255 = chybí. Stupnice je v JSONu jako tabulka, takže ji lze ladit (RR_REL) bez změny JS. Databáze má mřížku 0,5°, dekády od 2011 a podle názvu BEDD souboru je to model GFDL-ESM2M s RCP6.0. 2.

Zdroj: Copernicus C3S, model GFDL-ESM2M, scénář RCP6.0, mřížka 0,5°, dekádní krok. Vegetace, PET a vláhová bilance jsou orientační odhady.
        
Databáze stažena ze stránek [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu/) .
