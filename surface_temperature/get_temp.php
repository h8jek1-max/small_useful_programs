
<?php
header('Access-Control-Allow-Origin: *');
header('Content-Type: application/json; charset=utf-8');

// Kontrola hromadných parametrů z viditelného okna mapy
if (!isset($_GET['lats']) || !isset($_GET['lons'])) {
    http_response_code(400);
    echo json_encode(["error" => "Chybí hromadné parametry lats/lons"]);
    exit;
}

$lats = explode(',', $_GET['lats']);
$lons = explode(',', $_GET['lons']);
$pocet_bodu = count($lats);

if ($pocet_bodu !== count($lons)) {
    http_response_code(400);
    echo json_encode(["error" => "Počet Lat a Lon souřadnic nesouhlasí"]);
    exit;
}

// Nativní parametry Gaussovy mřížky teplotního modelu (256x128)
$ncols = 256; $nrows = 128;
$xdim = 1.40625; $ydim = 1.40625; // Krok sítě v stupních
$ulx = 0.0; $uly = 88.92774;      // Počáteční severozápadní body modelu
$tile_size = 32; $dlazdic_x = 8;
$celkem_mesicu = 1020; // 2015 - 2099 (12 měsíců v roce)

// Načtení společného JSON rejstříku teplotních dlaždic
$rejstrik_temp = json_decode(file_get_contents("registrace_teploty.json"), true);
$mapa_dlazdic = $rejstrik_temp['dlaždice'];

// Pomocná paměť (Cache) v RAM skriptu pro bleskové hromadné čtení
$cache = [];
$vysledky_balik = [];

// Otevřeme binární databázi pouze JEDNOU pro celý požadavek
$f_temp = fopen("celysvet_teplota.gz.bin", "rb");

for ($i = 0; $i < $pocet_bodu; $i++) {
    $lat = (float)$lats[$i];
    $lon = (float)$lons[$i];
    
    // Klíč pro přesné spárování s sub-pixely na frontendu
    $bod_klic = "{$lats[$i]}_{$lons[$i]}";

    // Ošetření záporných délek (NetCDF model běží v rozsahu 0 až 358.59 stupňů východně)
    if ($lon < 0) {
        $lon += 360.0;
    }

    // Výpočet pozice v Gaussově mřížce 256x128
    $col = (int)round(($lon - $ulx) / $xdim);
    $row = (int)round(($uly - $lat) / $ydim);
    $col = max(0, min($ncols - 1, $col));
    $row = max(0, min($nrows - 1, $row));

    // Zjištění indexu souborové dlaždice 32x32
    $tx = (int)($col / $tile_size);
    $ty = (int)($row / $tile_size);
    $tile_index = ($ty * $dlazdic_x) + $tx;

    if (!isset($mapa_dlazdic[$tile_index]) || $mapa_dlazdic[$tile_index]['delka'] == 0) {
        $vysledky_balik[$bod_klic] = ["error" => "Oceán"];
        continue;
    }

    // Načtení a dekomprese dlaždice do paměti RAM (pokud tam ještě není)
    if (!isset($cache[$tile_index])) {
        fseek($f_temp, $mapa_dlazdic[$tile_index]['offset']);
        $gzip_data = fread($f_temp, $mapa_dlazdic[$tile_index]['delka']);
        $cache[$tile_index] = gzdecode($gzip_data);
    }

    $raw_tile = $cache[$tile_index];
    $local_col = $col % $tile_size;
    $local_row = $row % $tile_size;
    
    // Výpočet přesného offsetu časové osy (každý pixel drží 1020 bajtů za sebou)
    $pixel_bytes_offset = ($local_row * $tile_size + $local_col) * $celkem_mesicu;
    $surove_teploty_string = substr($raw_tile, $pixel_bytes_offset, $celkem_mesicu);
    
    // ÚPRAVA PRO PŮLSTUPŇOVOU DEKOMPRESI:
    $historie_celsia = [];
    for ($j = 0; $j < $celkem_mesicu; $j++) {
        $kod_baju = ord($surove_teploty_string[$j]);
        if ($kod_baju === 255) {
            $historie_celsia[] = null;
        } else {
            // Odečteme posun 100 a vydělíme dvěma, čímž získáme přesný float na 0.5 °C
            $historie_celsia[] = ($kod_baju - 100) / 2.0;
        }
    }


    $vysledky_balik[$bod_klic] = [
        "smu_id" => $tile_index,
        "data" => $historie_celsia
    ];
}

fclose($f_temp);

// Odešleme hromadný balík všech souřadnic prohlížeči najednou v JSON formátu
echo json_encode($vysledky_balik);
exit;
