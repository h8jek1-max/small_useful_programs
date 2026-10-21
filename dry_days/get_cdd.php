
<?php
header('Access-Control-Allow-Origin: *');
header('Content-Type: application/json; charset=utf-8');

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
    echo json_encode(["error" => "Počet Lat a Lon nesouhlasí"]);
    exit;
}

// Nativní parametry mřížky CDD projekce (720x360, krok 0.5°)
$ncols = 720; $nrows = 360;
$xdim = 0.5; $ydim = 0.5;
$ulx = -179.75; $uly = 89.75;
$tile_size = 45; $dlazdic_x = 16;
$celkem_kroku = 353; // 89 let x 4 sezóny

$rejstrik_cdd = json_decode(file_get_contents("registrace_cdd.json"), true);
$mapa_dlazdic = $rejstrik_cdd['dlaždice'];

$cache = [];
$vysledky = [];

$f_cdd = fopen("celysvet_cdd.gz.bin", "rb");

for ($i = 0; $i < $pocet_bodu; $i++) {
    $lat = (float)$lats[$i];
    $lon = (float)$lons[$i];
    $bod_klic = "{$lats[$i]}_{$lons[$i]}";

    // Výpočet pozice v mřížce
    $col = (int)round(($lon - $ulx) / $xdim);
    $row = (int)round(($uly - $lat) / $ydim);
    $col = max(0, min($ncols - 1, $col));
    $row = max(0, min($nrows - 1, $row));

    $tx = (int)($col / $tile_size);
    $ty = (int)($row / $tile_size);
    $tile_index = ($ty * $dlazdic_x) + $tx;

    if (!isset($mapa_dlazdic[$tile_index]) || $mapa_dlazdic[$tile_index]['delka'] == 0) {
        $vysledky[$bod_klic] = ["error" => "Oceán"];
        continue;
    }

    // Načtení a dekomprese dlaždice do paměti RAM (pouze jednou)
    if (!isset($cache[$tile_index])) {
        fseek($f_cdd, $mapa_dlazdic[$tile_index]['offset']);
        $cache[$tile_index] = gzdecode(fread($f_cdd, $mapa_dlazdic[$tile_index]['delka']));
    }

    $raw_tile = $cache[$tile_index];
    $local_col = $col % $tile_size;
    $local_row = $row % $tile_size;
    
    // CDD ukládá 1 rok = 1 bajt, délka bloku jednoho pixelu je přesně $celkem_kroku bajtů
    $pixel_bytes_offset = ($local_row * $tile_size + $local_col) * $celkem_kroku;
    $surove_sezony_string = substr($raw_tile, $pixel_bytes_offset, $celkem_kroku);
    
    $historie_cdd = [];
    for ($j = 0; $j < $celkem_kroku; $j++) {
        $hodnota = ord($surove_sezony_string[$j]);
        $historie_cdd[] = ($hodnota === 255) ? null : $hodnota;
    }

    $vysledky[$bod_klic] = [
        "smu_id" => $tile_index,
        "data" => $historie_cdd
    ];
}

fclose($f_cdd);
echo json_encode($vysledky);
exit;
