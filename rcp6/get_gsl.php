<?php
header('Access-Control-Allow-Origin: *');
header('Content-Type: application/json; charset=utf-8');

if (!isset($_GET['lats']) || !isset($_GET['lons'])) {
    http_response_code(400); echo json_encode(["error" => "Chybí parametry"]); exit;
}

$lats = explode(',', $_GET['lats']);
$lons = explode(',', $_GET['lons']);
$pocet_bodu = count($lats);

$ncols = 720; $nrows = 360;
$xdim = 0.5; $ydim = 0.5;
$ulx = -179.75; $uly = 89.75;
$tile_size = 45; $dlazdic_x = 16;
$celkem_let = 89; 

$rejstrik_gsl = json_decode(file_get_contents("registrace_gsl.json"), true);
$mapa_dlazdic = $rejstrik_gsl['dlaždice'];
$cache = []; $vysledky = [];

$f_gsl = fopen("celysvet_gsl.gz.bin", "rb");

for ($i = 0; $i < $pocet_bodu; $i++) {
    $lat = (float)$lats[$i]; $lon = (float)$lons[$i];
    $bod_klic = "{$lats[$i]}_{$lons[$i]}";

    $col = (int)round(($lon - $ulx) / $xdim);
    $row = (int)round(($uly - $lat) / $ydim);
    $col = max(0, min($ncols - 1, $col));
    $row = max(0, min($nrows - 1, $row));

    $tx = (int)($col / $tile_size); $ty = (int)($row / $tile_size);
    $tile_index = ($ty * $dlazdic_x) + $tx;

    if (!isset($mapa_dlazdic[$tile_index]) || $mapa_dlazdic[$tile_index]['delka'] == 0) {
        $vysledky[$bod_klic] = ["error" => "Oceán"]; continue;
    }

    if (!isset($cache[$tile_index])) {
        fseek($f_gsl, $mapa_dlazdic[$tile_index]['offset']);
        $cache[$tile_index] = gzdecode(fread($f_gsl, $mapa_dlazdic[$tile_index]['delka']));
    }

    $raw_tile = $cache[$tile_index];
    $local_col = $col % $tile_size; $local_row = $row % $tile_size;
    
    // OPRAVA: Každý rok má nyní 2 bajty (uint16), proto násobíme celkovou délku osy dvěma!
    $pixel_bytes_offset = ($local_row * $tile_size + $local_col) * ($celkem_let * 2);
    $surove_bytes = substr($raw_tile, $pixel_bytes_offset, $celkem_let * 2);
    
    // Rozbalujeme po dvojicích bajtů (formát 'v' = unsigned 16-bit little endian)
    $historie = [];
    for ($j = 0; $j < $celkem_let; $j++) {
        $val_data = unpack("v", substr($surove_bytes, $j * 2, 2));
        $val = $val_data[1];
        $historie[] = ($val === 65535) ? null : $val;
    }

    $vysledky[$bod_klic] = ["smu_id" => $tile_index, "data" => $historie];
}

fclose($f_gsl);
echo json_encode($vysledky);
