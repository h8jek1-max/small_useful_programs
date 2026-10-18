<?php
header('Access-Control-Allow-Origin: *');
header('Content-Type: application/json; charset=utf-8');

if (!isset($_GET['lat']) || !isset($_GET['lon'])) {
    http_response_code(400);
    echo json_encode(["error" => "Chybí parametry lat/lon"]);
    exit;
}

$lat = (float)$_GET['lat'];
$lon = (float)$_GET['lon'];

if ($lat > 60.0 || $lat < -60.0) {
    http_response_code(404);
    echo json_encode(["error" => "Mimo rozsah"]);
    exit;
}

$ncols = 43200; $nrows = 21600;
$xdim = 0.00833333333333333; $ydim = 0.00833333333333333;
$ulx = -179.995833333333; $uly = 89.9958333333333;
$tile_size = 1080; $dlazdic_x = 40;

$col = (int)(($lon - $ulx) / $xdim);
$row = (int)(($uly - $lat) / $ydim);
$col = max(0, min($ncols - 1, $col));
$row = max(0, min($nrows - 1, $row));

$tx = (int)($col / $tile_size);
$ty = (int)($row / $tile_size);
$tile_index = ($ty * $dlazdic_x) + $tx;

$rejstrik_mapy = json_decode(file_get_contents("registrace_mapy.json"), true);
$dlaždice_mapy_info = $rejstrik_mapy['mapa'][$tile_index];

if ($dlaždice_mapy_info['delka'] == 0) {
    echo json_encode(["error" => "Oceán / Žádná data"]);
    exit;
}

$f_map = fopen("svet_mapa_precizni.gz.bin", "rb");
fseek($f_map, $dlaždice_mapy_info['offset']);
$gzip_map_data = fread($f_map, $dlaždice_mapy_info['delka']);
fclose($f_map);

$raw_map_data = gzdecode($gzip_map_data);
$local_col = $col % $tile_size;
$local_row = $row % $tile_size;
$pixel_offset = ($local_row * $tile_size + $local_col) * 2;

$smu_id_data = unpack("v", substr($raw_map_data, $pixel_offset, 2));
$smu_id = $smu_id_data[1];

if ($smu_id == 65535 || $smu_id == 0) {
    echo json_encode(["error" => "Mimo pevninu"]);
    exit;
}

// KLÍČOVÁ ZMĚNA: Odkaz na kompletní databázi a index
$index_pud = json_decode(file_get_contents("pudy_index_komplet.json"), true);

if (!isset($index_pud[$smu_id])) {
    echo json_encode(["error" => "Profil nenalezen", "smu_id" => $smu_id]);
    exit;
}

$profil_info = $index_pud[$smu_id];
$f_pudy = fopen("pudy_data_komplet.gz.bin", "rb");
fseek($f_pudy, $profil_info['offset']);
$gzip_puda_data = fread($f_pudy, $profil_info['delka']);
fclose($f_pudy);

echo gzdecode($gzip_puda_data);
