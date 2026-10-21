<?php
// Aktivujeme zachytávání výstupu, abychom mohli eliminovat případné mezery/chyby
ob_start();

header('Access-Control-Allow-Origin: *');
header('Content-Type: application/octet-stream');

if (!isset($_GET['lat']) || !isset($_GET['lon'])) {
    http_response_code(400);
    die("Chybí parametry lat/lon");
}

$lat = (float)$_GET['lat'];
$lon = (float)$_GET['lon'];

//if ($lat > 60.0 || $lat < -60.0) {
    //http_response_code(404);
    //die("Oblast mimo obyvatelnou zónu");
//}

$X = (int)round(($lon - (-179.875)) / 0.25);
$Y = (int)round(((89.875) - $lat) / 0.25);

$X = max(0, min(1439, $X));
$Y = max(0, min(719, $Y));

$tile_size = 45;
$tx = (int)($X / $tile_size);
$ty = (int)($Y / $tile_size);
$tile_index = ($ty * 32) + $tx;

$rejstrik_raw = file_get_contents("registrace_gzip.json");
$rejstrik = json_decode($rejstrik_raw, true);
$dlaždice_info = $rejstrik['dlaždice'][$tile_index];

if ($dlaždice_info['delka'] == 0) {
    http_response_code(404);
    die("Zvolená oblast je čistý oceán");
}

$f = fopen("celysvet_sucho_gzip.bin", "rb");
fseek($f, $dlaždice_info['offset']);
$gzip_data = fread($f, $dlaždice_info['delka']);
fclose($f);

// KLÍČOVÁ OPRAVA: Vymažeme případný textový balast (např. mezery za z jiných souborů)
ob_end_clean();

// Odešleme čistá binární data
echo $gzip_data;
exit; // Okamžitě ukončíme skript, aby se nic dalšího nepřipsalo
