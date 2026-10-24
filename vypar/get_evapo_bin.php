<?php
header('Access-Control-Allow-Origin: *');
header('Content-Type: application/octet-stream');

// Očekáváme indexy matice z klientského ořezu
if (!isset($_GET['r_min']) || !isset($_GET['r_max']) || !isset($_GET['c_min']) || !isset($_GET['c_max'])) {
    http_response_code(400);
    exit;
}

$r_min = max(0, min(127, (int)$_GET['r_min']));
$r_max = max(0, min(127, (int)$_GET['r_max']));
$c_min = max(0, min(255, (int)$_GET['c_min']));
$c_max = max(0, min(255, (int)$_GET['c_max']));

$ncols = 256;
$celkem_mesicu = 1020;
$soubor_bin = "celysvet_evapo.raw.bin";

if (!file_exists($soubor_bin)) {
    http_response_code(500);
    exit;
}

$f = fopen($soubor_bin, "rb");

// Sekvenčně vyřezáváme řádky matice v požadovaném Bounding-Boxu
for ($r = $r_min; $r <= $r_max; $r++) {
    // Spočítáme počáteční offset pro daný řádek a rozsah sloupců
    $start_pixel_idx = ($r * $ncols) + $c_min;
    $offset_bajtů = $start_pixel_idx * $celkem_mesicu;
    
    // Počet pixelů v řádku, které chceme vzít
    $pocet_pixelu_v_rade = ($c_max - $c_min) + 1;
    $delka_cteni_bajtu = $pocet_pixelu_v_rade * $celkem_mesicu;

    fseek($f, $offset_bajtů);
    echo fread($f, $delka_cteni_bajtu);
}

fclose($f);
exit;
