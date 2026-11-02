<?php
// Hromadné načtení dlaždic: vstup JSON {"tile_ids":[...]}, výstup binární stream:
// pro každou dlaždici v pořadí dotazu: uint32 BE délka + surový gzip blok (délka 0 = dlaždice bez dat)
header('Access-Control-Allow-Origin: *');
header('Access-Control-Allow-Headers: Content-Type');

if ($_SERVER['REQUEST_METHOD'] === 'OPTIONS') { http_response_code(204); exit; }
if ($_SERVER['REQUEST_METHOD'] !== 'POST') { http_response_code(405); exit; }

function chyba($kod, $zprava) {
    http_response_code($kod);
    header('Content-Type: text/plain; charset=utf-8');
    echo $zprava;
    exit;
}

$data_post = json_decode(file_get_contents('php://input'), true);
if (!isset($data_post['tile_ids']) || !is_array($data_post['tile_ids'])) chyba(400, 'Chybi tile_ids');
$tile_ids = $data_post['tile_ids'];
if (count($tile_ids) > 64) chyba(400, 'Prilis mnoho dlazdic');

$soubor_bin  = __DIR__ . "/celysvet_srazky.gz.bin";
$soubor_json = __DIR__ . "/registrace_srazky.json";
if (!is_file($soubor_bin))  chyba(500, 'Chybi celysvet_srazky.gz.bin');
if (!is_file($soubor_json)) chyba(500, 'Chybi registrace_srazky.json');

$rejstrik = json_decode(file_get_contents($soubor_json), true);
if (!is_array($rejstrik)) chyba(500, 'Rejstrik neni platny JSON');
$mapa = $rejstrik['dlazdice'] ?? $rejstrik['dlaždice'] ?? null;   // funguje s ASCII i původním klíčem
if (!is_array($mapa)) chyba(500, 'V rejstriku chybi seznam dlazdic');

$velikost_bin = filesize($soubor_bin);

// Nejdřív sestavit plán (kontrola mezí) a spočítat Content-Length
$plan = [];
$celkem = 0;
foreach ($tile_ids as $id) {
    $i = (int)$id;
    $delka = 0; $offset = 0;
    if (isset($mapa[$i])) {
        $offset = (int)$mapa[$i]['offset'];
        $delka  = (int)$mapa[$i]['delka'];
        if ($offset < 0 || $delka < 0 || $offset + $delka > $velikost_bin) { $delka = 0; $offset = 0; }
    }
    $plan[] = [$offset, $delka];
    $celkem += 4 + $delka;
}

// Čistý binární výstup: žádná kompresní vrstva PHP ani dřívější buffer, který by stream poškodil
@ini_set('zlib.output_compression', '0');
while (ob_get_level() > 0) ob_end_clean();
header('Content-Type: application/octet-stream');
header('Content-Length: ' . $celkem);
header('Cache-Control: no-store');

$f = fopen($soubor_bin, 'rb');
foreach ($plan as [$offset, $delka]) {
    echo pack('N', $delka);
    if ($delka > 0) {
        fseek($f, $offset);
        $zbyva = $delka;
        while ($zbyva > 0) {
            $kus = fread($f, min(262144, $zbyva));
            if ($kus === false || $kus === '') break;
            echo $kus;
            $zbyva -= strlen($kus);
        }
    }
}
fclose($f);
