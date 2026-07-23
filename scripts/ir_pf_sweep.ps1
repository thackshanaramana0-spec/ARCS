# IR-prefetch N-variant sweep — pinned to P-core 0 (affinity 0x1), 5 reps each.
# Usage: powershell -ExecutionPolicy Bypass -File ir_pf_sweep.ps1
# Requires: arcs.exe built with ARCS_DNA_PROFILE support; DS1 and GIAB FASTQ at $Env:ARCS_BENCH

param(
    [string]$ArcsExe  = "C:\Users\Lenovo\OneDrive\Desktop\nu256_upgrade\arcs-clean\build\arcs.exe",
    [string]$Ds1Fq    = "C:\Temp\audit_ds1_dec.fq",
    [string]$GiabFq   = "C:\Temp\audit_giab_dec.fq",
    [string]$TmpDir   = "C:\Temp\irpf_sweep",
    [int]   $Reps     = 5,
    [int[]] $Variants = @(0, 1, 2, 4, 5)   # ARCS_IR_PF_N values (0 = baseline)
)

New-Item -ItemType Directory -Force -Path $TmpDir | Out-Null

function Run-Pinned {
    param([string]$Cmd, [hashtable]$Env)
    # Build env-prefix string for cmd /c
    $envStr = ($Env.GetEnumerator() | ForEach-Object { "set $($_.Key)=$($_.Value)" }) -join " & "
    # start /affinity 1 pins to CPU 0 (P-core 0, HT0 on i5-1235U)
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = "cmd.exe"
    $psi.Arguments = "/c ($envStr) & $Cmd"
    $psi.UseShellExecute = $false
    $psi.RedirectStandardError  = $true
    $psi.RedirectStandardOutput = $true
    $psi.CreateNoWindow = $true
    $proc = [System.Diagnostics.Process]::Start($psi)
    $proc.ProcessorAffinity = [System.IntPtr]1    # CPU 0 only
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()
    return @{ Stdout = $stdout; Stderr = $stderr; Exit = $proc.ExitCode }
}

function Extract-Profile {
    param([string]$Stderr)
    # Parse [DNA_PROFILE] lines
    $cyb = [regex]::Match($Stderr, 'cy/base\s+select=\s*([\d.]+)\s+code=\s*([\d.]+)\s+update=\s*([\d.]+)\s+other=\s*([\d.]+)\s+total=\s*([\d.]+)')
    $wall = [regex]::Match($Stderr, 'wall=([\d.]+)\s+s')
    if ($cyb.Success -and $wall.Success) {
        return [pscustomobject]@{
            Select = [double]$cyb.Groups[1].Value
            Code   = [double]$cyb.Groups[2].Value
            Update = [double]$cyb.Groups[3].Value
            Other  = [double]$cyb.Groups[4].Value
            Total  = [double]$cyb.Groups[5].Value
            Wall   = [double]$wall.Groups[1].Value
        }
    }
    return $null
}

function Median {
    param([double[]]$arr)
    $s = $arr | Sort-Object
    $n = $s.Count
    if ($n % 2 -eq 1) { return $s[($n-1)/2] }
    return ($s[$n/2-1] + $s[$n/2]) / 2.0
}

$datasets = @(
    @{ Name = "DS1";  Fq = $Ds1Fq;   Arc = "$TmpDir\ds1" },
    @{ Name = "GIAB"; Fq = $GiabFq;  Arc = "$TmpDir\giab" }
)

$results = @()

foreach ($ds in $datasets) {
    foreach ($n in $Variants) {
        Write-Host "=== $($ds.Name) ENCODE  ARCS_IR_PF_N=$n ===" -ForegroundColor Cyan
        $selArr=[double[]]@(); $codeArr=[double[]]@(); $updArr=[double[]]@()
        $totArr=[double[]]@(); $wallArr=[double[]]@()

        for ($r = 1; $r -le $Reps; $r++) {
            $arc = "$($ds.Arc)_n${n}_r${r}.arcs"
            $envVars = @{
                ARCS_DNA_PROFILE = "1"
                ARCS_ENC_NOPAR   = "1"
                ARCS_IR_PF_N     = "$n"
            }
            $cmd = "`"$ArcsExe`" compress `"$($ds.Fq)`" `"$arc`""
            $res = Run-Pinned -Cmd $cmd -Env $envVars
            $p = Extract-Profile -Stderr $res.Stderr
            if ($p) {
                $selArr  += $p.Select; $codeArr += $p.Code
                $updArr  += $p.Update; $totArr  += $p.Total; $wallArr += $p.Wall
                Write-Host "  rep $r : total=$($p.Total.ToString('F1')) cy/b  update=$($p.Update.ToString('F1'))  wall=$($p.Wall.ToString('F3'))s"
            } else {
                Write-Host "  rep $r : PARSE ERROR" -ForegroundColor Red
                Write-Host $res.Stderr
            }
        }

        if ($selArr.Count -gt 0) {
            $row = [pscustomobject]@{
                Dataset  = $ds.Name
                Mode     = "ENCODE"
                N_IR_PF  = $n
                SelMed   = Median $selArr
                CodeMed  = Median $codeArr
                UpdMed   = Median $updArr
                TotMed   = Median $totArr
                WallMed  = Median $wallArr
                Reps     = $selArr.Count
            }
            $results += $row
            Write-Host "  MEDIAN: total=$($row.TotMed.ToString('F1')) cy/b  update=$($row.UpdMed.ToString('F1'))  wall=$($row.WallMed.ToString('F3'))s" -ForegroundColor Green
        }
    }
}

# ── Print summary table ────────────────────────────────────────────────────────
Write-Host "`n========== SUMMARY (median of $Reps reps, pinned CPU0 / affinity=1) ==========" -ForegroundColor Yellow
Write-Host ("{0,-6} {1,-6} {2,-8} {3,8} {4,8} {5,8} {6,8} {7,8} {8,8}" -f `
    "DS","Mode","N_IR_PF","Sel","Code","Upd","Tot","Wall","Upd%")
Write-Host ("-" * 78)
$baseline = @{}
foreach ($r in $results) {
    $key = "$($r.Dataset)_$($r.Mode)"
    if ($r.N_IR_PF -eq 0) { $baseline[$key] = $r }
}
foreach ($r in $results) {
    $key = "$($r.Dataset)_$($r.Mode)"
    $b = $baseline[$key]
    $updPct = if ($b) { (($r.UpdMed - $b.UpdMed) / $b.UpdMed * 100) } else { 0.0 }
    $totPct = if ($b) { (($r.TotMed - $b.TotMed) / $b.TotMed * 100) } else { 0.0 }
    $tag = if ($r.N_IR_PF -eq 0) { "" } else { "($('{0:+0.1f}' -f $totPct)%)" }
    Write-Host ("{0,-6} {1,-6} {2,-8} {3,8:F1} {4,8:F1} {5,8:F1} {6,8:F1} {7,8:F3} {8,8:F1} {9}" -f `
        $r.Dataset, $r.Mode, $r.N_IR_PF,
        $r.SelMed, $r.CodeMed, $r.UpdMed, $r.TotMed, $r.WallMed,
        ($r.UpdMed / $r.TotMed * 100), $tag)
}

# Save CSV
$csvPath = "$TmpDir\results.csv"
$results | Export-Csv -Path $csvPath -NoTypeInformation
Write-Host "`nResults saved to $csvPath"
