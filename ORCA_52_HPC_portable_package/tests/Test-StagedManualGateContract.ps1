[CmdletBinding()]
param(
    [switch]$RuntimeSmoke
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$Runner = Join-Path $ProjectRoot 'run_serial.ps1'
$JobsRoot = Join-Path $ProjectRoot 'jobs'

function Assert-Contract([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw "CONTRACT TEST FAILED: $Message" }
}

function Get-StageIndex([object[]]$Rows, [string]$Stage) {
    for ($i = 0; $i -lt $Rows.Count; $i++) {
        if ($Rows[$i].stage -eq $Stage) { return $i }
    }
    return -1
}

$tokens = $null
$parseErrors = $null
[System.Management.Automation.Language.Parser]::ParseFile(
    $Runner,
    [ref]$tokens,
    [ref]$parseErrors
) | Out-Null
Assert-Contract ($parseErrors.Count -eq 0) 'run_serial.ps1 has PowerShell parse errors.'

$phase1 = @(Import-Csv -LiteralPath (Join-Path $JobsRoot 'stage_manifest_phase1_goat_template.csv') -Encoding UTF8)
$full = @(Import-Csv -LiteralPath (Join-Path $JobsRoot 'stage_manifest_template.csv') -Encoding UTF8)
Assert-Contract ($phase1.Count -gt 0) 'The phase-1 manifest is empty.'
Assert-Contract ($phase1[-1].stage -eq 'goat_ensemble_selection') 'Phase 1 must end at goat_ensemble_selection.'
Assert-Contract (@($phase1 | Where-Object { $_.stage -match '^goat_start\d+$' }).Count -ge 10) 'Phase 1 has fewer than 10 GOAT starts.'
Assert-Contract (@($phase1 | Where-Object { $_.stage -match '^(opt_freq|shielding|j_coupling)_' }).Count -eq 0) 'Phase 1 contains a downstream ORCA job.'

$goatGateIndex = Get-StageIndex -Rows $full -Stage 'goat_ensemble_selection'
$optIndex = Get-StageIndex -Rows $full -Stage 'opt_freq_conf01'
$minimumGateIndex = Get-StageIndex -Rows $full -Stage 'minimum_boltzmann_audit'
$shieldIndex = Get-StageIndex -Rows $full -Stage 'shielding_conf01'
$jIndex = Get-StageIndex -Rows $full -Stage 'j_coupling_conf01'
Assert-Contract ($goatGateIndex -ge 0 -and $goatGateIndex -lt $optIndex) 'goat_ensemble_selection is not before opt/freq.'
Assert-Contract ($minimumGateIndex -gt $optIndex -and $minimumGateIndex -lt $shieldIndex -and $minimumGateIndex -lt $jIndex) 'minimum_boltzmann_audit is not between opt/freq and NMR properties.'
$pendingRows = @($full | Where-Object { $_.stage -match '^(opt_freq|shielding|j_coupling)_' })
Assert-Contract ($pendingRows.Count -gt 0 -and @($pendingRows | Where-Object { $_.readiness -ne 'PENDING_UPSTREAM' }).Count -eq 0) 'Downstream template rows are not all PENDING_UPSTREAM.'

$goatGate = Get-Content -LiteralPath (Join-Path $JobsRoot 'goat_ensemble_selection_gate_template.json') -Raw -Encoding UTF8 | ConvertFrom-Json
$minimumGate = Get-Content -LiteralPath (Join-Path $JobsRoot 'minimum_boltzmann_gate_template.json') -Raw -Encoding UTF8 | ConvertFrom-Json
Assert-Contract ($goatGate.stage -eq 'goat_ensemble_selection') 'GOAT gate template has the wrong stage identity.'
Assert-Contract ($goatGate.approval_status -eq 'DRAFT_NOT_APPROVED') 'GOAT gate template must be non-approved by default.'
Assert-Contract (@($goatGate.selected_conformers).Count -gt 0) 'GOAT gate template has no conformer-selection record.'
Assert-Contract ($minimumGate.stage -eq 'minimum_boltzmann_audit') 'Minimum/Boltzmann gate template was removed or renamed.'
Assert-Contract ($minimumGate.approval_status -eq 'DRAFT_NOT_APPROVED') 'Minimum/Boltzmann template must be non-approved by default.'

$stageMethods = Get-Content -LiteralPath (Join-Path $ProjectRoot 'config\stage_methods.json') -Raw -Encoding UTF8 | ConvertFrom-Json
$methodStages = @($stageMethods.stages | ForEach-Object { $_.stage })
$methodGoatGate = [array]::IndexOf($methodStages, 'goat_ensemble_selection')
$methodOpt = [array]::IndexOf($methodStages, 'opt_freq')
$methodMinimumGate = [array]::IndexOf($methodStages, 'minimum_boltzmann_audit')
$methodShielding = [array]::IndexOf($methodStages, 'shielding')
Assert-Contract ($methodGoatGate -ge 0 -and $methodGoatGate -lt $methodOpt) 'stage_methods.json does not gate GOAT before opt/freq.'
Assert-Contract ($methodMinimumGate -gt $methodOpt -and $methodMinimumGate -lt $methodShielding) 'stage_methods.json does not retain minimum_boltzmann_audit after opt/freq.'

$runnerText = Get-Content -LiteralPath $Runner -Raw -Encoding UTF8
foreach ($requiredText in @(
    "function Assert-GoatEnsembleSelectionGate",
    "`$expectedGate = 'goat_ensemble_selection'",
    "`$expectedGate = 'minimum_boltzmann_audit'",
    "selected_coordinate_sha256",
    "optimized_coordinate_sha256"
)) {
    Assert-Contract ($runnerText.Contains($requiredText)) "Runner is missing required gate contract text: $requiredText"
}
$mainLoopIndex = $runnerText.IndexOf("foreach (`$stageRow in `$stages)")
$pendingRefusalIndex = $runnerText.IndexOf('is PENDING_UPSTREAM', $mainLoopIndex)
$invokeIndex = $runnerText.IndexOf("Invoke-OrcaStage -Row `$row", $mainLoopIndex)
Assert-Contract ($mainLoopIndex -ge 0 -and $pendingRefusalIndex -gt $mainLoopIndex -and $invokeIndex -gt $pendingRefusalIndex) 'PENDING_UPSTREAM refusal is not before Invoke-OrcaStage.'

if ($RuntimeSmoke) {
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $Runner -Status | Out-Host
    if ($LASTEXITCODE -ne 0) { throw 'Status smoke test failed.' }
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $Runner -Compound F3_019 | Out-Host
    if ($LASTEXITCODE -ne 0) { throw 'Dry-run smoke test failed.' }
}

Write-Host 'PASS: staged GOAT/manual-gate/opt-freq/minimum-Boltzmann contract is intact.'
Write-Host 'PASS: no ORCA execution path was invoked by this test.'
