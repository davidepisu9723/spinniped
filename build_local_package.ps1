# Build and test Spinniped without creating generated folders in this Git repo.
#
# The external folder will contain:
#   source  - a copy of the package source
#   tests   - a copy of the test suite
#   venv    - the isolated Python environment
#   dist    - the wheel and source archive

$ErrorActionPreference = "Stop"

# Repository containing this script.
$Repo = $PSScriptRoot

# Dedicated workspace beside the repository, not inside it.
$Work = Join-Path (Split-Path $Repo -Parent) "spinniped-package-build"
$Source = Join-Path $Work "source"
$Tests = Join-Path $Work "tests"
$Venv = Join-Path $Work "venv"
$Dist = Join-Path $Work "dist"

Write-Host "Preparing external workspace: $Work"

# Start with an empty external workspace.
if (Test-Path $Work) {
    Remove-Item $Work -Recurse -Force
}

New-Item $Source -ItemType Directory | Out-Null
New-Item $Tests -ItemType Directory | Out-Null
New-Item $Dist -ItemType Directory | Out-Null

# Copy the files required to build the package.
Copy-Item (Join-Path $Repo "pyproject.toml") $Source
Copy-Item (Join-Path $Repo "README.md") $Source
Copy-Item (Join-Path $Repo "LICENSE") $Source
Copy-Item (Join-Path $Repo "spinniped") $Source -Recurse

# Copy tests separately so they use the installed wheel.
Copy-Item (Join-Path $Repo "tests\*") $Tests -Recurse

Write-Host "Creating virtual environment..."
python -m venv $Venv
$Python = Join-Path $Venv "Scripts\python.exe"

Write-Host "Installing build tools..."
& $Python -m pip install --upgrade pip build

Write-Host "Building package..."
& $Python -m build --outdir $Dist $Source

# Select the wheel produced by the build.
$Wheel = Get-ChildItem $Dist -Filter "*.whl" | Select-Object -First 1
if (-not $Wheel) {
    throw "No wheel was created."
}

Write-Host "Installing wheel and test dependencies..."
& $Python -m pip install "$($Wheel.FullName)[test]"

Write-Host "Checking installed package..."
Push-Location $Work
try {
    & $Python -c "import spinniped; print('Version:', spinniped.__version__); print('Installed from:', spinniped.__file__)"

    Write-Host "Running tests..."
    & $Python -m pytest $Tests -q
}
finally {
    Pop-Location
}

Write-Host ""
Write-Host "Done."
Write-Host "Wheel: $($Wheel.FullName)"
