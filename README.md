# FIRESTORM Wind Data

Automated GFS wind data pipeline for the FIRESTORM wildfire intelligence platform.

## What this does

Every 3 hours, a GitHub Action:

1. Downloads the latest GFS forecast from NOAA NOMADS
2. Extracts surface (10m) and 850mb wind U/V components  
3. Converts GRIB2 to JSON (nullschool-compatible format)
4. Commits the updated JSON files to this repo

## Data files

| File | Description | Use |
|------|-------------|-----|
| `data/current-wind-surface.json` | Surface wind (10m) | Wind flow animation |
| `data/current-wind-850mb.json` | 850mb wind (~5000ft) | Smoke transport direction |
| `data/meta.json` | Update timestamps | Data freshness tracking |

## Usage in FIRESTORM

Fetch the raw JSON from GitHub:

```javascript
fetch('https://raw.githubusercontent.com/YOUR_USER/firestorm-wind-data/main/data/current-wind-surface.json')
  .then(r => r.json())
  .then(data => {
    // data[0] = U-component, data[1] = V-component
    // Each has .header (grid info) and .data (values)
  });
```

## Data source

- **GFS** (Global Forecast System) from NOAA/NCEP
- 1° resolution, global coverage (360×181 grid)
- Updated 4x daily (00z, 06z, 12z, 18z)
- [NOAA NOMADS](https://nomads.ncep.noaa.gov)

## Manual trigger

Go to Actions tab → "Update GFS Wind Data" → "Run workflow"
