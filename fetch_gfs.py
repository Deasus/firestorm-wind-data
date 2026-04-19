"""
Fetch latest GFS wind data from NOAA NOMADS and convert to JSON.
Outputs JSON in the same format as nullschool/cambecc earth project.
Fetches surface (10m) and 850mb (smoke transport) wind data.
"""

import json
import struct
import os
import sys
from datetime import datetime, timedelta, timezone
import requests
import numpy as np

DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

# NOAA NOMADS GFS filter endpoint
NOMADS_BASE = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_1p00.pl"

def get_latest_gfs_cycle():
    """Determine the most recent GFS cycle that should be available.
    GFS runs at 00z, 06z, 12z, 18z. Data becomes available ~3.5 hours after.
    """
    now = datetime.now(timezone.utc)
    # Try cycles from most recent backwards
    for hours_ago in range(0, 24, 6):
        cycle_time = now - timedelta(hours=hours_ago + 4)  # 4hr processing delay
        cycle_hour = (cycle_time.hour // 6) * 6
        cycle_dt = cycle_time.replace(hour=cycle_hour, minute=0, second=0, microsecond=0)
        yield cycle_dt

def download_gfs_grib(cycle_dt, level_param, var_params, output_path):
    """Download a filtered GFS GRIB2 file from NOMADS."""
    date_str = cycle_dt.strftime("%Y%m%d")
    hour_str = cycle_dt.strftime("%H")
    
    params = {
        "file": f"gfs.t{hour_str}z.pgrb2.1p00.f000",
        level_param: "on",
        "dir": f"/gfs.{date_str}/{hour_str}/atmos"
    }
    for var in var_params:
        params[var] = "on"
    
    url = NOMADS_BASE
    print(f"  Downloading: cycle={date_str}/{hour_str}z level={level_param}")
    
    try:
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code == 200 and len(resp.content) > 100:
            with open(output_path, "wb") as f:
                f.write(resp.content)
            print(f"  Downloaded {len(resp.content)} bytes")
            return True
        else:
            print(f"  HTTP {resp.status_code}, size={len(resp.content)}")
            return False
    except Exception as e:
        print(f"  Error: {e}")
        return False

def parse_grib2_simple(filepath):
    """
    Simple GRIB2 parser that extracts the data values from a GFS 1-degree file.
    This handles the specific case of GFS 1x1 degree global grid.
    Returns list of dicts with header info and data arrays.
    """
    records = []
    with open(filepath, "rb") as f:
        data = f.read()
    
    pos = 0
    while pos < len(data) - 4:
        # Find GRIB marker
        if data[pos:pos+4] != b"GRIB":
            pos += 1
            continue
        
        # Section 0: Indicator
        edition = data[pos+7]
        if edition != 2:
            pos += 1
            continue
        
        msg_len = struct.unpack(">Q", data[pos+8:pos+16])[0]
        msg_data = data[pos:pos+msg_len]
        
        try:
            record = parse_grib2_message(msg_data)
            if record:
                records.append(record)
        except Exception as e:
            print(f"  Warning: failed to parse GRIB2 message at offset {pos}: {e}")
        
        pos += msg_len
    
    return records

def parse_grib2_message(msg):
    """Parse a single GRIB2 message. Extracts grid info and data values."""
    # We need sections 3 (grid), 4 (product), 5 (data rep), 6 (bitmap), 7 (data)
    pos = 16  # Skip section 0
    
    nx = ny = 0
    lo1 = la1 = lo2 = la2 = dx = dy = 0
    param_cat = param_num = level_type = level_value = 0
    ref_time = ""
    values = None
    
    while pos < len(msg) - 4:
        if msg[pos:pos+4] == b"7777":
            break
        
        sec_len = struct.unpack(">I", msg[pos:pos+4])[0]
        sec_num = msg[pos+4]
        
        if sec_num == 1:
            # Identification section
            year = struct.unpack(">H", msg[pos+12:pos+14])[0]
            month, day, hour, minute, second = msg[pos+14], msg[pos+15], msg[pos+16], msg[pos+17], msg[pos+18]
            ref_time = f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}.000Z"
        
        elif sec_num == 3:
            # Grid definition
            nx = struct.unpack(">I", msg[pos+30:pos+34])[0]
            ny = struct.unpack(">I", msg[pos+34:pos+38])[0]
            la1 = struct.unpack(">I", msg[pos+46:pos+50])[0] / 1e6
            lo1 = struct.unpack(">I", msg[pos+50:pos+54])[0] / 1e6
            la2 = struct.unpack(">I", msg[pos+55:pos+59])[0] / 1e6
            lo2 = struct.unpack(">I", msg[pos+59:pos+63])[0] / 1e6
            dx = struct.unpack(">I", msg[pos+63:pos+67])[0] / 1e6
            dy = struct.unpack(">I", msg[pos+67:pos+71])[0] / 1e6
        
        elif sec_num == 4:
            # Product definition
            param_cat = msg[pos+9]
            param_num = msg[pos+10]
            level_type = msg[pos+22]
            level_value = struct.unpack(">I", msg[pos+23:pos+27])[0]
        
        elif sec_num == 5:
            # Data representation
            num_points = struct.unpack(">I", msg[pos+5:pos+9])[0]
            template = struct.unpack(">H", msg[pos+9:pos+11])[0]
            
            if template == 0:  # Simple packing
                ref_val = struct.unpack(">f", msg[pos+11:pos+15])[0]
                bin_scale = struct.unpack(">h", msg[pos+15:pos+17])[0]
                dec_scale = struct.unpack(">h", msg[pos+17:pos+19])[0]
                nbits = msg[pos+19]
                
                # Store for section 7
                packing = {"ref": ref_val, "bscale": bin_scale, "dscale": dec_scale, "nbits": nbits, "npoints": num_points}
        
        elif sec_num == 7:
            # Data section
            if 'packing' in dir() and packing and packing["nbits"] > 0:
                raw = msg[pos+5:pos+sec_len]
                nbits = packing["nbits"]
                ref = packing["ref"]
                bs = 2.0 ** packing["bscale"]
                ds = 10.0 ** (-packing["dscale"])
                
                # Unpack bit-packed integers
                values = []
                bit_pos = 0
                for _ in range(packing["npoints"]):
                    byte_idx = bit_pos // 8
                    bit_offset = bit_pos % 8
                    
                    # Read enough bytes
                    raw_val = 0
                    bits_remaining = nbits
                    while bits_remaining > 0:
                        if byte_idx >= len(raw):
                            raw_val = 0
                            break
                        available = 8 - bit_offset
                        take = min(available, bits_remaining)
                        mask = ((1 << take) - 1) << (available - take)
                        raw_val = (raw_val << take) | ((raw[byte_idx] & mask) >> (available - take))
                        bits_remaining -= take
                        bit_offset = 0
                        byte_idx += 1
                    
                    values.append((ref + raw_val * bs) * ds)
                    bit_pos += nbits
        
        pos += sec_len
    
    if values and nx > 0:
        # Determine parameter name
        param_name = "unknown"
        if param_cat == 2 and param_num == 2:
            param_name = "U-component_of_wind"
        elif param_cat == 2 and param_num == 3:
            param_name = "V-component_of_wind"
        
        return {
            "header": {
                "discipline": 0,
                "parameterCategory": param_cat,
                "parameterNumber": param_num,
                "parameterNumberName": param_name,
                "parameterUnit": "m.s-1",
                "refTime": ref_time,
                "surface1Type": level_type,
                "surface1Value": level_value,
                "nx": nx,
                "ny": ny,
                "lo1": lo1,
                "la1": la1,
                "lo2": lo2,
                "la2": la2,
                "dx": dx,
                "dy": dy
            },
            "data": [round(v, 2) for v in values]
        }
    
    return None

def fetch_level(level_param, var_params, output_name, level_label):
    """Fetch and convert one level of wind data."""
    print(f"\n{'='*50}")
    print(f"Fetching {level_label} wind...")
    
    grib_path = os.path.join(DATA_DIR, f"{output_name}.grib2")
    json_path = os.path.join(DATA_DIR, f"{output_name}.json")
    
    for cycle_dt in get_latest_gfs_cycle():
        if download_gfs_grib(cycle_dt, level_param, var_params, grib_path):
            # Parse GRIB2
            print(f"  Parsing GRIB2...")
            records = parse_grib2_simple(grib_path)
            
            if len(records) >= 2:
                # Find U and V components
                u_rec = next((r for r in records if r["header"]["parameterNumber"] == 2), None)
                v_rec = next((r for r in records if r["header"]["parameterNumber"] == 3), None)
                
                if u_rec and v_rec:
                    output = [u_rec, v_rec]
                    with open(json_path, "w") as f:
                        json.dump(output, f, separators=(",", ":"))
                    
                    size_kb = os.path.getsize(json_path) / 1024
                    print(f"  ✓ Saved {json_path} ({size_kb:.0f} KB)")
                    print(f"    Grid: {u_rec['header']['nx']}×{u_rec['header']['ny']}")
                    print(f"    Ref time: {u_rec['header']['refTime']}")
                    
                    # Clean up GRIB
                    os.remove(grib_path)
                    return True
                else:
                    print(f"  Could not find U/V components in {len(records)} records")
            else:
                print(f"  Only {len(records)} records parsed (need 2)")
        
        print(f"  Trying older cycle...")
    
    print(f"  ✗ Failed to fetch {level_label} wind data")
    return False

def write_metadata(results):
    """Write a metadata file with timestamps."""
    meta = {
        "updated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "levels": results
    }
    meta_path = os.path.join(DATA_DIR, "meta.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"\n✓ Metadata written to {meta_path}")

def main():
    print("=" * 60)
    print("FIRESTORM Wind Data Updater")
    print(f"Time: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    print("=" * 60)
    
    results = {}
    
    # Surface wind (10m above ground)
    if fetch_level("lev_10_m_above_ground", ["var_UGRD", "var_VGRD"],
                    "current-wind-surface", "Surface (10m)"):
        results["surface"] = "ok"
    else:
        results["surface"] = "failed"
    
    # 850mb wind (for smoke transport)
    if fetch_level("lev_850_mb", ["var_UGRD", "var_VGRD"],
                    "current-wind-850mb", "850mb (smoke transport)"):
        results["850mb"] = "ok"
    else:
        results["850mb"] = "failed"
    
    write_metadata(results)
    
    print("\n" + "=" * 60)
    print("Done!")
    for level, status in results.items():
        icon = "✓" if status == "ok" else "✗"
        print(f"  {icon} {level}: {status}")

if __name__ == "__main__":
    main()
