"""FastAPI entry: positioning (/position) and dual-frequency TEC (/ionosphere)."""

from __future__ import annotations

import json

import numpy as np
from fastapi import FastAPI, File, Form, HTTPException, UploadFile

from .errors import GnssError
from .ionosphere import Site, compute_ionosphere
from .rinex import parse_rinex
from .solver import solve_all
from .sp3 import parse_sp3

app = FastAPI(title="GNSS Positioning & Dual-frequency Ionosphere TEC",
              version="0.2.0")


async def _read_text(up: UploadFile, tag: str) -> str:
    if up.filename and up.filename.endswith((".gz", ".zip", ".Z")):
        raise HTTPException(400, f"{tag}: compressed files not accepted")
    try:
        return (await up.read()).decode("ascii")
    except UnicodeDecodeError:
        raise HTTPException(400, f"{tag}: file is not plain ASCII text")


async def _decode_inputs(rinex: UploadFile, sp3: UploadFile):
    rinex_text = await _read_text(rinex, "rinex")
    sp3_text = await _read_text(sp3, "sp3")
    try:
        obs = parse_rinex(rinex_text)
    except GnssError as e:
        raise HTTPException(422, f"RINEX rejected: {e}")
    try:
        eph = parse_sp3(sp3_text)
    except GnssError as e:
        raise HTTPException(422, f"SP3 rejected: {e}")
    return obs, eph


def _finite_form(name: str, raw: str) -> float:
    try:
        v = float(raw)
    except (TypeError, ValueError):
        raise HTTPException(400, f"{name} must be a finite number, got {raw!r}")
    if v != v or v in (float("inf"), float("-inf")):
        raise HTTPException(400, f"{name} must be finite")
    return v


@app.post("/position")
async def position(rinex: UploadFile = File(...), sp3: UploadFile = File(...)):
    obs, eph = await _decode_inputs(rinex, sp3)

    results = solve_all(obs.epochs, eph, obs.approx_position)
    return {
        "n_epochs": len(results),
        "n_ok": sum(1 for r in results if r.status == "ok"),
        "notes": [
            "GPS time, GPS satellites, C1C only; no atmosphere/carrier corrections.",
            "Single-epoch independent solutions; approx position used as initial guess only.",
            "Accuracy boundary: standalone C1C pseudorange positioning, "
            "meter-level at best without ionosphere/troposphere modeling.",
        ],
        "epochs": [
            {
                "time": r.time,
                "status": r.status,
                "reason": r.reason,
                "ecef_m": r.ecef_m,
                "geodetic": r.geodetic,
                "receiver_clock_s": r.receiver_clock_s,
                "used_satellites": r.used_satellites,
                "excluded_satellites": r.excluded_satellites,
                "residuals_m": r.residuals_m,
                "rms_m": r.rms_m,
            }
            for r in results
        ],
    }


@app.post("/ionosphere")
async def ionosphere(
    rinex: UploadFile = File(...),
    sp3: UploadFile = File(...),
    lat_deg: str = Form(..., description="station WGS84 latitude"),
    lon_deg: str = Form(...),
    height_m: str = Form(...),
    biases_ns: str = Form(
        ..., description='JSON object, per-satellite merged code bias ns, e.g. {"G01":2.1}'),
):
    obs, eph = await _decode_inputs(rinex, sp3)
    lat = _finite_form("lat_deg", lat_deg)
    lon = _finite_form("lon_deg", lon_deg)
    hgt = _finite_form("height_m", height_m)
    if not -90.0 <= lat <= 90.0:
        raise HTTPException(400, "lat_deg must be within [-90, 90]")
    if not -180.0 <= lon <= 180.0:
        raise HTTPException(400, "lon_deg must be within [-180, 180]")
    if not -1000.0 <= hgt <= 12000.0:
        raise HTTPException(400, "height_m outside plausible range [-1000, 12000]")
    try:
        biases = json.loads(biases_ns)
    except json.JSONDecodeError as e:
        raise HTTPException(400, f"biases_ns must be a JSON object: {e}")
    if not isinstance(biases, dict) or not biases:
        raise HTTPException(400, "biases_ns must be a non-empty JSON object")
    site = Site(lat_deg=lat, lon_deg=lon, height_m=hgt)
    try:
        records = compute_ionosphere(obs, eph, site, biases)
    except ValueError as e:
        raise HTTPException(400, str(e))
    n_ok = sum(1 for r in records if r.status == "ok")
    return {
        "n_epochs": len(obs.epochs),
        "n_records": len(records),
        "n_ok": n_ok,
        "site": {"lat_deg": lat, "lon_deg": lon, "height_m": hgt},
        "shell_radius_km": 6821,
        "min_elevation_deg": 10.0,
        "notes": [
            "code difference C2W-C1C corrected by c*bias (bias sign: TGD-style, "
            "positive delay added at P1/P2 merge convention; applied P2-P1-c*b)",
            "phase leveling per arc (no shared offsets); arcs break on missing data, "
            "non-zero LLI, gaps >120 s or phase jumps >1 m; need >=3 samples",
            "slant TEC = (L1C*lambda1 - L2W*lambda2 + arc_offset)/K, K=40.3e16*(1/f2^2-1/f1^2)",
            "vertical TEC uses a 6821 km spherical thin shell; negative TEC retained",
            "satellite positions: 8-node SP3 Lagrange interpolation at receive epoch, "
            "never across missing nodes or extrapolated",
        ],
        "records": [
            {
                "time": r.time,
                "prn": r.prn,
                "arc_id": r.arc_id,
                "status": r.status,
                "reason": r.reason,
                "slant_tec_tceu": r.slant_tec_tceu,
                "vertical_tec_tceu": r.vertical_tec_tceu,
                "pierce_lat_deg": r.pierce_lat_deg,
                "pierce_lon_deg": r.pierce_lon_deg,
                "elevation_deg": r.elevation_deg,
            }
            for r in records
        ],
    }


@app.get("/health")
async def health():
    return {"status": "ok"}
