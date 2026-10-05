"""FastAPI entry point: upload RINEX 3.04 obs + SP3-c ephemeris, get positions."""

from __future__ import annotations

from fastapi import FastAPI, File, HTTPException, UploadFile

from .errors import GnssError
from .rinex import parse_rinex
from .solver import solve_all
from .sp3 import parse_sp3

app = FastAPI(title="GNSS Pseudorange Positioning", version="0.1.0")


@app.post("/position")
async def position(rinex: UploadFile = File(...), sp3: UploadFile = File(...)):
    if rinex.filename and rinex.filename.endswith((".gz", ".zip", ".Z")):
        raise HTTPException(400, "compressed RINEX not accepted; upload uncompressed file")
    try:
        rinex_text = (await rinex.read()).decode("ascii")
    except UnicodeDecodeError:
        raise HTTPException(400, "rinex: file is not plain ASCII text")
    try:
        sp3_text = (await sp3.read()).decode("ascii")
    except UnicodeDecodeError:
        raise HTTPException(400, "sp3: file is not plain ASCII text")
    try:
        obs = parse_rinex(rinex_text)
    except GnssError as e:
        raise HTTPException(422, f"RINEX rejected: {e}")
    try:
        eph = parse_sp3(sp3_text)
    except GnssError as e:
        raise HTTPException(422, f"SP3 rejected: {e}")

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


@app.get("/health")
async def health():
    return {"status": "ok"}
