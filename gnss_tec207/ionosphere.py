"""Thin-shell pierce geometry and the dual-frequency TEC pipeline.

Receiver-to-satellite rays intersect a spherical shell of radius 6821 km.
The station radial defines zenith; elevations below 10 deg are excluded.
Vertical TEC = slant TEC * cos(z') with z' the angle between the ray and
the shell normal at the pierce point (single-layer mapping).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .arcs import IonoSample, calibrate
from .geodesy import C_LIGHT, OMEGA_E, ecef_to_geodetic, epoch_to_seconds
from .interp import satellite_position
from .rinex import RinexData
from .sp3 import Sp3Data

SHELL_R = 6821e3
MIN_ELEVATION_DEG = 10.0


@dataclass
class Site:
    lat_deg: float
    lon_deg: float
    height_m: float


@dataclass
class IonoRecord:
    time: str
    prn: str
    arc_id: str | None
    status: str
    reason: str | None
    slant_tec_tceu: float | None = None
    vertical_tec_tceu: float | None = None
    pierce_lat_deg: float | None = None
    pierce_lon_deg: float | None = None
    elevation_deg: float | None = None


def geodetic_to_ecef(lat_deg: float, lon_deg: float, h: float) -> np.ndarray:
    from .geodesy import WGS84_A, WGS84_E2
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    return np.array([
        (n + h) * math.cos(lat) * math.cos(lon),
        (n + h) * math.cos(lat) * math.sin(lon),
        (n * (1.0 - WGS84_E2) + h) * math.sin(lat),
    ])


def pierce_point(rx: np.ndarray, sat: np.ndarray) -> tuple[np.ndarray, float, float] | None:
    """Intersect the receiver->sat ray with the SHELL_R sphere.

    Returns (pierce_ecef, elevation_deg, cos_zenith) or None when below
    10 deg or when no forward intersection exists.
    """
    rx_n = np.linalg.norm(rx)
    d = sat - rx
    d_unit = d / np.linalg.norm(d)
    rx_u = rx / rx_n
    sin_elev = float(rx_u @ d_unit)
    elevation = math.degrees(math.asin(max(-1.0, min(1.0, sin_elev))))
    if elevation < MIN_ELEVATION_DEG:
        return None
    # |rx + k*d|^2 = SHELL_R^2, k > 0
    b = 2.0 * float(rx @ d)
    c = rx_n * rx_n - SHELL_R * SHELL_R
    disc = b * b - 4.0 * c
    if disc < 0:
        return None
    roots = [(-b + math.sqrt(disc)) / 2.0, (-b - math.sqrt(disc)) / 2.0]
    forward = [k for k in roots if k > 0.0]
    if not forward:
        return None
    pierce = rx + min(forward) * d
    # angle between ray direction and shell outward normal at pierce
    cos_z = float((pierce / SHELL_R) @ d_unit)
    return pierce, elevation, cos_z


def compute_ionosphere(rinex: RinexData, sp3: Sp3Data, site: Site,
                       biases_ns: dict[str, float]) -> list[IonoRecord]:
    rx = geodetic_to_ecef(site.lat_deg, site.lon_deg, site.height_m)
    samples_by_prn: dict[str, list[IonoSample]] = {}
    # geometry resolution keyed the same way arc outputs are
    geometry: dict[tuple[int, str], dict] = {}

    # validate bias table: finite floats, GPS prn keys
    norm_bias: dict[str, float] = {}
    for key, val in biases_ns.items():
        if not isinstance(key, str) or len(key) != 3 or key[0] != "G" \
                or not key[1:].isdigit() or not 1 <= int(key[1:]) <= 32:
            raise ValueError(f"invalid bias satellite id {key!r}; use G01..G32")
        xf = float(val)
        if xf != xf or xf in (float("inf"), float("-inf")):
            raise ValueError(f"bias for {key} is not finite")
        norm_bias[key] = xf

    for ei, epoch in enumerate(rinex.epochs):
        t_rx = epoch_to_seconds(epoch.time)
        for prn in sorted(epoch.signals):
            sig = epoch.signals[prn]
            bias = norm_bias.get(prn)
            samples_by_prn.setdefault(prn, []).append(
                IonoSample(ei, epoch.time.isoformat(), t_rx, prn,
                           dict(sig.values), dict(sig.lli), bias))
            if bias is None:
                continue
            rec = sp3.sats.get(prn)
            if rec is None:
                geometry[(ei, prn)] = {"status": "invalid",
                                       "reason": "satellite not present in SP3 file"}
                continue
            pos, reason = satellite_position(rec, t_rx)
            if pos is None:
                geometry[(ei, prn)] = {"status": "invalid", "reason": reason}
                continue
            # Sagnac-corrected satellite ECEF for geometric direction at rx time
            tau = float(np.linalg.norm(pos - rx) / C_LIGHT)
            theta = OMEGA_E * tau
            rot = np.array([[math.cos(theta), math.sin(theta), 0.0],
                            [-math.sin(theta), math.cos(theta), 0.0],
                            [0.0, 0.0, 1.0]])
            sat_ecef = rot @ pos
            hit = pierce_point(rx, sat_ecef)
            if hit is None:
                geometry[(ei, prn)] = {
                    "status": "excluded",
                    "reason": f"elevation below {MIN_ELEVATION_DEG:.0f} deg or no shell intersection"}
                continue
            pierce, elev, cos_z = hit
            plat, plon, _ = ecef_to_geodetic(pierce)
            geometry[(ei, prn)] = {"status": "ok", "pierce_lat": plat,
                                   "pierce_lon": plon, "cos_z": cos_z,
                                   "elevation": elev}

    arc_out = calibrate(samples_by_prn)

    records: list[IonoRecord] = []
    for ei, epoch in enumerate(rinex.epochs):
        for prn in sorted(epoch.signals):
            a = arc_out[(ei, prn)]
            g = geometry.get((ei, prn))
            rec = IonoRecord(time=epoch.time.isoformat(), prn=prn,
                             arc_id=a.arc_id, status=a.status, reason=a.reason,
                             slant_tec_tceu=a.slant_tec)
            if a.status == "ok" and g is not None and g["status"] == "ok":
                rec.vertical_tec_tceu = a.slant_tec * g["cos_z"]
                rec.pierce_lat_deg = g["pierce_lat"]
                rec.pierce_lon_deg = g["pierce_lon"]
                rec.elevation_deg = g["elevation"]
            elif a.status == "ok":
                # valid arc but geometry unusable at this epoch: must not vanish
                rec.status = g["status"] if g else "invalid"
                rec.reason = g["reason"] if g else "missing ephemeris/geometry"
                rec.arc_id = a.arc_id
                rec.slant_tec_tceu = a.slant_tec
            elif g is not None and g["status"] != "ok" and rec.reason is None:
                rec.reason = g["reason"]
            records.append(rec)
    return records
