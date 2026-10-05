"""Generate a synthetic RINEX 3.04 obs file and SP3-c ephemeris for testing.

Simulates 8 GPS satellites on circular orbits, a fixed receiver with a clock
offset, dual-frequency ionosphere delays and per-satellite merged code
biases; writes examples/obs.rnx + examples/eph.sp3 + examples/biases.json.

Injected features for arc testing:
  * G03: cycle slip (phase jump ~2 m + LLI=1) after epoch 4
  * G07: C2W missing at epoch 3 (mid-arc break)
  * G08: omitted from biases.json -> every sample marked invalid
"""

import datetime as dt
import math
import os

import numpy as np

GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)
C = 299792458.0
OMEGA = 7.2921151467e-5
F1 = 1575.42e6
F2 = 1227.60e6
LAM1, LAM2 = C / F1, C / F2
M_PER_TECU = 40.3e16 * (1.0 / F2 ** 2 - 1.0 / F1 ** 2)

# receiver: lat 30N, lon 114E, h 50 m
LAT, LON, H = math.radians(30.0), math.radians(114.0), 50.0
A = 6378137.0
F = 1 / 298.257223563
E2 = F * (2 - F)
N = A / math.sqrt(1 - E2 * math.sin(LAT) ** 2)
RX = np.array([(N + H) * math.cos(LAT) * math.cos(LON),
               (N + H) * math.cos(LAT) * math.sin(LON),
               (N * (1 - E2) + H) * math.sin(LAT)])
RX_CLK = 1.5e-4  # seconds

ORB_R = 26560e3
N_SATS = 8


def _enu_rotation(lat: float, lon: float):
    return np.array([
        [-math.sin(lon), math.cos(lon), 0.0],
        [-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)],
        [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)],
    ])


def _pick_orbit_params(t0_s: float):
    """Choose (raan, u0) pairs giving high, azimuth-spread elevation at site."""
    enu = _enu_rotation(LAT, LON)
    chosen: list[tuple[float, float, float]] = []
    grid = [(math.radians(ra), math.radians(u))
            for ra in range(0, 360, 10) for u in range(0, 360, 5)]
    for s in range(N_SATS):
        target_az = math.radians(s * 360.0 / N_SATS)
        best = None
        for raan, u0 in grid:
            p = sat_pos_raw(raan, u0, t0_s)
            e = enu @ (p - RX)
            horiz = math.hypot(e[0], e[1])
            elev = math.degrees(math.atan2(e[2], horiz))
            az = math.degrees(math.atan2(e[0], e[1])) % 360.0
            # require decent elevation; score by azimuth match + separation
            if elev < 25.0:
                continue
            daz = min(abs(az - math.degrees(target_az)),
                      360 - abs(az - math.degrees(target_az)))
            sep = min(
                (math.degrees(math.acos(
                    np.clip(p @ q / (ORB_R * ORB_R), -1, 1)))
                for q in [sat_pos_raw(r, u, t0_s) for r, u, _ in chosen]),
                default=180.0)
            score = -daz + 0.2 * sep
            if best is None or score > best[0]:
                best = (score, raan, u0)
        if best is None:
            raise RuntimeError("no visible orbit candidate")
        chosen.append((best[1], best[2], 0.0))
    return [(r, u) for r, u, _ in chosen]


def sat_pos_raw(raan: float, u0: float, t: float) -> np.ndarray:
    inc = math.radians(55.0)
    n = 2 * math.pi / 43082.0
    u = u0 + n * t
    x_o = ORB_R * math.cos(u)
    y_o = ORB_R * math.sin(u)
    x1 = x_o
    y1 = y_o * math.cos(inc)
    z1 = y_o * math.sin(inc)
    x2 = x1 * math.cos(raan) - y1 * math.sin(raan)
    y2 = x1 * math.sin(raan) + y1 * math.cos(raan)
    theta = OMEGA * t
    x3 = x2 * math.cos(theta) + y2 * math.sin(theta)
    y3 = -x2 * math.sin(theta) + y2 * math.cos(theta)
    return np.array([x3, y3, z1])


# orbit parameters are chosen in main() once t0 is known; default keeps import safe
ORBIT_PARAMS: list[tuple[float, float]] | None = None


def sat_pos(prn_idx: int, t: float) -> np.ndarray:
    """Circular orbit position at t seconds since GPS epoch (inertial-ish frame
    rotated to ECEF at time t)."""
    raan, u0 = ORBIT_PARAMS[prn_idx]
    return sat_pos_raw(raan, u0, t)


def sat_clk(prn_idx: int, t: float) -> float:
    return 2e-5 * math.sin(2 * math.pi * t / 86400.0 + prn_idx) + prn_idx * 1e-7


def fmt_dt(t: dt.datetime) -> str:
    return (f"{t.year:5d}{t.month:3d}{t.day:3d}{t.hour:3d}{t.minute:3d}"
            f"{t.second + t.microsecond / 1e6:11.7f}")


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    t0 = dt.datetime(2024, 6, 1, 0, 10, 0, tzinfo=dt.timezone.utc)
    global ORBIT_PARAMS
    ORBIT_PARAMS = _pick_orbit_params((t0 - GPS_EPOCH).total_seconds())
    n_epochs = 8
    step = 30  # seconds between obs epochs

    # --- SP3: 900 s spacing, wide margin for 8-node windows ---
    sp3_start = t0 - dt.timedelta(seconds=3600)
    sp3_n = 12
    lines = []
    first = sp3_start
    lines.append(f"#cP{first.year:4d} {first.month:2d} {first.day:2d} "
                 f"{first.hour:2d} {first.minute:2d} {first.second:11.8f} "
                 f"{sp3_n:7d}    u+U  IGS14 FIT  AIUB")
    lines.append("## 0000      0.00000000    900.00000000   00000   0.0000000000000")
    lines.append("+   10   G01G02G03G04G05G06G07G08G09G10  0  0  0  0  0  0  0")
    lines.append("++          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0")
    lines.append("%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc")
    lines.append("%f  0.0000000  0.000000000  0.00000000000  0.000000000000000")
    lines.append("%i    0    0    0    0         0         0         0         0")
    lines.append("/* SYNTHETIC TEST EPHEMERIS - NOT REAL NAVIGATION DATA")
    for k in range(sp3_n):
        t = sp3_start + dt.timedelta(seconds=900 * k)
        ts = (t - GPS_EPOCH).total_seconds()
        lines.append(f"*  {t.year:4d} {t.month:2d} {t.day:2d} {t.hour:2d} "
                     f"{t.minute:2d} {t.second:11.8f}")
        for s in range(N_SATS):
            p = sat_pos(s, ts) / 1000.0
            c_us = sat_clk(s, ts) * 1e6
            lines.append(f"PG{s + 1:02d}{p[0]:14.6f}{p[1]:14.6f}{p[2]:14.6f}"
                         f"{c_us:14.6f}")
    lines.append("EOF")
    with open(os.path.join(outdir, "eph.sp3"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # --- RINEX obs ---
    rng = np.random.default_rng(42)
    # merged code bias per satellite in nanoseconds (G08 intentionally absent)
    biases_ns = {f"G{s + 1:02d}": 0.5 + 0.25 * s for s in range(N_SATS - 1)}
    lines = []
    lines.append(f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE")
    lines.append(f"{'SYNTH':<20}{'TEST':<20}{'20240601 000000 UTC':<20}PGM / RUN BY / DATE")
    lines.append(f"{RX[0]:14.4f}{RX[1]:14.4f}{RX[2]:14.4f}{'':18}APPROX POSITION XYZ")
    lines.append(f"{0:6d}{'':54}RCV CLOCK OFFS APPL")
    lines.append(f"{'G':<1}{4:5d} {'C1C C2W L1C L2W':<53}SYS / # / OBS TYPES")
    lines.append(f"{t0.year:6d}{t0.month:6d}{t0.day:6d}{t0.hour:6d}{t0.minute:6d}"
                 f"{t0.second:12.7f}GPS{'':9}TIME OF FIRST OBS")
    lines.append(f"{'':60}END OF HEADER")
    # fixed integer ambiguities per satellite
    amb1 = {s: 100000.0 + s * 137.0 for s in range(N_SATS)}
    amb2 = {s: 70000.0 + s * 91.0 for s in range(N_SATS)}
    for e in range(n_epochs):
        t = t0 + dt.timedelta(seconds=step * e)
        ts = (t - GPS_EPOCH).total_seconds()
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{t.second:11.7f}  0{N_SATS:3d}{'':6}")
        for s in range(N_SATS):
            # transmit time iteration matching the solver
            pr_guess = 22000e3
            for _ in range(4):
                t_tx = ts - (pr_guess / C + sat_clk(s, ts - pr_guess / C))
                p = sat_pos(s, t_tx)
                # Earth rotation during flight
                tau = np.linalg.norm(p - RX) / C
                th = OMEGA * tau
                rot = np.array([[math.cos(th), math.sin(th), 0],
                                [-math.sin(th), math.cos(th), 0],
                                [0, 0, 1]])
                pr_guess = (np.linalg.norm(rot @ p - RX)
                            + C * RX_CLK - C * sat_clk(s, t_tx))
            geo = pr_guess - C * RX_CLK + C * sat_clk(s, t_tx)
            # smoothly varying slant TEC (can be negative for one sat)
            tec = 12.0 + 4.0 * math.sin(2 * math.pi * e / n_epochs + s * 0.9) \
                + (-25.0 if s == N_SATS - 2 else 0.0)
            iono1 = 40.3e16 / F1 ** 2 * tec
            iono2 = 40.3e16 / F2 ** 2 * tec
            clk_m = C * RX_CLK - C * sat_clk(s, t_tx)
            p1 = geo + clk_m + iono1 + rng.normal(0, 0.2)
            p2 = geo + clk_m + iono2 + rng.normal(0, 0.2)
            # merged code bias enters the observed P2-P1 as +c*b
            if s < N_SATS - 1:
                p2 += C * biases_ns[f"G{s + 1:02d}"] * 1e-9
            # G07 (index 6): C2W missing at epoch 3 -> mid-arc break
            if s == 6 and e == 3:
                p2 = None
            # phase in cycles; phase advance = -iono
            slip = 2.0 / (LAM1 - LAM2) if (s == 2 and e >= 4) else 0.0
            l1 = (geo + clk_m - iono1) / LAM1 + amb1[s] + slip
            l2 = (geo + clk_m - iono2) / LAM2 + amb2[s]
            lli1 = 1 if (s == 2 and e == 4) else 0
            def field(v, lli=0):
                if v is None:
                    return " " * 16
                return f"{v:14.3f}{lli:1d} "
            lines.append(
                f"G{s + 1:02d}" + field(p1) + field(p2) +
                field(l1, lli1) + field(l2))
    with open(os.path.join(outdir, "obs.rnx"), "w") as f:
        f.write("\n".join(lines) + "\n")
    import json
    with open(os.path.join(outdir, "biases.json"), "w") as f:
        json.dump(biases_ns, f, indent=2)
    print("wrote examples/obs.rnx, examples/eph.sp3 and examples/biases.json")


if __name__ == "__main__":
    main()
