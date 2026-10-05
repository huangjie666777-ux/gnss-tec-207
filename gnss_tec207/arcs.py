"""Dual-frequency arc segmentation, phase leveling and slant TEC.

Code geometry-free combination:  P = C2W - C1C - c * bias
Phase geometry-free combination (meters): Lphi = L1C*lambda1 - L2W*lambda2
Slant TEC uses the arc-mean leveled offset so that P and Lphi agree:
    b_arc = mean(P - Lphi) over complete samples of the arc
    STEC  = (Lphi + b_arc) / K
K = 40.3e16 * (1/f2^2 - 1/f1^2)  [m per TECU].  Negative values are kept.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .geodesy import C_LIGHT

F1 = 1575.42e6  # GPS L1 Hz
F2 = 1227.60e6  # GPS L2 Hz
LAMBDA1 = C_LIGHT / F1
LAMBDA2 = C_LIGHT / F2
# meters of geometry-free delay per 1 TECU
M_PER_TECU = 40.3e16 * (1.0 / (F2 * F2) - 1.0 / (F1 * F1))

MAX_GAP_S = 120.0
JUMP_M = 1.0
MIN_ARC_SAMPLES = 3
PHASE_TYPES = ("L1C", "L2W")
CODE_TYPES = ("C1C", "C2W", "L1C", "L2W")


@dataclass
class IonoSample:
    epoch_index: int
    time_iso: str
    t: float                 # seconds since GPS epoch
    prn: str
    values: dict[str, float]  # observed quantities in m / cycles
    lli: dict[str, int]
    bias_ns: float | None     # merged per-satellite code bias, nanoseconds


@dataclass
class ArcOut:
    epoch_index: int
    prn: str
    arc_id: str | None = None
    status: str = "invalid"   # ok | short_arc | invalid
    reason: str | None = None
    slant_tec: float | None = None


def _completeness(s: IonoSample) -> str | None:
    """Return a break/invalid reason, or None when all four quantities are good."""
    missing = [name for name in CODE_TYPES if name not in s.values]
    if missing:
        return "missing observation(s): " + ",".join(missing)
    bad_lli = [name for name in PHASE_TYPES if s.lli.get(name, 0) != 0]
    if bad_lli:
        return "non-zero phase LLI: " + ",".join(f"{n}={s.lli[n]}" for n in bad_lli)
    if s.bias_ns is None:
        return "merged code bias missing for satellite; marked invalid"
    for name in CODE_TYPES:
        v = s.values[name]
        if v != v or v in (float("inf"), float("-inf")):
            return f"non-finite value in {name}"
    return None


def _combinations(s: IonoSample) -> tuple[float, float]:
    """Bias-corrected code difference (m) and phase difference in meters."""
    code = s.values["C2W"] - s.values["C1C"] - C_LIGHT * s.bias_ns * 1e-9
    phase = s.values["L1C"] * LAMBDA1 - s.values["L2W"] * LAMBDA2
    return code, phase


def calibrate(samples_by_prn: dict[str, list[IonoSample]]) -> dict[tuple[int, str], ArcOut]:
    """Segment per-satellite arcs and level phase against code.

    Arcs break on: any missing quantity, non-zero phase LLI, epoch gap > 120 s,
    or a consecutive phase-difference jump > 1 m.  Offsets are never shared
    between arcs.  Arcs need >= 3 complete samples; their members stay visible
    with status short_arc.  Every input sample produces one output record.
    """
    out: dict[tuple[int, str], ArcOut] = {}
    for prn, samples in samples_by_prn.items():
        arc_no = 0
        cur: list[IonoSample] = []
        prev_t: float | None = None
        prev_phi: float | None = None

        def flush(arc: list[IonoSample]) -> None:
            nonlocal arc_no
            if not arc:
                return
            arc_no += 1
            aid = f"{prn}-{arc_no}"
            combos = [_combinations(x) for x in arc]
            if len(arc) < MIN_ARC_SAMPLES:
                for x in arc:
                    out[(x.epoch_index, prn)] = ArcOut(
                        x.epoch_index, prn, aid, "short_arc",
                        f"arc has only {len(arc)} complete sample(s) (need >= 3)")
                return
            offset = sum(p - phi for p, phi in combos) / len(combos)
            for x, (_, phi) in zip(arc, combos):
                stec = (phi + offset) / M_PER_TECU
                out[(x.epoch_index, prn)] = ArcOut(
                    x.epoch_index, prn, aid, "ok", None, stec)

        for s in samples:
            reason = _completeness(s)
            if reason is not None:
                flush(cur)
                cur = []
                prev_t = prev_phi = None
                out[(s.epoch_index, prn)] = ArcOut(
                    s.epoch_index, prn, None, "invalid", reason)
                continue
            _, phi = _combinations(s)
            if prev_t is not None:
                if s.t - prev_t > MAX_GAP_S:
                    flush(cur)
                    cur = []
                elif abs(phi - prev_phi) > JUMP_M:
                    flush(cur)
                    cur = []
            cur.append(s)
            prev_t = s.t
            prev_phi = phi
        flush(cur)
    return out
