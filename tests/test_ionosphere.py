import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gnss_tec207.arcs import M_PER_TECU
from gnss_tec207.errors import RejectedContentError
from gnss_tec207.ionosphere import Site, compute_ionosphere
from gnss_tec207.rinex import parse_rinex
from gnss_tec207.sp3 import parse_sp3


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_synthetic.py")],
                   check=True, cwd=ROOT)


@pytest.fixture(scope="module")
def iono():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    return compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), biases)


def _by_prn(iono, prn):
    return [r for r in iono if r.prn == prn]


def test_every_observation_surfaces(iono):
    # 8 epochs * 8 satellites, nothing silently dropped
    assert len(iono) == 64


def test_missing_bias_satellite_invalid(iono):
    g08 = _by_prn(iono, "G08")
    assert len(g08) == 8
    assert all(r.status == "invalid" and "bias" in r.reason for r in g08)
    assert all(r.arc_id is None and r.slant_tec_tceu is None for r in g08)


def test_lli_cycle_slip_starts_new_arc(iono):
    g03 = _by_prn(iono, "G03")
    arcs = [r.arc_id for r in g03 if r.status == "ok"]
    assert "G03-1" in arcs and "G03-2" in arcs
    # the epoch carrying LLI=1 is itself invalid and breaks the arc
    lli_epochs = [r for r in g03 if r.status == "invalid"]
    assert lli_epochs and "LLI" in lli_epochs[0].reason


def test_missing_c2w_breaks_arc(iono):
    g07 = _by_prn(iono, "G07")
    invalid = [r for r in g07 if r.status == "invalid"]
    assert len(invalid) == 1
    assert "C2W" in invalid[0].reason
    arc_ids = {r.arc_id for r in g07 if r.status == "ok"}
    assert arc_ids == {"G07-1", "G07-2"}  # offset not shared


def test_negative_tec_retained(iono):
    g07 = [r for r in _by_prn(iono, "G07") if r.status == "ok"]
    assert all(r.slant_tec_tceu < -8.0 for r in g07)
    assert all(r.slant_tec_tceu <= r.vertical_tec_tceu < 0 for r in g07)


def test_vertical_mapping_and_pierce(iono):
    ok = [r for r in iono if r.status == "ok"]
    for r in ok:
        assert r.elevation_deg >= 10.0
        assert -90.0 <= r.pierce_lat_deg <= 90.0
        assert -180.0 <= r.pierce_lon_deg <= 180.0
        # vertical = slant * cos(z'), |cos| <= 1
        ratio = r.vertical_tec_tceu / r.slant_tec_tceu
        assert 0.0 < ratio <= 1.0 + 1e-9


def test_tec_amplitude_reconstructed(iono):
    # injected values range roughly -13..17 TECU plus noise; sanity envelope
    vals = [r.slant_tec_tceu for r in iono if r.status == "ok"]
    assert min(vals) < -10.0 and max(vals) < 20.0


def test_short_arc_after_big_gap():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    # keep only 2 epochs separated by >120 s for G01 via deleting middle epochs
    obs.epochs = [obs.epochs[0], obs.epochs[5]]  # 150 s gap > 120 s
    recs = compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), biases)
    g01 = [r for r in recs if r.prn == "G01"]
    assert all(r.status == "short_arc" and r.reason for r in g01)
    assert {r.arc_id for r in g01} == {"G01-1", "G01-2"}


def test_phase_jump_breaks_arc():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    # inject a 2 m phase jump on G01 at epoch 3 without LLI
    sig = obs.epochs[3].signals["G01"]
    sig.values["L1C"] += 2.0 / (M_PER_TECU) * 0 + 2.0 / (299792458.0 / 1575.42e6)
    recs = compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), biases)
    g01 = [r for r in recs if r.prn == "G01" and r.status == "ok"]
    assert len({r.arc_id for r in g01}) == 2


def test_ephemeris_missing_node_excludes():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    sp3_text = (ROOT / "examples" / "eph.sp3").read_text()
    lines = sp3_text.splitlines()
    # zero the coordinates of G01 on the SP3 node nearest the first rx epoch
    hit = 0
    for i, ln in enumerate(lines):
        if ln.startswith("PG01"):
            hit += 1
            if hit == 5:  # middle node within any 8-window around reception
                lines[i] = "PG01" + f"{0.0:14.6f}" * 3 + f"{1.0:14.6f}"
                break
    eph = parse_sp3("\n".join(lines) + "\n")
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    recs = compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), biases)
    g01 = _by_prn(recs, "G01")
    # record kept with reason instead of vanishing
    bad = [r for r in g01 if r.status != "ok"]
    assert bad and any("coordinates" in (r.reason or "") for r in bad)


def test_reject_unsupported_correction_header():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    extra = ("G  0  0" + " " * 46 + "IONOSPHERIC CORR")[:60] + "IONOSPHERIC CORR\n"
    bad = text.replace(" " * 60 + "END OF HEADER",
                       extra + " " * 60 + "END OF HEADER", 1)
    with pytest.raises(RejectedContentError):
        parse_rinex(bad)


def test_invalid_bias_table_rejected():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    with pytest.raises(ValueError):
        compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), {"R01": 1.0})
    with pytest.raises(ValueError):
        compute_ionosphere(obs, eph, Site(30.0, 114.0, 50.0), {"G01": float("nan")})
