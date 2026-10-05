import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import app


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_synthetic.py")],
                   check=True, cwd=ROOT)


def test_position_endpoint():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/position",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)})
    f1.close()
    f2.close()
    assert r.status_code == 200
    body = r.json()
    assert body["n_epochs"] == 8 and body["n_ok"] == 8
    ep = body["epochs"][0]
    assert ep["status"] == "ok"
    assert abs(ep["geodetic"]["lat_deg"] - 30.0) < 1e-4
    assert abs(ep["receiver_clock_s"] - 1.5e-4) < 1e-6
    assert len(ep["residuals_m"]) == 8


def test_reject_garbage():
    client = TestClient(app)
    r = client.post("/position",
                    files={"rinex": ("a.rnx", b"garbage"),
                           "sp3": ("b.sp3", b"garbage")})
    assert r.status_code == 422
    assert "rinex" in r.json()["detail"]


def test_ionosphere_endpoint():
    import json
    client = TestClient(app)
    biases = (ROOT / "examples" / "biases.json").read_text()
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/ionosphere",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data={"lat_deg": "30", "lon_deg": "114", "height_m": "50",
                          "biases_ns": biases})
    f1.close(); f2.close()
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_epochs"] == 8 and body["n_records"] == 64
    assert body["n_ok"] > 0
    ok = [x for x in body["records"] if x["status"] == "ok"]
    for x in ok:
        assert x["pierce_lat_deg"] is not None and x["elevation_deg"] >= 10.0
    g08 = [x for x in body["records"] if x["prn"] == "G08"]
    assert all(x["status"] == "invalid" and "bias" in x["reason"] for x in g08)
    negs = [x for x in body["records"]
            if x["prn"] == "G07" and x["slant_tec_tceu"] is not None]
    assert any(x["slant_tec_tceu"] < 0 for x in negs)


def test_ionosphere_bad_bias_json():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/ionosphere",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data={"lat_deg": "30", "lon_deg": "114", "height_m": "50",
                          "biases_ns": "not-json"})
    f1.close(); f2.close()
    assert r.status_code == 400


def test_ionosphere_bad_coordinates():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/ionosphere",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data={"lat_deg": "999", "lon_deg": "0", "height_m": "50",
                          "biases_ns": "{}\"G01\":1}"})
    f1.close(); f2.close()
    assert r.status_code == 400
