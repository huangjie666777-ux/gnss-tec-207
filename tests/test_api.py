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
    assert body["n_epochs"] == 5 and body["n_ok"] == 5
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
