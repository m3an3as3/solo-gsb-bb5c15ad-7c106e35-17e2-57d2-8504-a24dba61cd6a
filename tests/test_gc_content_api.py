"""固定窗口 GC 含量范围接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.gccontent import gc_content_windows, window_starts


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def _compute(client, sequence, window_length, step=1, topology="linear"):
    return client.post(
        "/gc-content",
        json={
            "sequence": sequence,
            "topology": topology,
            "window_length": window_length,
            "step": step,
        },
    )


# ---------- 确定碱基与基本滑窗 ----------


def test_determinate_linear_sequence_min_equals_max(client):
    # 全确定碱基：上下界必然相同
    resp = _compute(client, "ATGCATGC", 4, step=2)
    assert resp.status_code == 200
    data = resp.json()
    assert data["sequence_length"] == 8
    assert data["window_length"] == 4
    assert data["step"] == 2
    assert data["window_count"] == 3
    assert data["windows"] == [
        {"index": 0, "start": 0, "end": 4, "length": 4, "wraps_origin": False,
         "min_gc_count": 2, "max_gc_count": 2,
         "min_gc_percent": 50.0, "max_gc_percent": 50.0},
        {"index": 1, "start": 2, "end": 6, "length": 4, "wraps_origin": False,
         "min_gc_count": 2, "max_gc_count": 2,
         "min_gc_percent": 50.0, "max_gc_percent": 50.0},
        {"index": 2, "start": 4, "end": 8, "length": 4, "wraps_origin": False,
         "min_gc_count": 2, "max_gc_count": 2,
         "min_gc_percent": 50.0, "max_gc_percent": 50.0},
    ]


def test_linear_windows_keep_only_full_windows(client):
    # L=5, W=3, step=2：起点 0、2 均完整；末尾位置 4 剩不足一步的短尾，不补窗
    resp = _compute(client, "GCATG", 3, step=2)
    assert resp.status_code == 200
    windows = resp.json()["windows"]
    assert [(w["start"], w["end"]) for w in windows] == [(0, 3), (2, 5)]
    assert all(not w["wraps_origin"] for w in windows)
    assert [w["min_gc_count"] for w in windows] == [2, 1]  # GCA / ATG


def test_step_one_visits_every_start(client):
    resp = _compute(client, "AAAAC", 3, step=1)
    data = resp.json()
    assert data["window_count"] == 3  # L-W+1
    assert [w["start"] for w in data["windows"]] == [0, 1, 2]
    # 只有最后一个窗口含末尾的 C
    assert [w["min_gc_count"] for w in data["windows"]] == [0, 0, 1]


def test_window_equals_sequence_length_single_window_ignores_step(client):
    resp = _compute(client, "GCGC", 4, step=2)
    assert resp.status_code == 200
    data = resp.json()
    assert data["window_count"] == 1
    assert data["windows"][0]["start"] == 0
    assert data["windows"][0]["end"] == 4
    assert data["windows"][0]["min_gc_count"] == 4
    assert data["windows"][0]["min_gc_percent"] == 100.0


def test_all_at_window_zero_percent(client):
    resp = _compute(client, "ATAT", 2)
    windows = resp.json()["windows"]
    assert all(w["min_gc_count"] == w["max_gc_count"] == 0 for w in windows)
    assert all(w["min_gc_percent"] == w["max_gc_percent"] == 0.0 for w in windows)


# ---------- 简并符号的上下界 ----------


def test_degenerate_bases_split_bounds(client):
    # 窗口 S W N：S={G,C} 必然 GC；W={A,T} 必不 GC；N 可能 GC
    # -> 下界 1、上界 2，即 33.33..% ~ 66.66..%
    resp = _compute(client, "SWN", 3)
    w = resp.json()["windows"][0]
    assert w["min_gc_count"] == 1
    assert w["max_gc_count"] == 2
    assert w["min_gc_percent"] == pytest.approx(100.0 / 3, abs=1e-6)
    assert w["max_gc_percent"] == pytest.approx(200.0 / 3, abs=1e-6)
    # 百分比四舍五入到 6 位小数
    assert w["min_gc_percent"] == round(100.0 / 3, 6)
    assert w["max_gc_percent"] == round(200.0 / 3, 6)


@pytest.mark.parametrize(
    "base,certain,possible",
    [
        ("A", False, False),
        ("T", False, False),
        ("C", True, True),
        ("G", True, True),
        ("S", True, True),    # GC：必然
        ("W", False, False),  # AT：必不
        ("R", False, True),   # AG：仅可能
        ("Y", False, True),   # CT：仅可能
        ("K", False, True),   # GT：仅可能
        ("M", False, True),   # AC：仅可能
        ("B", False, True),   # CGT：仅可能
        ("D", False, True),   # AGT：仅可能
        ("H", False, True),   # ACT：仅可能
        ("V", False, True),   # ACG：仅可能
        ("N", False, True),   # ACGT：仅可能
    ],
)
def test_every_iupac_code_contributes_to_correct_bound(base, certain, possible):
    # 混入一个确定的 C 做背景，窗口长度恒为 3
    result = gc_content_windows(f"C{base}A", "linear", 3, 1)
    w = result["windows"][0]
    assert w["min_gc_count"] == 1 + int(certain)
    assert w["max_gc_count"] == 1 + int(possible)


def test_lower_bound_never_exceeds_upper_bound(client):
    resp = _compute(client, "RYSWKMBDHVNACGT", 7, step=3)
    for w in resp.json()["windows"]:
        assert w["min_gc_count"] <= w["max_gc_count"]
        assert w["min_gc_percent"] <= w["max_gc_percent"]


# ---------- 环形序列 ----------


def test_circular_window_wraps_origin(client):
    # L=5, W=3, step=4：起点 0 与 4；起点 4 的窗口取位置 4,0,1 = G,A,T
    resp = _compute(client, "ATATG", 3, step=4, topology="circular")
    assert resp.status_code == 200
    windows = resp.json()["windows"]
    assert len(windows) == 2
    assert windows[0] == {
        "index": 0, "start": 0, "end": 3, "length": 3, "wraps_origin": False,
        "min_gc_count": 0, "max_gc_count": 0,
        "min_gc_percent": 0.0, "max_gc_percent": 0.0,
    }
    assert windows[1] == {
        "index": 1, "start": 4, "end": 2, "length": 3, "wraps_origin": True,
        "min_gc_count": 1, "max_gc_count": 1,
        "min_gc_percent": round(100.0 / 3, 6), "max_gc_percent": round(100.0 / 3, 6),
    }


def test_circular_step_covers_ring_with_full_length_ceiling(client):
    # L=10, W=3, step=4：起点 0,4,8（ceil(10/4)=3 个），最后一个跨原点
    resp = _compute(client, "AAAAAAAAAA", 3, step=4, topology="circular")
    windows = resp.json()["windows"]
    assert [w["start"] for w in windows] == [0, 4, 8]
    assert [w["wraps_origin"] for w in windows] == [False, False, True]
    assert windows[-1]["end"] == 1


def test_circular_step_greater_than_length_single_window(client):
    # step >= L 时环形也只在起点 0 取一个窗口（不取一整圈以上）
    resp = _compute(client, "GCAT", 2, step=10, topology="circular")
    windows = resp.json()["windows"]
    assert len(windows) == 1
    assert windows[0]["start"] == 0
    assert not windows[0]["wraps_origin"]


def test_circular_full_length_window_covers_each_base_once(client):
    # W == L：唯一窗口必须覆盖每个物理位置恰好一次，即使含简并码
    sequence = "N" * 10
    resp = _compute(client, sequence, 10, step=3, topology="circular")
    data = resp.json()
    assert data["window_count"] == 1
    w = data["windows"][0]
    assert w["min_gc_count"] == 0
    assert w["max_gc_count"] == 10  # 不是 20：复制串的窗口切片长度仅 L
    assert not w["wraps_origin"]
    assert w["end"] == 0


def test_circular_degenerate_window_content(client):
    # L=6, W=4, step=5：起点 0 与 5；窗口 5 覆盖位置 5(N),0(G),1(C),2(S)
    # 必然 GC：G,C,S = 3；可能：再加 N = 4
    resp = _compute(client, "GCSATN", 4, step=5, topology="circular")
    by_start = {w["start"]: w for w in resp.json()["windows"]}
    wrapped = by_start[5]
    assert wrapped["wraps_origin"] is True
    assert wrapped["end"] == 3
    assert wrapped["min_gc_count"] == 3
    assert wrapped["max_gc_count"] == 4
    # 起点 0 的窗口 G,C,S,A：必然/可能都是 3
    assert by_start[0]["min_gc_count"] == 3
    assert by_start[0]["max_gc_count"] == 3


def test_linear_and_circular_window_starts_differ():
    assert window_starts(8, 4, 3, circular=False) == [0, 3]
    assert window_starts(8, 4, 3, circular=True) == [0, 3, 6]
    assert window_starts(8, 8, 3, circular=True) == [0]
    assert window_starts(8, 8, 3, circular=False) == [0]


# ---------- 校验与状态码 ----------


def test_window_longer_than_sequence_returns_422(client):
    resp = _compute(client, "ACGT", 5, step=1)
    assert resp.status_code == 422
    assert "window_length 5 exceeds sequence length 4" in resp.text


def test_non_positive_window_or_step_returns_422(client):
    for payload in (
        {"sequence": "ACGT", "window_length": 0, "step": 1},
        {"sequence": "ACGT", "window_length": 4, "step": 0},
        {"sequence": "ACGT", "window_length": -1, "step": 1},
        {"sequence": "ACGT", "window_length": 4, "step": -2},
    ):
        assert client.post("/gc-content", json=payload).status_code == 422, payload


def test_empty_and_illegal_sequence_returns_422(client):
    assert client.post(
        "/gc-content",
        json={"sequence": "", "window_length": 1, "step": 1},
    ).status_code == 422
    resp = client.post(
        "/gc-content",
        json={"sequence": "ACGU", "window_length": 2, "step": 1},
    )
    assert resp.status_code == 422
    assert "non-IUPAC" in resp.text


def test_illegal_topology_returns_422(client):
    resp = client.post(
        "/gc-content",
        json={"sequence": "ACGT", "topology": "plasmid",
              "window_length": 2, "step": 1},
    )
    assert resp.status_code == 422


def test_lowercase_sequence_normalized(client):
    resp = _compute(client, "gcswn", 5)
    assert resp.status_code == 200
    w = resp.json()["windows"][0]
    assert w["min_gc_count"] == 3  # G,C,S
    assert w["max_gc_count"] == 4  # +N；W 不含 G/C


def test_default_topology_is_linear(client):
    resp = client.post(
        "/gc-content",
        json={"sequence": "GCAT", "window_length": 2, "step": 2},
    )
    assert resp.status_code == 200
    assert resp.json()["topology"] == "linear"


# ---------- 只读 ----------


def test_gc_content_is_read_only_and_needs_no_records(client):
    # 先产生一条已保存记录，确认 GC 计算不新增提交
    client.post(
        "/digests",
        json={
            "sequence": "AAAAGAATTCTTTT",
            "topology": "linear",
            "enzymes": [
                {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
            ],
        },
    )
    before = client.get("/digests").json()
    resp = _compute(client, "ATGCSWNATG", 4, step=3, topology="circular")
    assert resp.status_code == 200
    after = client.get("/digests").json()
    assert after == before
