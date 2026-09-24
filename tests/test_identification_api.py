"""未知样本识别接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 回文模式 GC 的一位切口酶（cut 1/1 平端），只在 GC 位点产生切口，
# 便于按目标片段长度精确构造记录
GC_ENZYME = {"name": "GC-cutter", "pattern": "GC", "cut_top": 1, "cut_bottom": 1}
# 识别 GAATTC，在全 A 序列上无位点，用于构造"不切割"的记录
ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}


def _seq_with_cuts(length, cut_positions):
    """构造线性序列：全 A 背景，在每个切口位置 p 处放置 GC（切口在 p-1 与 p 之间）。

    相邻切口至少相隔 2（即片段长度 >= 2），否则 GC 放置会互相覆盖。
    """
    seq = ["A"] * length
    for p in cut_positions:
        seq[p - 1] = "G"
        seq[p] = "C"
    return "".join(seq)


def _create_with_fragments(client, fragment_lengths):
    """创建一条线性记录，使其片段长度恰为给定清单（按物理顺序，每段 >= 2）。"""
    cuts = []
    position = 0
    for frag_len in fragment_lengths[:-1]:
        position += frag_len
        cuts.append(position)
    resp = client.post(
        "/digests",
        json={
            "sequence": _seq_with_cuts(sum(fragment_lengths), cuts),
            "topology": "linear",
            "enzymes": [GC_ENZYME],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _create_uncut(client, length):
    """创建一条无切口的记录，唯一片段为全长分子。"""
    resp = client.post(
        "/digests",
        json={
            "sequence": "A" * length,
            "topology": "linear",
            "enzymes": [ECORI],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _identify(client, digest_ids, observed, tolerance, max_unpaired):
    return client.post(
        "/digest-identification",
        json={
            "digest_ids": digest_ids,
            "observed_fragment_lengths": observed,
            "tolerance": tolerance,
            "max_unpaired": max_unpaired,
        },
    )


# ---------- 命中、排序与回显 ----------


def test_best_candidate_ranked_first_and_marked(client):
    # 候选 A：片段 [100, 200, 300]，观测 305 在容差 10 内 -> 3 配对、差 5、0 未配对
    rec_a = _create_with_fragments(client, [100, 200, 300])
    # 候选 B：片段 [100, 200, 400]，305 无法与 400 配对 -> 2 配对、两侧各 1 未配对
    rec_b = _create_with_fragments(client, [100, 200, 400])
    resp = _identify(
        client, [rec_a["id"], rec_b["id"]], [100, 200, 305], 10, 4
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["digest_ids"] == [rec_a["id"], rec_b["id"]]
    assert data["observed_fragment_lengths"] == [100, 200, 305]
    assert data["tolerance"] == 10
    assert data["max_unpaired"] == 4
    assert data["candidate_count"] == 2
    assert data["hit"] is True
    assert data["ambiguous"] is False
    assert data["best_digest_id"] == rec_a["id"]

    cand_a, cand_b = data["candidates"]
    assert cand_a["digest_id"] == rec_a["id"]
    assert cand_a["fragment_count"] == 3
    assert cand_a["pair_count"] == 3
    assert cand_a["unpaired_count"] == 0
    assert cand_a["diff_sum"] == 5
    assert cand_a["eligible"] is True
    assert cand_a["best"] is True
    assert cand_a["pairs"] == [
        {"observed_index": 0, "expected_index": 0,
         "observed_length": 100, "expected_length": 100, "length_diff": 0},
        {"observed_index": 1, "expected_index": 1,
         "observed_length": 200, "expected_length": 200, "length_diff": 0},
        {"observed_index": 2, "expected_index": 2,
         "observed_length": 305, "expected_length": 300, "length_diff": 5},
    ]
    assert cand_a["unpaired_observed_indices"] == []
    assert cand_a["unpaired_expected_indices"] == []

    assert cand_b["digest_id"] == rec_b["id"]
    assert cand_b["pair_count"] == 2
    assert cand_b["unpaired_count"] == 2
    assert cand_b["diff_sum"] == 0
    assert cand_b["eligible"] is True
    assert cand_b["best"] is False
    assert cand_b["unpaired_observed_indices"] == [2]
    assert cand_b["unpaired_expected_indices"] == [2]


def test_candidates_sorted_by_unpaired_then_diff_then_request_order(client):
    # A：未配对 2、差值 0；B：完美匹配（0 未配对）；C：0 未配对但差值 1
    rec_a = _create_with_fragments(client, [100, 200, 400])
    rec_b = _create_uncut(client, 100)
    rec_c = _create_uncut(client, 99)
    resp = _identify(
        client, [rec_a["id"], rec_b["id"], rec_c["id"]], [100], 2, 6
    )
    assert resp.status_code == 200
    order = [c["digest_id"] for c in resp.json()["candidates"]]
    # B(0,0) -> C(0,1) -> A(2,0)
    assert order == [rec_b["id"], rec_c["id"], rec_a["id"]]


def test_request_order_breaks_score_tie_and_picks_best_representative(client):
    rec_x = _create_uncut(client, 100)
    rec_y = _create_uncut(client, 100)
    resp = _identify(client, [rec_x["id"], rec_y["id"]], [100], 0, 0)
    assert resp.status_code == 200
    data = resp.json()
    assert data["hit"] is True
    assert data["ambiguous"] is True
    # 请求顺序最靠前的并列候选作为代表
    assert data["best_digest_id"] == rec_x["id"]
    assert [c["digest_id"] for c in data["candidates"]] == [
        rec_x["id"], rec_y["id"]
    ]
    assert all(c["best"] is True for c in data["candidates"])

    # 调换请求顺序：代表随之变为另一条
    resp = _identify(client, [rec_y["id"], rec_x["id"]], [100], 0, 0)
    data = resp.json()
    assert data["ambiguous"] is True
    assert data["best_digest_id"] == rec_y["id"]


def test_tolerance_boundary_is_inclusive(client):
    rec = _create_uncut(client, 110)
    # 差值恰好 10 == tolerance，可配对
    resp = _identify(client, [rec["id"]], [100], 10, 0)
    assert resp.status_code == 200
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 1
    assert cand["unpaired_count"] == 0
    assert cand["pairs"][0]["length_diff"] == 10

    # 差值 11 > tolerance，不能配对
    resp = _identify(client, [rec["id"]], [99], 10, 0)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 0
    assert cand["unpaired_count"] == 2
    assert resp.json()["hit"] is False


# ---------- 配对优化与字典序择优 ----------


def test_more_pairs_beats_smaller_cost(client):
    # 观测 [100, 110]、预期 [105, 106]、容差 10：
    # 贪心看似 o0->e0(差5) 后 o1 无法配 e1(差4? 110-106=4 可配 => 2 对差9)
    # 构造迫使"多一对 > 差更小"的场景：
    # 观测 [100, 106]、预期 [100, 106] 容差 6：
    #   o0->e1(差6)+o1 无配对 = 1 对差 6；o0->e0 + o1->e1 = 2 对差 0
    rec = _create_with_fragments(client, [100, 106])
    resp = _identify(client, [rec["id"]], [100, 106], 6, 0)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 2
    assert cand["diff_sum"] == 0

    # 真正的"多配对优先于少差值"：
    # 观测 [100, 110]、预期 [104]、容差 10：
    #   只有一个预期，1 对：o0->e0 差 4 优于 o1->e0 差 6（最小化差值）
    rec2 = _create_uncut(client, 104)
    resp = _identify(client, [rec2["id"]], [100, 110], 10, 2)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 1
    assert cand["diff_sum"] == 4
    assert cand["pairs"][0]["observed_index"] == 0
    assert cand["unpaired_observed_indices"] == [1]

def test_equal_pair_count_minimizes_total_diff(client):
    # 观测 [100, 110]、预期 [102, 108]、容差 10：
    # 顺配差 2+2=4；叉配差 8+8=16 -> 选顺配
    rec = _create_with_fragments(client, [102, 108])
    resp = _identify(client, [rec["id"]], [100, 110], 10, 0)
    cand = resp.json()["candidates"][0]
    assert [(p["observed_index"], p["expected_index"]) for p in cand["pairs"]] == [
        (0, 0), (1, 1)
    ]
    assert cand["diff_sum"] == 4


def test_lexicographic_tiebreak_prefers_smaller_observed_paired(client):
    # 观测 [100, 110]、预期 [105]、容差 10：只有一个预期，两种配对都是
    # 1 对、差值均为 5（100/105、110/105）。配对序列字典序最小意味着
    # 更小编号的观测优先被配出 -> 配对为 o0->e0，观测 1 未配对。
    rec = _create_uncut(client, 105)
    resp = _identify(client, [rec["id"]], [100, 110], 10, 2)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 1
    assert cand["diff_sum"] == 5
    assert cand["pairs"][0]["observed_index"] == 0
    assert cand["pairs"][0]["expected_index"] == 0
    assert cand["unpaired_observed_indices"] == [1]


def test_lexicographic_tiebreak_prefers_smaller_expected_index(client):
    # 观测 [105]、预期 [105, 105]、容差 0：必配时占用预期索引 0
    rec = _create_with_fragments(client, [105, 105])
    resp = _identify(client, [rec["id"]], [105], 0, 2)
    cand = resp.json()["candidates"][0]
    assert cand["pairs"][0]["expected_index"] == 0
    assert cand["unpaired_expected_indices"] == [1]


def test_repeated_lengths_counted_as_distinct_fragments(client):
    # 观测两个 100，预期两个 100：分别配对，2 对而非按多重集折叠后 1 对
    rec = _create_with_fragments(client, [100, 100])
    resp = _identify(client, [rec["id"]], [100, 100], 0, 0)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 2
    assert cand["unpaired_count"] == 0
    assert [p["observed_index"] for p in cand["pairs"]] == [0, 1]
    assert [p["expected_index"] for p in cand["pairs"]] == [0, 1]


def test_repeated_observed_one_more_than_expected_leaves_highest_index(client):
    # 观测三个 100、预期两个 100：字典序最小方案让观测索引最大者（2）未配对
    rec = _create_with_fragments(client, [100, 100])
    resp = _identify(client, [rec["id"]], [100, 100, 100], 0, 2)
    cand = resp.json()["candidates"][0]
    assert cand["pair_count"] == 2
    assert cand["unpaired_observed_indices"] == [2]
    assert cand["unpaired_expected_indices"] == []


def test_unpaired_indices_returned_sorted(client):
    # 观测多于预期：4 个观测对 2 个预期（全部等长）-> 观测索引 2、3 未配对
    rec = _create_with_fragments(client, [50, 60])
    resp = _identify(client, [rec["id"]], [50, 60, 51, 61], 0, 4)
    cand = resp.json()["candidates"][0]
    assert cand["pairs"] and cand["pair_count"] == 2
    assert cand["unpaired_observed_indices"] == [2, 3]
    assert cand["unpaired_expected_indices"] == []

    # 预期多于观测：预期索引 2、3 未配对
    rec2 = _create_with_fragments(client, [50, 60, 70, 80])
    resp = _identify(client, [rec2["id"]], [50, 60], 0, 4)
    cand = resp.json()["candidates"][0]
    assert cand["unpaired_observed_indices"] == []
    assert cand["unpaired_expected_indices"] == [2, 3]


def test_single_candidate_and_single_observed_supported(client):
    rec = _create_uncut(client, 42)
    resp = _identify(client, [rec["id"]], [42], 0, 0)
    assert resp.status_code == 200
    data = resp.json()
    assert data["candidate_count"] == 1
    assert data["hit"] is True
    assert data["ambiguous"] is False
    assert data["best_digest_id"] == rec["id"]
    assert data["candidates"][0]["unpaired_count"] == 0


def test_eight_candidates_and_twelve_fragments_supported(client):
    recs = [_create_uncut(client, 50) for _ in range(8)]
    resp = _identify(client, [r["id"] for r in recs], [50] * 12, 0, 24)
    assert resp.status_code == 200
    data = resp.json()
    assert data["candidate_count"] == 8
    # 每个候选 1 个预期对 12 个观测：1 配对、11 观测未配对
    assert all(c["pair_count"] == 1 for c in data["candidates"])
    assert all(c["unpaired_count"] == 11 for c in data["candidates"])


# ---------- 命中门槛与歧义 ----------


def test_no_hit_when_all_candidates_exceed_max_unpaired(client):
    rec_a = _create_with_fragments(client, [100, 200, 400])
    rec_b = _create_uncut(client, 900)
    resp = _identify(
        client, [rec_a["id"], rec_b["id"]], [100, 200, 305], 10, 0
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["hit"] is False
    assert data["ambiguous"] is False
    assert data["best_digest_id"] is None
    # 候选仍全部返回并按分数排序，但均不计最优
    assert [c["digest_id"] for c in data["candidates"]] == [
        rec_a["id"], rec_b["id"]
    ]
    assert all(c["eligible"] is False for c in data["candidates"])
    assert all(c["best"] is False for c in data["candidates"])


def test_max_unpaired_boundary_is_inclusive(client):
    rec = _create_with_fragments(client, [100, 200, 400])
    # 观测 [100, 200, 305]、容差 10：2 配对，两侧各 1 未配对 = 2
    resp = _identify(client, [rec["id"]], [100, 200, 305], 10, 2)
    data = resp.json()
    assert data["hit"] is True
    assert data["candidates"][0]["eligible"] is True

    # 上限 1 时同一结果不命中
    resp = _identify(client, [rec["id"]], [100, 200, 305], 10, 1)
    data = resp.json()
    assert data["hit"] is False
    assert data["candidates"][0]["eligible"] is False


def test_ambiguity_only_counts_eligible_candidates(client):
    # A、B 并列最优且均命中；C 分数相同但未配对超门槛 -> 不参与歧义
    rec_a = _create_uncut(client, 100)
    rec_b = _create_uncut(client, 100)
    rec_c = _create_with_fragments(client, [100, 100])
    resp = _identify(
        client, [rec_a["id"], rec_b["id"], rec_c["id"]], [100], 0, 0
    )
    data = resp.json()
    assert data["hit"] is True
    assert data["ambiguous"] is True
    eligible = [c for c in data["candidates"] if c["eligible"]]
    assert {c["digest_id"] for c in eligible} == {rec_a["id"], rec_b["id"]}
    assert all(c["best"] is True for c in eligible)
    assert not any(
        c["best"] for c in data["candidates"] if not c["eligible"]
    )


def test_non_ambiguous_when_scores_differ(client):
    rec_a = _create_uncut(client, 100)
    rec_b = _create_uncut(client, 101)
    resp = _identify(client, [rec_a["id"], rec_b["id"]], [100], 5, 0)
    data = resp.json()
    assert data["hit"] is True
    assert data["ambiguous"] is False
    assert data["best_digest_id"] == rec_a["id"]
    best_flags = [c["best"] for c in data["candidates"]]
    assert best_flags == [True, False]


# ---------- 校验与状态码 ----------


def test_duplicate_digest_ids_return_422(client):
    rec = _create_uncut(client, 50)
    resp = _identify(client, [rec["id"], rec["id"]], [50], 0, 0)
    assert resp.status_code == 422
    assert "duplicate digest ids" in resp.text


def test_candidate_count_limits_return_422(client):
    rec = _create_uncut(client, 50)
    # 少于 1 个
    resp = _identify(client, [], [50], 0, 0)
    assert resp.status_code == 422
    # 多于 8 个（校验先于记录存在性检查）
    resp = _identify(client, list(range(1, 10)), [50], 0, 0)
    assert resp.status_code == 422


def test_observed_length_limits_return_422(client):
    rec = _create_uncut(client, 50)
    # 空观测列表
    assert _identify(client, [rec["id"]], [], 0, 0).status_code == 422
    # 超过 12 个观测
    assert _identify(
        client, [rec["id"]], [50] * 13, 0, 0
    ).status_code == 422
    # 非正长度
    assert _identify(client, [rec["id"]], [0], 0, 0).status_code == 422
    assert _identify(client, [rec["id"]], [-5], 0, 0).status_code == 422


def test_negative_tolerance_and_max_unpaired_return_422(client):
    rec = _create_uncut(client, 50)
    assert _identify(client, [rec["id"]], [50], -1, 0).status_code == 422
    assert _identify(client, [rec["id"]], [50], 0, -1).status_code == 422


def test_missing_fields_return_422(client):
    rec = _create_uncut(client, 50)
    base = {
        "digest_ids": [rec["id"]],
        "observed_fragment_lengths": [50],
        "tolerance": 0,
        "max_unpaired": 0,
    }
    for field in ("digest_ids", "observed_fragment_lengths", "tolerance", "max_unpaired"):
        payload = {k: v for k, v in base.items() if k != field}
        resp = client.post("/digest-identification", json=payload)
        assert resp.status_code == 422, field


def test_missing_record_returns_404(client):
    rec = _create_uncut(client, 50)
    resp = _identify(client, [rec["id"], 999], [50], 0, 0)
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_first_missing_id_in_request_order_is_reported(client):
    resp = _identify(client, [999, 998], [50], 0, 0)
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_validation_precedes_existence_check(client):
    # 重复 id 直接 422，即使 id 不存在也不报 404
    assert _identify(client, [999, 999], [50], 0, 0).status_code == 422


def test_record_with_more_than_12_fragments_returns_422(client):
    # 13 个片段（每段 >= 2，全长 26）
    rec = _create_with_fragments(client, [2] * 13)
    ok = _create_uncut(client, 50)
    resp = _identify(client, [ok["id"], rec["id"]], [2], 0, 0)
    assert resp.status_code == 422
    assert str(rec["id"]) in resp.json()["detail"]
    assert "12" in resp.json()["detail"]

    # 恰好 12 个片段合法
    rec12 = _create_with_fragments(client, [2] * 12)
    resp = _identify(client, [rec12["id"]], [2], 0, 22)
    assert resp.status_code == 200


def test_404_precedes_fragment_limit_check(client):
    # 片段超限检查在记录取到之后：缺失记录先报 404
    missing = 999
    rec = _create_with_fragments(client, [2] * 13)
    resp = _identify(client, [missing, rec["id"]], [2], 0, 0)
    assert resp.status_code == 404


# ---------- 只读 ----------


def test_identification_is_read_only(client):
    rec_a = _create_with_fragments(client, [100, 200, 300])
    rec_b = _create_with_fragments(client, [100, 200, 400])
    snapshot = client.get("/digests").json()
    resp = _identify(
        client, [rec_a["id"], rec_b["id"]], [100, 200, 305], 10, 4
    )
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{rec_a['id']}").json() == rec_a
    assert client.get(f"/digests/{rec_b['id']}").json() == rec_b
    assert client.get("/digests").json() == snapshot
