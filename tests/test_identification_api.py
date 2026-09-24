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


def _create(client, sequence, enzymes=None):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": "linear",
            "enzymes": enzymes if enzymes is not None else [GC_ENZYME],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _create_with_fragments(client, fragment_lengths):
    """创建一条线性记录，使其片段长度恰为给定清单（按物理顺序，每段 >= 2）。"""
    cuts = []
    position = 0
    for frag_len in fragment_lengths[:-1]:
        position += frag_len
        cuts.append(position)
    return _create(client, _seq_with_cuts(sum(fragment_lengths), cuts))


def _create_uncut(client, length):
    """创建一条无切口的记录，唯一片段为全长分子。"""
    return _create(client, "A" * length, enzymes=[ECORI])


def _identify(client, digest_ids, observed, tolerance=0, max_unmatched=0):
    return client.post(
        "/digest-identification",
        json={
            "digest_ids": digest_ids,
            "observed_lengths": observed,
            "tolerance": tolerance,
            "max_unmatched": max_unmatched,
        },
    )


# ---------- 配对与得分 ----------


def test_exact_full_match_is_unique_best_hit(client):
    record = _create_with_fragments(client, [10, 20, 30])
    other = _create_with_fragments(client, [10, 20, 40])
    resp = _identify(client, [record["id"], other["id"]], [10, 20, 30])
    assert resp.status_code == 200
    data = resp.json()
    assert data["digest_ids"] == [record["id"], other["id"]]
    assert data["observed_lengths"] == [10, 20, 30]
    assert data["tolerance"] == 0
    assert data["max_unmatched"] == 0
    assert data["observed_count"] == 3
    assert data["candidate_count"] == 2

    best = data["candidates"][0]
    assert best["digest_id"] == record["id"]
    assert best["fragment_lengths"] == [10, 20, 30]
    assert best["matched_count"] == 3
    assert best["unmatched_count"] == 0
    assert best["total_difference"] == 0
    assert best["pairs"] == [
        {"observed_index": 0, "fragment_index": 0, "observed_length": 10,
         "expected_length": 10, "difference": 0},
        {"observed_index": 1, "fragment_index": 1, "observed_length": 20,
         "expected_length": 20, "difference": 0},
        {"observed_index": 2, "fragment_index": 2, "observed_length": 30,
         "expected_length": 30, "difference": 0},
    ]
    assert best["unmatched_observed_indices"] == []
    assert best["unmatched_fragment_indices"] == []

    second = data["candidates"][1]
    assert second["digest_id"] == other["id"]
    assert second["matched_count"] == 2
    assert second["unmatched_count"] == 2  # 观测 30 与预期 40 均未配对
    assert second["unmatched_observed_indices"] == [2]
    assert second["unmatched_fragment_indices"] == [2]

    assert data["hit_digest_ids"] == [record["id"]]
    assert data["hit_count"] == 1
    assert data["best_digest_ids"] == [record["id"]]
    assert data["ambiguous"] is False


def test_tolerance_boundary_is_inclusive(client):
    record = _create_uncut(client, 10)
    # 差值恰好等于容差：配对成立
    resp = _identify(client, [record["id"]], [13], tolerance=3, max_unmatched=0)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert candidate["pairs"] == [
        {"observed_index": 0, "fragment_index": 0, "observed_length": 13,
         "expected_length": 10, "difference": 3}
    ]
    assert candidate["unmatched_count"] == 0
    assert resp.json()["hit_digest_ids"] == [record["id"]]

    # 差值超出容差 1：不配对，两侧各计一个未配对
    resp = _identify(client, [record["id"]], [14], tolerance=3, max_unmatched=2)
    candidate = resp.json()["candidates"][0]
    assert candidate["pairs"] == []
    assert candidate["matched_count"] == 0
    assert candidate["unmatched_count"] == 2
    assert candidate["unmatched_observed_indices"] == [0]
    assert candidate["unmatched_fragment_indices"] == [0]


def test_more_pairs_beat_smaller_difference_sum(client):
    # 观测 [10, 12]、预期 [10, 13]、容差 1：
    # 方案 A：10->10（差 0）且 12->13（差 1），2 对、差值总和 1；
    # 方案 B：仅 10->10，1 对、差值总和 0。先最大化配对数，A 胜出。
    record = _create_with_fragments(client, [10, 13])
    resp = _identify(client, [record["id"]], [10, 12], tolerance=1)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert candidate["matched_count"] == 2
    assert candidate["total_difference"] == 1
    assert candidate["unmatched_count"] == 0
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (0, 0),
        (1, 1),
    ]


def test_min_difference_sum_breaks_pair_count_tie(client):
    # 观测 [11, 10]、预期 [10]、容差 1：两种配对都是 1 对，
    # 11->10（差 1）与 10->10（差 0），差值总和最小者胜出
    record = _create_uncut(client, 10)
    resp = _identify(client, [record["id"]], [11, 10], tolerance=1)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert candidate["matched_count"] == 1
    assert candidate["total_difference"] == 0
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (1, 0)
    ]
    assert candidate["unmatched_observed_indices"] == [0]
    assert candidate["unmatched_count"] == 1


def test_lexicographic_tie_break_prefers_earlier_observed(client):
    # 观测 [3, 5]、预期 [4]、容差 1：3->4 与 5->4 的配对数、差值总和
    # 均相同，按（观测索引, 片段索引）字典序固定为 3->4
    record = _create_uncut(client, 4)
    resp = _identify(client, [record["id"]], [3, 5], tolerance=1)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (0, 0)
    ]
    assert candidate["unmatched_observed_indices"] == [1]


def test_lexicographic_tie_break_prefers_earlier_fragment(client):
    # 观测 [15]、预期 [14, 16]、容差 1：15->14 与 15->16 同分，
    # 字典序固定为片段索引更小者
    record = _create_with_fragments(client, [14, 16])
    resp = _identify(client, [record["id"]], [15], tolerance=1)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (0, 0)
    ]
    assert candidate["unmatched_fragment_indices"] == [1]
    assert candidate["unmatched_count"] == 1


def test_duplicate_lengths_count_as_distinct_fragments(client):
    # 预期 [5, 5]、观测 [5]：重复长度按不同片段计数，只有一个 5 能配对；
    # 同分按字典序固定为片段 0，片段 1 未配对
    record = _create_with_fragments(client, [5, 5])
    resp = _identify(client, [record["id"]], [5], tolerance=0, max_unmatched=1)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert candidate["matched_count"] == 1
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (0, 0)
    ]
    assert candidate["unmatched_fragment_indices"] == [1]

    # 观测 [5, 5] 时两个重复长度各自配对，同分按字典序固定为 (0,0) (1,1)
    resp = _identify(client, [record["id"]], [5, 5], tolerance=0)
    candidate = resp.json()["candidates"][0]
    assert candidate["matched_count"] == 2
    assert [(p["observed_index"], p["fragment_index"]) for p in candidate["pairs"]] == [
        (0, 0),
        (1, 1),
    ]
    assert candidate["unmatched_count"] == 0


def test_unmatched_total_counts_both_sides(client):
    record = _create_with_fragments(client, [10, 20, 30])
    resp = _identify(client, [record["id"]], [10, 99], tolerance=0, max_unmatched=3)
    assert resp.status_code == 200
    candidate = resp.json()["candidates"][0]
    assert candidate["matched_count"] == 1
    # 观测侧 1 个未配对（99），预期侧 2 个未配对（20、30）
    assert candidate["unmatched_observed_indices"] == [1]
    assert candidate["unmatched_fragment_indices"] == [1, 2]
    assert candidate["unmatched_count"] == 3


# ---------- 候选排序、命中与歧义 ----------


def test_candidates_sorted_by_unmatched_then_difference(client):
    exact = _create_with_fragments(client, [10, 20])      # 得分 (0, 0)
    near = _create_with_fragments(client, [10, 21])       # 得分 (0, 1)
    short = _create_with_fragments(client, [10, 25])      # 得分 (2, 0)
    ids = [short["id"], near["id"], exact["id"]]
    resp = _identify(client, ids, [10, 20], tolerance=1, max_unmatched=2)
    assert resp.status_code == 200
    data = resp.json()
    assert [c["digest_id"] for c in data["candidates"]] == [
        exact["id"],
        near["id"],
        short["id"],
    ]
    assert [(c["unmatched_count"], c["total_difference"]) for c in data["candidates"]] == [
        (0, 0),
        (0, 1),
        (2, 0),
    ]
    assert data["hit_digest_ids"] == [exact["id"], near["id"], short["id"]]
    assert data["best_digest_ids"] == [exact["id"]]
    assert data["ambiguous"] is False


def test_equal_scores_follow_request_order(client):
    first = _create_with_fragments(client, [10, 20])
    second = _create_with_fragments(client, [10, 20])
    # 两个候选得分完全相同：候选顺序与 best 列表均按请求 ID 顺序
    resp = _identify(client, [second["id"], first["id"]], [10, 20])
    assert resp.status_code == 200
    data = resp.json()
    assert [c["digest_id"] for c in data["candidates"]] == [
        second["id"],
        first["id"],
    ]
    assert data["best_digest_ids"] == [second["id"], first["id"]]
    assert data["ambiguous"] is True


def test_ambiguous_only_when_best_score_tied(client):
    exact = _create_with_fragments(client, [10, 20])   # 得分 (0, 0)
    near = _create_with_fragments(client, [10, 21])    # 得分 (0, 1)
    resp = _identify(
        client, [near["id"], exact["id"]], [10, 20], tolerance=1, max_unmatched=0
    )
    assert resp.status_code == 200
    data = resp.json()
    # 两个候选都是命中，但最优得分唯一 -> 不歧义
    assert data["hit_count"] == 2
    assert data["best_digest_ids"] == [exact["id"]]
    assert data["ambiguous"] is False


def test_no_candidate_within_max_unmatched_returns_no_hit(client):
    record = _create_with_fragments(client, [10, 20])
    # 观测 [10]：预期 20 未配对，未配对总数 1 > max_unmatched 0
    resp = _identify(client, [record["id"]], [10], tolerance=0, max_unmatched=0)
    assert resp.status_code == 200
    data = resp.json()
    assert data["candidate_count"] == 1
    assert data["candidates"][0]["unmatched_count"] == 1
    assert data["hit_digest_ids"] == []
    assert data["hit_count"] == 0
    assert data["best_digest_ids"] == []
    assert data["ambiguous"] is False

    # 放宽到 1 后同一候选成为命中
    resp = _identify(client, [record["id"]], [10], tolerance=0, max_unmatched=1)
    data = resp.json()
    assert data["hit_digest_ids"] == [record["id"]]
    assert data["best_digest_ids"] == [record["id"]]


def test_single_record_and_eight_records_supported(client):
    single = _create_uncut(client, 50)
    resp = _identify(client, [single["id"]], [50])
    assert resp.status_code == 200
    assert resp.json()["candidate_count"] == 1
    assert resp.json()["best_digest_ids"] == [single["id"]]

    records = [_create_uncut(client, 50) for _ in range(8)]
    resp = _identify(client, [r["id"] for r in records], [50])
    assert resp.status_code == 200
    data = resp.json()
    assert data["candidate_count"] == 8
    assert data["hit_count"] == 8
    assert data["ambiguous"] is True
    assert data["best_digest_ids"] == [r["id"] for r in records]


# ---------- 校验与状态码 ----------


def test_missing_record_returns_404(client):
    record = _create_uncut(client, 50)
    resp = _identify(client, [record["id"], 999], [50])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_first_missing_id_in_request_order_is_reported(client):
    resp = _identify(client, [999, 998], [50])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_duplicate_digest_ids_return_422(client):
    record = _create_uncut(client, 50)
    resp = _identify(client, [record["id"], record["id"]], [50])
    assert resp.status_code == 422
    assert "duplicate digest ids" in resp.text


def test_digest_id_count_limits_return_422(client):
    record = _create_uncut(client, 50)
    # 空列表
    assert _identify(client, [], [50]).status_code == 422
    # 多于 8 条（校验先于记录存在性检查）
    assert _identify(client, list(range(1, 10)), [50]).status_code == 422


def test_observed_lengths_validation_returns_422(client):
    record = _create_uncut(client, 50)
    # 空观测列表
    assert _identify(client, [record["id"]], []).status_code == 422
    # 多于 12 个观测片段
    assert _identify(client, [record["id"]], [10] * 13).status_code == 422
    # 非正整数长度
    assert _identify(client, [record["id"]], [0]).status_code == 422
    assert _identify(client, [record["id"]], [-5]).status_code == 422
    # 非整数长度
    resp = client.post(
        "/digest-identification",
        json={
            "digest_ids": [record["id"]],
            "observed_lengths": ["abc"],
            "tolerance": 0,
            "max_unmatched": 0,
        },
    )
    assert resp.status_code == 422


def test_negative_tolerance_or_max_unmatched_return_422(client):
    record = _create_uncut(client, 50)
    assert _identify(client, [record["id"]], [50], tolerance=-1).status_code == 422
    assert (
        _identify(client, [record["id"]], [50], max_unmatched=-1).status_code == 422
    )


def test_record_with_more_than_12_fragments_returns_422(client):
    record = _create_with_fragments(client, [2] * 13)  # 13 个片段
    resp = _identify(client, [record["id"]], [2], max_unmatched=24)
    assert resp.status_code == 422
    assert "13" in resp.json()["detail"]


def test_validation_precedes_existence_check(client):
    # 重复 id 直接 422，即使 id 不存在也不报 404
    assert _identify(client, [999, 999], [50]).status_code == 422


# ---------- 只读 ----------


def test_identification_is_read_only(client):
    record = _create_with_fragments(client, [10, 20])
    other = _create_with_fragments(client, [10, 21])
    resp = _identify(
        client, [record["id"], other["id"]], [10, 20], tolerance=1, max_unmatched=0
    )
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{record['id']}").json() == record
    assert client.get(f"/digests/{other['id']}").json() == other
    assert len(client.get("/digests").json()) == 2
