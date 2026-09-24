"""PCR 扩增产物预测接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


# 长引物 + 纯 T 填充背景可避免 1-错配位点干扰，便于构造唯一位点模板
FWD = "ATGCGATCGC"          # 10
REV = "GCAGCTAA"            # 8，其反向互补 motif = TTAGCTGC
MOTIF = "TTAGCTGC"
KF, KR = 10, 8


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def _create(client, sequence, topology="linear"):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": topology,
            "enzymes": [
                {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
            ],
        },
    )
    assert resp.status_code == 201
    return resp.json()["id"]


def _predict(client, digest_id, forward, reverse, max_len, expect=200):
    # 测试便捷语义：max_len=None 时取模板长度作为上限（即不限制产物长度）
    if max_len is None:
        record = client.get(f"/digests/{digest_id}").json()
        max_len = record["sequence_length"]
    resp = client.post(
        f"/digests/{digest_id}/pcr-products",
        json={
            "forward_primer": forward,
            "reverse_primer": reverse,
            "max_product_length": max_len,
        },
    )
    assert resp.status_code == expect, resp.text
    return resp.json() if expect == 200 else resp


# ---------- 线性基本场景 ----------


def test_linear_exact_pair_single_product(client):
    # fwd@0，反向 motif@20，产物 [0,28)
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    data = _predict(client, digest_id, FWD, REV, None)
    assert data["id"] == digest_id
    assert data["topology"] == "linear"
    assert data["sequence_length"] == 28
    assert data["forward_primer"] == FWD
    assert data["reverse_primer"] == REV
    assert data["max_product_length"] == 28
    assert data["product_count"] == 1
    assert data["products"] == [
        {
            "sequence": seq,
            "start": 0,
            "end": 28,
            "length": 28,
            "wraps_origin": False,
            "forward_binding": {"start": 0, "end": 10, "mismatches": 0},
            "reverse_binding": {"start": 20, "end": 28, "mismatches": 0},
        }
    ]


def test_product_includes_primer_sequences_at_both_ends(client):
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    product = _predict(client, digest_id, FWD, REV, None)["products"][0]
    # 产物 5' 端即正引物；3' 端为反引物反向互补后的 motif
    assert product["sequence"].startswith(FWD)
    assert product["sequence"].endswith(MOTIF)


def test_multiple_pairs_sorted_by_start_then_length(client):
    # 两个正向位点（0、10），一个反向位点（20）；产物按起点升序
    seq = FWD + FWD + MOTIF
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert [(p["start"], p["length"]) for p in products] == [(0, 28), (10, 18)]


def test_max_product_length_filters(client):
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, FWD, REV, 27)["products"] == []
    products = _predict(client, digest_id, FWD, REV, 28)["products"]
    assert len(products) == 1 and products[0]["length"] == 28


def test_no_facing_pair_returns_empty_array(client):
    # 反向位点在正向位点之前（背向）：不产生产物
    seq = MOTIF + "T" * 10 + FWD
    digest_id = _create(client, seq)
    data = _predict(client, digest_id, FWD, REV, None)
    assert data["products"] == []
    assert data["product_count"] == 0


def test_linear_product_cannot_extend_past_boundary(client):
    # motif 在末尾只露出 7 个碱基（不完整），线性扫描无完整窗口 -> 无产物
    seq = FWD + "T" * 11 + MOTIF[:7]
    digest_id = _create(client, seq)
    data = _predict(client, digest_id, FWD, REV, None)
    assert data["products"] == []


def test_overlapping_binding_regions_excluded(client):
    # 反向位点起点 r < f + kf 时结合区重叠：fwd@0 [0,10)，motif@8 [8,16)
    seq = FWD[:8] + MOTIF  # 0..7 + motif@8，共 16
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, FWD, REV, None)["products"] == []


def test_adjacent_non_overlapping_regions_allowed(client):
    # r == f + kf：结合区恰好相邻不重叠，产物长度 = kf + kr = 18
    seq = FWD + MOTIF
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert len(products) == 1
    assert products[0]["start"] == 0
    assert products[0]["end"] == 18
    assert products[0]["length"] == 18
    assert products[0]["forward_binding"] == {"start": 0, "end": 10, "mismatches": 0}
    assert products[0]["reverse_binding"] == {"start": 10, "end": 18, "mismatches": 0}


# ---------- 错配与 3' 端规则 ----------


def test_single_internal_mismatch_allowed(client):
    # idx4 G->T：与引物 G 无交集的一处内部错配（idx4 在末两位 idx8/9 之外）
    site_mm = FWD[:4] + "T" + FWD[5:]
    seq = site_mm + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert len(products) == 1
    assert products[0]["forward_binding"]["mismatches"] == 1
    assert products[0]["reverse_binding"]["mismatches"] == 0
    # 产物序列取模板原始碱基，不做替换
    assert products[0]["sequence"][:KF] == site_mm


def test_two_mismatches_rejected(client):
    # idx4 G->T、idx5 A->C：两处均与引物碱基（G、A）无交集
    site_2mm = FWD[:4] + "TC" + FWD[6:]
    seq = site_2mm + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, FWD, REV, None)["products"] == []


def test_forward_terminal_two_bases_must_match(client):
    # idx9（3' 末位）C->T 不兼容
    site_last = FWD[:9] + "T"
    seq = site_last + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, FWD, REV, None)["products"] == []
    # idx8（3' 倒数第二位）G->T 不兼容
    site_penult = FWD[:8] + "TT"
    seq2 = site_penult + "T" * 10 + MOTIF
    digest_id2 = _create(client, seq2)
    assert _predict(client, digest_id2, FWD, REV, None)["products"] == []


def test_reverse_terminal_two_bases_must_match(client):
    # 反引物 3' 端对应反向互补 motif 的开头两位（motif = TTAGCTGC）
    # motif[0] T->G（引物 3' 末位不兼容）
    site_last = "G" + MOTIF[1:]
    seq = FWD + "T" * 10 + site_last
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, FWD, REV, None)["products"] == []
    # motif[1] T->G（引物 3' 倒数第二位不兼容）
    site_penult = MOTIF[0] + "G" + MOTIF[2:]
    seq2 = FWD + "T" * 10 + site_penult
    digest_id2 = _create(client, seq2)
    assert _predict(client, digest_id2, FWD, REV, None)["products"] == []


def test_reverse_primer_five_prime_mismatch_allowed(client):
    # motif 末位（= 引物 5' 端）C->A 一处错配允许
    site_5p = MOTIF[:7] + "A"
    seq = FWD + "T" * 10 + site_5p
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert len(products) == 1
    assert products[0]["reverse_binding"]["mismatches"] == 1


# ---------- IUPAC ----------


def test_iupac_degenerate_matching(client):
    # 引物 idx4 用 R={A,G}，模板该位为 G：集合有交集，0 错配
    primer = FWD[:4] + "R" + FWD[5:]
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, primer, REV, None)["products"]
    assert len(products) == 1
    assert products[0]["forward_binding"]["mismatches"] == 0


def test_iupac_disjoint_sets_count_as_mismatch(client):
    # 引物 idx4 用 W={A,T}，模板该位为 G：集合无交集，记一处内部错配（仍允许）
    primer_one = FWD[:4] + "W" + FWD[5:]
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, primer_one, REV, None)["products"]
    assert products[0]["forward_binding"]["mismatches"] == 1


def test_iupac_two_disjoint_positions_rejected(client):
    # 引物 idx4 W 对模板 G、idx5 A 对模板 C：两处无交集（idx4/5 均为内部位置）
    primer_two = FWD[:4] + "WA" + FWD[6:]
    template_site = FWD[:5] + "C" + FWD[6:]  # idx5 A->C
    seq = template_site + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    assert _predict(client, digest_id, primer_two, REV, None)["products"] == []


def test_n_primer_matches_anywhere_with_zero_mismatches(client):
    # 正引物 NNNN 可与任意四聚体按 0 错配结合；产物仍需相向且结合区不重叠
    seq = "ATGC" + "TTTT" + "GCAT"  # L=12
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, "NNNN", "NNNN", None)["products"]
    keys = [(p["start"], p["length"]) for p in products]
    assert products
    assert keys == sorted(keys)
    assert len(keys) == len(set(keys))
    assert all(p["forward_binding"]["mismatches"] == 0 for p in products)
    assert all(p["reverse_binding"]["mismatches"] == 0 for p in products)


# ---------- 环形 ----------


def _place(chars, segment, at, length):
    for i, base in enumerate(segment):
        chars[(at + i) % length] = base


def test_circular_product_wraps_origin(client):
    # L=24：fwd@18（18..27 跨原点，取模 end=4），motif@6；d=(6-18)%24=12，长度 20
    chars = ["T"] * 24
    _place(chars, FWD, 18, 24)
    _place(chars, MOTIF, 6, 24)
    seq = "".join(chars)
    digest_id = _create(client, seq, topology="circular")
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert len(products) == 1
    p = products[0]
    assert p["wraps_origin"] is True
    assert p["start"] == 18
    assert p["end"] == 14
    assert p["length"] == 20
    assert p["sequence"] == seq[18:] + seq[:14]
    assert p["forward_binding"] == {"start": 18, "end": 4, "mismatches": 0}
    assert p["reverse_binding"] == {"start": 6, "end": 14, "mismatches": 0}


def test_circular_nonwrapping_product(client):
    # fwd@0、motif@10，d=10 不跨原点，长度 18；模板尾部额外填充使 L=28
    seq = FWD + MOTIF + "T" * 10
    digest_id = _create(client, seq, topology="circular")
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert len(products) == 1
    p = products[0]
    assert p["wraps_origin"] is False
    assert p["end"] == 18
    assert p["length"] == 18


def test_circular_product_never_longer_than_one_turn(client):
    chars = ["T"] * 24
    _place(chars, FWD, 18, 24)
    _place(chars, MOTIF, 6, 24)
    seq = "".join(chars)
    digest_id = _create(client, seq, topology="circular")
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    assert all(p["length"] <= 24 for p in products)


def test_circular_full_length_product(client):
    # L=18：fwd(8 聚)@8（8..15），motif@0（0..7）；d=10，长度恰为 18（整圈）
    fwd8 = FWD[:8]
    chars = ["T"] * 18
    _place(chars, fwd8, 8, 18)
    _place(chars, MOTIF, 0, 18)
    seq = "".join(chars)
    digest_id = _create(client, seq, topology="circular")
    products = _predict(client, digest_id, fwd8, REV, None)["products"]
    assert len(products) == 1
    p = products[0]
    assert p["length"] == 18
    assert p["wraps_origin"] is True
    assert p["end"] == p["start"]  # 整圈回到起点
    assert len(p["sequence"]) == 18
    assert set(p["sequence"]) == set(seq)  # 整圈覆盖


def test_circular_binding_coordinates_are_modular(client):
    # fwd@18 在 L=24 上跨原点结合，end=(18+10)%24=4 < start
    chars = ["T"] * 24
    _place(chars, FWD, 18, 24)
    _place(chars, MOTIF, 6, 24)
    seq = "".join(chars)
    digest_id = _create(client, seq, topology="circular")
    products = _predict(client, digest_id, FWD, REV, None)["products"]
    binding = products[0]["forward_binding"]
    assert binding["start"] == 18
    assert binding["end"] == 4


def test_circular_linear_results_differ_for_same_sequence(client):
    # 同一 24 长模板：环形允许跨原点配对，产物严格更多
    chars = ["T"] * 24
    _place(chars, FWD, 18, 24)
    _place(chars, MOTIF, 6, 24)
    seq = "".join(chars)
    lin = _create(client, seq, topology="linear")
    circ = _create(client, seq, topology="circular")
    lin_products = _predict(client, lin, FWD, REV, None)["products"]
    circ_products = _predict(client, circ, FWD, REV, None)["products"]
    # 线性下 fwd(18..27) 越界无位点；环形有 1 个跨原点产物
    assert lin_products == []
    assert len(circ_products) == 1
    assert circ_products[0]["wraps_origin"] is True


# ---------- 去重与排序 ----------


def test_products_dedup_and_sort_by_start_then_length(client):
    seq = "ATGC" + "TTTT" + "GCAT"
    digest_id = _create(client, seq)
    products = _predict(client, digest_id, "NNNN", "NNNN", None)["products"]
    keys = [(p["start"], p["length"]) for p in products]
    assert keys == sorted(keys)
    assert len(keys) == len(set(keys))


# ---------- 状态码与校验 ----------


def test_missing_record_returns_404(client):
    resp = client.post(
        "/digests/999/pcr-products",
        json={"forward_primer": "ATGC", "reverse_primer": "GCAT", "max_product_length": 10},
    )
    assert resp.status_code == 404


def test_max_product_length_exceeding_template_returns_422(client):
    digest_id = _create(client, FWD + MOTIF)
    resp = client.post(
        f"/digests/{digest_id}/pcr-products",
        json={"forward_primer": "ATGC", "reverse_primer": "GCAT",
              "max_product_length": 19},
    )
    assert resp.status_code == 422
    assert "max_product_length 19 exceeds template length 18" in resp.text


def test_max_product_length_equal_template_allowed(client):
    digest_id = _create(client, FWD + MOTIF)
    data = _predict(client, digest_id, FWD, REV, 18)
    assert data["products"][0]["length"] == 18


def test_non_positive_max_length_returns_422(client):
    digest_id = _create(client, FWD + MOTIF)
    for value in (0, -1):
        resp = client.post(
            f"/digests/{digest_id}/pcr-products",
            json={"forward_primer": "ATGC", "reverse_primer": "TTAGCTGC",
                  "max_product_length": value},
        )
        assert resp.status_code == 422


@pytest.mark.parametrize("field", ["forward_primer", "reverse_primer"])
def test_short_or_empty_primer_returns_422(client, field):
    digest_id = _create(client, FWD + MOTIF)
    for primer in ("", "A", " A", "  "):
        payload = {"forward_primer": "ATGCGATCGC", "reverse_primer": "TTAGCTGC",
                   "max_product_length": 18}
        payload[field] = primer
        resp = client.post(f"/digests/{digest_id}/pcr-products", json=payload)
        assert resp.status_code == 422, primer


def test_non_iupac_primer_returns_422(client):
    digest_id = _create(client, FWD + MOTIF)
    for bad in ("ATGCGATCGU", "ACGT1ACGT", "ACG ACGT"):
        payload = {"forward_primer": bad, "reverse_primer": "TTAGCTGC",
                   "max_product_length": 18}
        resp = client.post(f"/digests/{digest_id}/pcr-products", json=payload)
        assert resp.status_code == 422, bad


def test_missing_required_fields_returns_422(client):
    digest_id = _create(client, FWD + MOTIF)
    resp = client.post(f"/digests/{digest_id}/pcr-products", json={})
    assert resp.status_code == 422


def test_lowercase_primers_normalized(client):
    digest_id = _create(client, FWD + MOTIF)
    data = _predict(client, digest_id, FWD.lower(), REV.lower(), 18)
    assert data["forward_primer"] == FWD
    assert data["reverse_primer"] == REV


def test_primer_longer_than_template_yields_empty_not_422(client):
    digest_id = _create(client, "ATGC")
    data = _predict(client, digest_id, "ATGCATGCAA", "ATGCATGCAA", 4)
    assert data["products"] == []


# ---------- 只读 ----------


def test_prediction_is_read_only(client):
    digest_id = _create(client, FWD + "T" * 10 + MOTIF)
    before = client.get("/digests").json()
    record_before = client.get(f"/digests/{digest_id}").json()
    _predict(client, digest_id, FWD, REV, 28)
    after = client.get("/digests").json()
    record_after = client.get(f"/digests/{digest_id}").json()
    assert after == before
    assert record_after == record_before


def test_prediction_ignores_digest_fragments(client):
    # 酶 EcoRI 在模板上没有任何位点：PCR 预测仍基于原始序列正常工作
    seq = FWD + "T" * 10 + MOTIF
    digest_id = _create(client, seq)
    record = client.get(f"/digests/{digest_id}").json()
    assert record["cuts"] == []
    assert len(record["fragments"]) == 1  # 未酶切的整条分子
    products = _predict(client, digest_id, FWD, REV, 28)["products"]
    assert len(products) == 1
