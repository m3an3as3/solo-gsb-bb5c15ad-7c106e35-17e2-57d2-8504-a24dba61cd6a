"""开放阅读框（ORF）检索接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.iupac import reverse_complement
from app.main import app
from app.orf import find_orfs


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


def _search(client, digest_id, min_aa_length, expect=200):
    resp = client.post(
        f"/digests/{digest_id}/orfs",
        json={"min_aa_length": min_aa_length},
    )
    assert resp.status_code == expect, resp.text
    return resp.json() if expect == 200 else resp


# ---------- 正链线性 ----------


def test_basic_forward_orf(client):
    # ATG AAA TTT TAG：3 个氨基酸（不含终止）
    digest_id = _create(client, "ATGAAATTTTAG")
    data = _search(client, digest_id, 1)
    assert data["id"] == digest_id
    assert data["topology"] == "linear"
    assert data["sequence_length"] == 12
    assert data["min_aa_length"] == 1
    assert data["orf_count"] == 1
    assert data["orfs"] == [
        {
            "strand": "+",
            "start": 0,
            "end": 12,
            "wraps_origin": False,
            "sequence": "ATGAAATTTTAG",
            "length": 12,
            "aa_length": 3,
        }
    ]


def test_orf_sequence_starts_with_atg_and_ends_with_literal_stop(client):
    digest_id = _create(client, "ATGAAATAG")
    orf = _search(client, digest_id, 1)["orfs"][0]
    assert orf["sequence"][:3] == "ATG"
    assert orf["sequence"][-3:] == "TAG"
    assert orf["length"] == 3 * (orf["aa_length"] + 1)


def test_all_three_stop_codons_close_orfs(client):
    # ATG + 终止密码子：恰为 1 个氨基酸
    for stop in ("TAA", "TAG", "TGA"):
        digest_id = _create(client, "ATG" + stop)
        orfs = _search(client, digest_id, 1)["orfs"]
        assert len(orfs) == 1, stop
        assert orfs[0]["sequence"][-3:] == stop
        assert orfs[0]["aa_length"] == 1
        assert orfs[0]["length"] == 6


def test_nested_orfs_all_retained(client):
    # ATG ATG AAT TAG：起点 0 给 3 aa，起点 3（同框嵌套 ATG）给 2 aa，
    # 两个都保留
    digest_id = _create(client, "ATGATGAATTAG")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert [(o["start"], o["aa_length"]) for o in orfs] == [(0, 3), (3, 2)]
    assert [o["sequence"] for o in orfs] == ["ATGATGAATTAG", "ATGAATTAG"]


def test_min_aa_length_filters_inclusive(client):
    # ATG AAT TAG：恰为 2 个氨基酸
    digest_id = _create(client, "ATGAATTAG")
    assert _search(client, digest_id, 2)["orfs"]
    assert _search(client, digest_id, 3)["orfs"] == []
    assert _search(client, digest_id, 3)["orf_count"] == 0


def test_every_base_position_checked_for_atg_not_just_three_frames(client):
    # 从位置 2 开始的同框 ORF：逐碱基扫描才能找到（并非按阅读框分桶）
    digest_id = _create(client, "TTATGAATTAGC")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    assert (orfs[0]["start"], orfs[0]["end"], orfs[0]["aa_length"]) == (2, 11, 2)


def test_out_of_frame_stop_is_ignored(client):
    # 位置 1..3 的 TGA 与起点 0 错框，不得作为终止；同框终止是位置 12
    # 的 TAA：ATG ACT GGT GGT TAA 共 4 个氨基酸
    seq = "ATGACTGGTGGTTAA"
    assert seq[1:4] == "TGA"  # 错框终止不得截断
    digest_id = _create(client, seq)
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    assert (orfs[0]["start"], orfs[0]["end"], orfs[0]["aa_length"]) == (0, 15, 4)
    assert orfs[0]["sequence"] == seq


def test_linear_orf_cannot_extend_past_boundary(client):
    # 边界前没有同框终止密码子：不产生 ORF
    digest_id = _create(client, "ATGAAATTTCCC")
    assert _search(client, digest_id, 1)["orfs"] == []


def test_truncated_triplet_at_boundary_is_not_a_stop(client):
    # 末端只剩一个/两个碱基：不完整密码子不是终止，无 ORF
    digest_id = _create(client, "ATGAAATAA")
    assert _search(client, digest_id, 1)["orfs"]


def test_no_atg_returns_empty_array(client):
    digest_id = _create(client, "CCCGGGAAATTT")
    data = _search(client, digest_id, 1)
    assert data["orfs"] == []
    assert data["orf_count"] == 0


# ---------- 反向互补链 ----------


def test_reverse_strand_orf(client):
    # top = rc("ATGAAATAG")：反向互补链上存在 ATG AAA TAG
    seq = reverse_complement("ATGAAATAG")
    digest_id = _create(client, seq)
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["strand"] == "-"
    assert orf["sequence"] == "ATGAAATAG"  # ORF 自身 5'->3' 定向序列
    assert orf["start"] == 0
    assert orf["end"] == 9
    assert orf["wraps_origin"] is False
    assert orf["aa_length"] == 2


def test_reverse_strand_coordinates_reference_plus_strand(client):
    # 反链 ORF 起始（ATG）位于 top 链右端、终止位于左端
    # rc(ATGAATTAG) = CTAATTCAT：ATG 对应 top 位置 6..8，TAG 对应 0..2
    seq = reverse_complement("ATGAATTAG")
    digest_id = _create(client, seq)
    orf = _search(client, digest_id, 1)["orfs"][0]
    assert orf["strand"] == "-"
    assert orf["start"] == 0
    assert orf["end"] == 9
    assert seq[6:9] == "CAT"  # 反链 5' 端 ATG 在 top 上为 CAT


def test_both_strands_scanned_and_sorted_plus_first(client):
    # 正链：ATG AAA TAA（9，2 aa）；反链 ORF 落在另一处
    # top = "ATGAAATAA" 同时其 rc = TTATTTTCAT 无 ATG，故另构造：
    # top 前 9 位正链 ORF；后接 rc(ATGTAA)="TTACAT" 提供反链 ATG TAA（1 aa）
    seq = "ATGAAATAA" + reverse_complement("ATGTAA")
    digest_id = _create(client, seq)
    orfs = _search(client, digest_id, 1)["orfs"]
    strands = [o["strand"] for o in orfs]
    assert strands == sorted(strands)
    assert orfs[0]["strand"] == "+"
    assert any(o["strand"] == "-" and o["sequence"] == "ATGTAA" for o in orfs)


# ---------- 简并符号 ----------


def test_degenerate_codon_is_not_start(client):
    # 位置 0 的 NTA 不是字面量 ATG；真正的 ATG 在位置 2
    digest_id = _create(client, "NTATGAAATAG")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    assert orfs[0]["start"] == 2
    assert orfs[0]["aa_length"] == 2


def test_degenerate_codon_is_not_stop(client):
    # TAR（可能为 TAA/TAG）是简并密码子：不充当终止，读取继续到位置 9 的 TAG
    digest_id = _create(client, "ATGAAATARTAG")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    assert orfs[0]["aa_length"] == 3
    assert orfs[0]["sequence"] == "ATGAAATARTAG"


def test_degenerate_codon_inside_orf_does_not_break_frame(client):
    # 简并密码子按普通三联体处理：ATG AAN TTT TAA -> 3 aa
    digest_id = _create(client, "ATGAANTTTTAA")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    assert (orfs[0]["aa_length"], orfs[0]["length"]) == (3, 12)


# ---------- 环形 ----------


def test_circular_orf_wraps_origin(client):
    # AAA TAG ATG（L=9）：位置 6 的 ATG 跨原点先读位置 0..2 的 AAA（1 aa），
    # 同框终止 TAG 位于位置 3..5；完整 ORF 长 9（= 一圈），2 个氨基酸
    seq = "AAATAGATG"
    digest_id = _create(client, seq, topology="circular")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["strand"] == "+"
    assert orf["start"] == 6
    assert orf["end"] == 6
    assert orf["wraps_origin"] is True
    assert orf["sequence"] == "ATGAAATAG"
    assert orf["length"] == 9
    assert orf["aa_length"] == 2


def test_circular_atg_split_across_origin(client):
    # L=10：位置 9,0,1 拼出 ATG；同框终止 TAA 在位置 2..4
    chars = ["C"] * 10
    for i, base in enumerate("ATG"):
        chars[(9 + i) % 10] = base
    for i, base in enumerate("TAA"):
        chars[2 + i] = base
    seq = "".join(chars)
    assert seq == "TGTAACCCCA"
    digest_id = _create(client, seq, topology="circular")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["start"] == 9
    assert orf["end"] == 5
    assert orf["wraps_origin"] is True
    assert orf["sequence"] == "ATGTAA"
    assert orf["aa_length"] == 1


def test_circular_full_turn_allowed(client):
    # L=12：ATG AAA AAA TAA 恰好完整一圈，wraps_origin 为 False
    seq = "ATGAAAAAATAA"
    digest_id = _create(client, seq, topology="circular")
    orf = _search(client, digest_id, 1)["orfs"][0]
    assert orf["length"] == 12
    assert orf["start"] == orf["end"] == 0
    assert orf["wraps_origin"] is False


def test_circular_orf_must_close_within_one_turn(client):
    # L=13：起点 0 的同框终止密码子第一次出现在第二圈（超出一圈）-> 无 ORF
    seq = "ATGGGTAACGGGG"
    digest_id = _create(client, seq, topology="circular")
    assert _search(client, digest_id, 1)["orfs"] == []
    # 线性下同样越界无 ORF
    lin_id = _create(client, seq, topology="linear")
    assert _search(client, lin_id, 1)["orfs"] == []


def test_circular_wrapping_orf_not_found_on_linear(client):
    seq = "AAATAGATG"
    circ = _create(client, seq, topology="circular")
    lin = _create(client, seq, topology="linear")
    circ_orfs = _search(client, circ, 1)["orfs"]
    lin_orfs = _search(client, lin, 1)["orfs"]
    assert len(circ_orfs) == 1 and circ_orfs[0]["wraps_origin"] is True
    assert lin_orfs == []


def test_reverse_strand_circular_orf_wraps_origin(client):
    # 反链定向序列 GTT ATA GTA T（L=10）：位置 8,9,0 拼出 ATG，
    # 同框终止 TAT... 实际同框终止为位置 2..4 的 TAT 非终止；
    # 构造：ATG + TAT + TAG 跨原点，1 个间隔密码子 -> 2 aa
    oriented = "GTTATAGTAT"
    seq = reverse_complement(oriented)  # ATACTATAAC
    digest_id = _create(client, seq, topology="circular")
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["strand"] == "-"
    assert orf["start"] == 3
    assert orf["end"] == 2
    assert orf["wraps_origin"] is True
    assert orf["sequence"] == "ATGTTATAG"
    assert orf["length"] == 9
    assert orf["aa_length"] == 2
    # 同一序列线性下跨原点 ORF 消失
    lin_id = _create(client, seq, topology="linear")
    assert _search(client, lin_id, 1)["orfs"] == []


def test_circular_nested_orfs_wrapping(client):
    # L=12：ATG TAA CCC ATG。起点 0 给 1 aa（不跨原点）；起点 9 的 ATG
    # 同框跨原点后读到位置 3 的 TAA（中间还经过位置 0 的同框嵌套 ATG），
    # 形成 2 aa 的跨原点 ORF
    seq = "ATGTAACCCATG"
    digest_id = _create(client, seq, topology="circular")
    orfs = _search(client, digest_id, 1)["orfs"]
    keys = sorted((o["start"], o["aa_length"], o["wraps_origin"]) for o in orfs)
    assert keys == [(0, 1, False), (9, 2, True)]
    wrap = next(o for o in orfs if o["wraps_origin"])
    assert wrap["sequence"] == "ATGATGTAA"
    assert wrap["end"] == 6


# ---------- 排序 ----------


def test_sorted_by_strand_start_length(client):
    digest_id = _create(client, "ATGATGAATTAG")
    orfs = _search(client, digest_id, 1)["orfs"]
    keys = [(o["strand"], o["start"], o["length"]) for o in orfs]
    assert keys == sorted(keys)


def test_sort_stable_across_calls(client):
    digest_id = _create(client, "ATGATGAATTAG" + reverse_complement("ATGAATTAG"))
    first = _search(client, digest_id, 1)
    second = _search(client, digest_id, 1)
    assert first["orfs"] == second["orfs"]


# ---------- 状态码与校验 ----------


def test_missing_record_returns_404(client):
    resp = client.post("/digests/999/orfs", json={"min_aa_length": 20})
    assert resp.status_code == 404


def test_non_positive_min_aa_length_returns_422(client):
    digest_id = _create(client, "ATGAAATAG")
    for value in (0, -1):
        resp = client.post(
            f"/digests/{digest_id}/orfs", json={"min_aa_length": value}
        )
        assert resp.status_code == 422, value


def test_missing_required_field_returns_422(client):
    digest_id = _create(client, "ATGAAATAG")
    resp = client.post(f"/digests/{digest_id}/orfs", json={})
    assert resp.status_code == 422


def test_non_integer_min_aa_length_returns_422(client):
    digest_id = _create(client, "ATGAAATAG")
    resp = client.post(
        f"/digests/{digest_id}/orfs", json={"min_aa_length": 1.5}
    )
    assert resp.status_code == 422


# ---------- 只读 / 不依赖酶切片段 ----------


def test_search_is_read_only(client):
    digest_id = _create(client, "ATGAAATAG")
    before = client.get("/digests").json()
    record_before = client.get(f"/digests/{digest_id}").json()
    _search(client, digest_id, 1)
    after = client.get("/digests").json()
    record_after = client.get(f"/digests/{digest_id}").json()
    assert after == before
    assert record_after == record_before


def test_search_ignores_digest_fragments(client):
    # 酶 EcoRI 在模板上无任何位点：ORF 检索仍基于原始序列正常工作
    seq = "ATGAAATAG"
    digest_id = _create(client, seq)
    record = client.get(f"/digests/{digest_id}").json()
    assert record["cuts"] == []
    assert len(record["fragments"]) == 1
    orfs = _search(client, digest_id, 1)["orfs"]
    assert len(orfs) == 1


# ---------- 算法级直接用例 ----------


def test_empty_result_fields_directly():
    assert find_orfs("CCCCGGGG", "linear", 1) == []


def test_short_sequence_has_no_room_for_orf():
    # 不足 6 个碱基（至少一个编码密码子 + 终止）不可能形成 ORF
    assert find_orfs("ATGTAA", "linear", 2) == []
    assert find_orfs("ATGTAA", "linear", 1)[0]["aa_length"] == 1
