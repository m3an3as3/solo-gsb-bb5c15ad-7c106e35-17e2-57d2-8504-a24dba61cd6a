"""开放阅读框（ORF）检索（只读，不依赖酶切片段）。

输入
----
已保存记录的原始 IUPAC 序列与线性/环形拓扑，以及最短氨基酸数
``min_aa_length``（正整数）。

检索口径
--------
* 在**正链**（top 链 5'->3'）与**反向互补链**上分别**逐碱基**寻找
  **字面量** ``ATG``：简并符号（N/R/...）组成的密码子不充当起始信号。
* 从每个起点按三联体逐个密码子读取，遇到**同框**的第一个**字面量**
  ``TAA`` / ``TAG`` / ``TGA`` 才形成 ORF：不同框的终止密码子忽略；
  含简并符号的密码子既不是起始也不是终止，读取继续。
* **线性**：读取不得越出模板边界，边界前没有同框终止密码子的起点不
  产生 ORF。
* **环形**：读取可跨原点（在复制串上扫描），但**完整 ORF（含终止
  密码子）不得超过模板一圈**；读取超过一圈仍未见终止密码子的起点不
  产生 ORF。
* 保留满足 ``(密码子数 - 1) >= min_aa_length``（氨基酸数不含终止）的
  **全部嵌套 ORF**：每个 ATG 各自独立寻找其同框第一个终止密码子，
  同框内部 ATG 给出的更短嵌套 ORF 全部保留。

坐标口径
--------
所有返回坐标均为 top 链（原始正链）上的 0-based 半开区间端点，与
酶切片段 / PCR 产物的坐标约定一致：

* 正链 ORF：``start`` 为 ATG 首碱基坐标，``end`` 为终止密码子之后
  的位置；环形跨原点时 ``end`` 为取模后的坐标（可能小于 ``start``）。
* 反链 ORF：``start`` / ``end`` 同样给出其在 top 链上覆盖区间的
  左、右端点（与 ORF 自身 5'->3' 方向相反）；``sequence`` 始终为
  ORF 自身 5'->3' 的定向核苷酸序列（即反向互补链上的序列）。
* ``wraps_origin`` 为 true 表示该链 5'->3' 定向读取跨过模板原点。
"""

from __future__ import annotations

from .iupac import reverse_complement

START_CODON = "ATG"
STOP_CODONS = frozenset(("TAA", "TAG", "TGA"))


def _next_inframe_stops(oriented: str, circular: bool) -> dict[int, int]:
    """返回 ``{同框密码子起点: 该框右侧第一个字面量终止密码子的原始末端}``。

    键、值均取在扫描串（线性为原序列、环形为复制串）上的**原始下标**：
    值为严格位于该位置之后、同一阅读框（下标 mod 3 相同）的第一个字面量
    TAA/TAG/TGA 的结束位置（半开，即末碱基下标 + 1），不存在则不入表。
    三个阅读框分别自右向左扫描一遍，整体 O(L)。
    """
    length = len(oriented)
    haystack = oriented + oriented if circular else oriented
    limit = len(haystack) - 3
    next_stop: dict[int, int] = {}
    for frame in (0, 1, 2):
        pending: int | None = None
        for position in range(limit - ((limit - frame) % 3), frame - 1, -3):
            codon = haystack[position : position + 3]
            if codon in STOP_CODONS:
                pending = position + 3
            if pending is not None:
                next_stop[position] = pending
    return next_stop


def _scan_strand(oriented: str, circular: bool) -> list[tuple[int, int, int]]:
    """在一条 5'->3' 定向序列上枚举全部 (起点, 长度, 氨基酸数)。

    长度含终止密码子；氨基酸数不含终止。逐碱基寻找字面量 ATG，每个
    起点独立取其同框第一个字面量终止密码子（由 ``_next_inframe_stops``
    预先索引）。

    环形在 ``oriented + oriented`` 复制串上扫描：起点范围仍为
    ``[0, L)``，完整 ORF 的定向原始末端必须 ``<= start + L``
    （自起点起至多一圈）。
    """
    length = len(oriented)
    haystack = oriented + oriented if circular else oriented
    next_stop = _next_inframe_stops(oriented, circular)
    orfs: list[tuple[int, int, int]] = []
    for start in range(length):
        # 环形起点 L-2/L-1 的 ATG 跨原点，须在复制串上判定
        if haystack[start : start + 3] != START_CODON:
            continue
        end_raw = next_stop.get(start)
        if end_raw is None:
            continue  # 线性越界前无同框终止，或环形超过一圈仍无终止
        if circular and end_raw > start + length:
            continue  # 完整 ORF（含终止密码子）超过模板一圈
        orf_length = end_raw - start
        orfs.append((start, orf_length, orf_length // 3 - 1))
    return orfs


def find_orfs(sequence: str, topology: str, min_aa_length: int) -> list[dict]:
    """检索满足最短氨基酸数的全部 ORF，稳定排序后返回。

    调用方负责输入校验（序列为非空 IUPAC、拓扑合法、
    ``min_aa_length >= 1``）。
    """
    length = len(sequence)
    circular = topology == "circular"

    orfs: list[dict] = []

    for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
        # 环形跨原点 ORF 直接在复制串上切片：start < L 且长度 <= L，
        # 故 end_raw <= start + L < 2L
        haystack = oriented + oriented if circular else oriented
        for start, orf_length, aa_length in _scan_strand(oriented, circular):
            if aa_length < min_aa_length:
                continue
            end_raw = start + orf_length
            # 定向读取跨原点的判定与 PCR 产物/GC 窗口一致：定向原始末端
            # 严格超过 L 才算跨原点；恰好终止于原点边界不算
            wraps = circular and end_raw > length
            # ORF 自身 5'->3' 的定向核苷酸序列（环形跨原点时拼接环首尾）
            orf_sequence = haystack[start:end_raw]
            if strand == "+":
                top_start = start
                top_end = end_raw % length if circular else end_raw
            else:
                # 定向序列上 [start, end_raw) 对应 top 链上
                # [L - end_raw, L - start)；ORF 的 5'->3' 序列即定向序列
                top_start = (
                    (length - end_raw) % length if circular else length - end_raw
                )
                top_end = (
                    (length - start) % length if circular else length - start
                )
            orfs.append(
                {
                    "strand": strand,
                    "start": top_start,
                    "end": top_end,
                    "wraps_origin": wraps,
                    "sequence": orf_sequence,
                    "length": orf_length,
                    "aa_length": aa_length,
                }
            )

    # 按链方向（+ 在 - 前）、原始正链起点、长度稳定排序
    orfs.sort(key=lambda o: (o["strand"], o["start"], o["length"]))
    return orfs
