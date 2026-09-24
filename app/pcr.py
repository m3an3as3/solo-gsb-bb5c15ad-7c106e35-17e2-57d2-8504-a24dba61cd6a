"""PCR 扩增产物预测（只读，不依赖酶切片段）。

输入
----
已保存记录的原始序列与线性/环形拓扑，以及两条 IUPAC 引物（均按 5'->3' 书写）：

* **正引物** ``forward_primer``：结合 top 链，沿模板正向（向右）延伸；
* **反引物** ``reverse_primer``：结合 bottom 链，以其**反向互补序列**
  ``reverse_complement(reverse_primer)`` 在 top 链上匹配，沿模板反向（向左）延伸。

匹配口径
--------
* 引物碱基与模板碱基的 IUPAC 可能集合**有交集**即视为该位置匹配；
* 每条引物的整个结合区允许**至多一处错配**；
* 但各自 **3' 端最后两位必须匹配**（该两位任意一处不兼容即整条结合不成立）。
  正引物与 top 链同向，3' 端位于匹配串末尾；反引物先取反向互补再匹配，
  其 3' 端位于反向互补串的开头两位。

产物口径
--------
枚举正、反结合位点中**相向**（正向位点在前、反向位点在后）且**结合区互不重叠**
的位点对：

* **线性**：正向位点 ``[f, f+kf)`` 完全位于反向位点 ``[r, r+kr)`` 之前
  （``f + kf <= r``），产物区间 ``[f, r+kr)`` 不得越出 ``[0, L)``；
  产物长度为 ``r + kr - f``。
* **环形**：以 ``d = (r - f) mod L`` 为正向位点起点沿环向前到反向位点起点的
  有向距离，要求 ``d >= kf``（相向且结合区不重叠）；产物长度为 ``d + kr``，
  允许跨原点（``f + 长度 > L`` 即 ``wraps_origin``），但**不得超过一圈**
  （长度恒 ``<= L``）。

长度超过请求给定 ``max_product_length`` 的产物不返回。产物按**起点、长度**
去重并升序排列；产物序列取 top 链 5'->3' 方向（跨原点时拼接环首尾）。
"""

from __future__ import annotations

from .iupac import IUPAC_SETS, reverse_complement

# 每条引物整个结合区允许的最大错配数
MAX_MISMATCHES = 1

# 引物 3' 端必须精确匹配（IUPAC 集合有交集）的末尾碱基数
TERMINAL_MATCH_LENGTH = 2


def _match_sites(sequence: str, motif: str, circular: bool, terminal_at_end: bool):
    """枚举 motif 在 top 链上的全部结合起点与错配数。

    IUPAC 集合有交集即匹配；除 3' 端最后两位外，全结合区允许至多
    ``MAX_MISMATCHES`` 处错配。

    引物 5'->3' 书写：正引物与 top 链同向，其 3' 端对应 motif 的**末尾**
    两位（``terminal_at_end=True``）；反引物先取反向互补再与 top 链匹配，
    其 3' 端位于 motif 的**开头**两位（``terminal_at_end=False``）。

    环形序列在 ``seq + seq[:k-1]`` 上扫描，起点范围仍为 ``[0, L)``；
    motif 长于序列时无任何位点。生成 ``(start, mismatch_count)``。
    """
    length = len(sequence)
    k = len(motif)
    if k == 0 or k > length:
        return
    haystack = sequence + sequence[: k - 1] if circular else sequence
    last_start = length if circular else length - k + 1
    # 3' 端最后两位在 motif 中的下标：正引物在末尾，反引物在开头
    terminal = {k - TERMINAL_MATCH_LENGTH, k - 1} if terminal_at_end else {0, 1}
    for start in range(last_start):
        mismatches = 0
        for j in range(k):
            if not (IUPAC_SETS[motif[j]] & IUPAC_SETS[haystack[start + j]]):
                if j in terminal:
                    # 3' 端最后两位不得错配：该位点直接不成立
                    mismatches = MAX_MISMATCHES + 1
                    break
                mismatches += 1
                if mismatches > MAX_MISMATCHES:
                    break
        else:
            yield start, mismatches


def predict_pcr_products(
    sequence: str,
    topology: str,
    forward_primer: str,
    reverse_primer: str,
    max_product_length: int,
) -> list[dict]:
    """预测全部满足约束的 PCR 产物，按 (起点, 长度) 升序去重返回。

    调用方负责输入校验（引物为非空 IUPAC 序列、长度至少为 2、
    ``1 <= max_product_length <= len(sequence)``）。
    """
    length = len(sequence)
    circular = topology == "circular"
    kf = len(forward_primer)
    kr = len(reverse_primer)

    # 每个起点至多一个结合（起点即位点身份）；dict 便于后续按起点枚举
    forward_sites = dict(
        _match_sites(sequence, forward_primer, circular, terminal_at_end=True)
    )
    reverse_motif = reverse_complement(reverse_primer)
    reverse_sites = dict(
        _match_sites(sequence, reverse_motif, circular, terminal_at_end=False)
    )

    # 环形跨原点产物直接在复制串上切片：f < L、长度 <= L，故 f+长度 <= 2L-1
    doubled = sequence + sequence if circular else None

    # 不同位点对（尤其简并引物）可能给出相同的 (起点, 长度)，按此键去重
    products: dict[tuple[int, int], dict] = {}
    for f, forward_mismatches in forward_sites.items():
        for r, reverse_mismatches in reverse_sites.items():
            if circular:
                distance = (r - f) % length
                if distance < kf:
                    # 反向位点落在正向位点结合区内或其后不足一个结合区：
                    # 二者不相向或结合区重叠
                    continue
                product_length = distance + kr
                if product_length > length:
                    continue  # 环形产物不得超过一圈
                end_raw = f + product_length
                if product_length > max_product_length:
                    continue
                products[(f, product_length)] = {
                    "sequence": doubled[f:end_raw],
                    "start": f,
                    "end": end_raw % length,
                    "length": product_length,
                    "wraps_origin": end_raw > length,
                    "forward_binding": {
                        "start": f,
                        "end": (f + kf) % length,
                        "mismatches": forward_mismatches,
                    },
                    "reverse_binding": {
                        "start": r,
                        "end": (r + kr) % length,
                        "mismatches": reverse_mismatches,
                    },
                }
            else:
                if r < f + kf:
                    continue  # 不相向或结合区重叠
                end_raw = r + kr
                product_length = end_raw - f
                if end_raw > length:
                    continue  # 线性产物不得越界
                if product_length > max_product_length:
                    continue
                products[(f, product_length)] = {
                    "sequence": sequence[f:end_raw],
                    "start": f,
                    "end": end_raw,
                    "length": product_length,
                    "wraps_origin": False,
                    "forward_binding": {
                        "start": f,
                        "end": f + kf,
                        "mismatches": forward_mismatches,
                    },
                    "reverse_binding": {
                        "start": r,
                        "end": r + kr,
                        "mismatches": reverse_mismatches,
                    },
                }

    return [products[key] for key in sorted(products)]
