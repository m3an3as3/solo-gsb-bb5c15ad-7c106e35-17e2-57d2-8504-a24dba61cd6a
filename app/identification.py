"""未知样本的酶切片段长度识别（digest identification，只读）。

识别口径
--------
* 每条已保存记录是一个**候选**：取其完整酶切结果的片段长度（按片段索引
  顺序）作为**预期片段**，与请求给出的**观测片段长度**做二分图配对；
* 每个观测片段与每个预期片段至多参与一次配对，且仅当长度差绝对值
  ``|observed - expected| <= tolerance``（容差含端点）时才允许配对；
  重复出现的长度按不同片段计数（多重集语义，不做去重）；
* 在所有合法配对中：**先最大化配对数，再最小化配对差值总和**；仍并列时
  按配对列表 ``(观测索引, 片段索引)`` 的字典序取最小者——观测索引小的
  优先配对、同一观测优先配向片段索引小者，保证相同请求的结果唯一确定；
* 候选得分 =（未配对片段总数 = 观测侧与预期侧未配对数之和, 配对差值
  总和），候选按（未配对片段总数, 差值总和, 记录请求顺序）升序稳定排列；
* 未配对片段总数不超过 ``max_unmatched`` 的候选为**命中**；最优得分
  （排序后首个命中的得分）并列的命中多于一个时标记歧义；没有任何候选
  满足 ``max_unmatched`` 时无命中；
* 整个识别只在内存中进行，不读写数据库。

配对以观测顺序动态规划（状态 = 观测下标 × 预期片段占用掩码）求得最优
（配对数, 差值总和），再按观测索引升序贪心重建字典序最小的配对方案。
每侧至多 :data:`MAX_IDENTIFICATION_FRAGMENTS` 个片段，规模有界。
"""

from __future__ import annotations

# 每侧（观测 / 单条记录预期）允许的最大片段数：
# 状态空间 = 观测数 × 2^预期数，超出该上限无法在保证耗时的前提下
# 求精确最优配对，路由以 422 拒绝。
MAX_IDENTIFICATION_FRAGMENTS = 12


def optimal_pairing(observed, expected, tolerance):
    """求观测长度与预期长度之间的最优配对。

    ``observed`` / ``expected`` 为正整数长度列表（下标即对外返回的
    观测索引 / 片段索引）；``tolerance`` 为非负整数容差（含端点）。

    返回 ``(pairs, unmatched_observed, unmatched_expected)``：
    ``pairs`` 为 ``(观测索引, 片段索引, 差值)`` 列表，按观测索引升序；
    两个未配对列表分别为未配对的观测索引与片段索引，各自升序。
    """
    n, m = len(observed), len(expected)
    # 每个观测的兼容预期片段 (片段索引, 差值)，按片段索引升序
    compatible = []
    for i in range(n):
        row = []
        for j in range(m):
            diff = abs(observed[i] - expected[j])
            if diff <= tolerance:
                row.append((j, diff))
        compatible.append(row)

    # 得分编码为整数：配对数 * SCALE - 差值总和。每对差值 <= tolerance，
    # 至多 n 对，差值总和恒 < SCALE，因此整数比较等价于
    # （配对数, -差值总和）的字典序比较
    scale = n * tolerance + 1
    size = 1 << m
    # dp[i][mask] = 观测 i..n-1 在预期占用掩码 mask 下的最优得分；
    # 自底向上递推：跳过观测 i，或将其配给某个未占用且兼容的 j
    dp = [[0] * size for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        row, nxt = dp[i], dp[i + 1]
        compat_i = compatible[i]
        for mask in range(size):
            best = nxt[mask]  # 观测 i 不配对
            for j, diff in compat_i:
                if not mask >> j & 1:
                    candidate = nxt[mask | 1 << j] + scale - diff
                    if candidate > best:
                        best = candidate
            row[mask] = best

    # 贪心重建字典序最小的配对：观测索引升序，能配则配向可用的最小
    # 片段索引（该选择只要与全局最优得分一致即被固定），否则不配对
    pairs = []
    mask = 0
    for i in range(n):
        for j, diff in compatible[i]:
            if mask >> j & 1:
                continue
            if dp[i + 1][mask | 1 << j] + scale - diff == dp[i][mask]:
                pairs.append((i, j, diff))
                mask |= 1 << j
                break

    matched_observed = {i for i, _, _ in pairs}
    matched_expected = {j for _, j, _ in pairs}
    unmatched_observed = [i for i in range(n) if i not in matched_observed]
    unmatched_expected = [j for j in range(m) if j not in matched_expected]
    return pairs, unmatched_observed, unmatched_expected


def identify_candidates(records, observed, tolerance):
    """对每条记录计算最优配对与得分，返回按得分升序稳定排列的候选列表。

    ``records`` 为 ``{"digest_id": int, "fragment_lengths": [int, ...]}``
    序列，顺序即请求中 digest_ids 的顺序（同分时按该顺序稳定排列）。
    调用方须已完成全部校验（记录数 1~8、观测为 1~12 个正整数、容差
    非负、每条记录片段数不超过 :data:`MAX_IDENTIFICATION_FRAGMENTS`）。

    每个候选：{"digest_id", "fragment_lengths", "matched_count",
    "unmatched_count", "total_difference", "pairs",
    "unmatched_observed_indices", "unmatched_fragment_indices"}。
    """
    candidates = []
    for record in records:
        expected = record["fragment_lengths"]
        pairs, unmatched_observed, unmatched_expected = optimal_pairing(
            observed, expected, tolerance
        )
        candidates.append(
            {
                "digest_id": record["digest_id"],
                "fragment_lengths": list(expected),
                "matched_count": len(pairs),
                "unmatched_count": len(unmatched_observed)
                + len(unmatched_expected),
                "total_difference": sum(diff for _, _, diff in pairs),
                "pairs": [
                    {
                        "observed_index": i,
                        "fragment_index": j,
                        "observed_length": observed[i],
                        "expected_length": expected[j],
                        "difference": diff,
                    }
                    for i, j, diff in pairs
                ],
                "unmatched_observed_indices": unmatched_observed,
                "unmatched_fragment_indices": unmatched_expected,
            }
        )
    # 稳定排序：未配对片段总数升序、差值总和升序，再按记录请求顺序
    candidates.sort(key=lambda c: (c["unmatched_count"], c["total_difference"]))
    return candidates
