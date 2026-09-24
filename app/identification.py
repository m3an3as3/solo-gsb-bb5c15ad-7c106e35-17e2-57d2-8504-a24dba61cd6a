"""未知样本识别：观测片段长度与已保存酶切片段的最优二分图配对（只读）。

识别口径
--------
* 未知样本给出一份**观测片段长度列表**（凝胶上读得的条带长度，重复长度
  按不同片段计数），每条已保存记录的酶切片段长度作为一份**预期清单**，
  两侧都以**列表索引**标识片段；
* 一个观测片段与一个预期片段**至多配对一次**（任一侧都不能被两个配对
  复用），配对当且仅当长度差绝对值不超过整数容差 ``tolerance``
  （闭区间，``|observed - expected| <= tolerance``）；
* 对每条候选记录求唯一的**最优配对方案**，优化目标依次为：
  1. **最大化配对数**；
  2. 配对数相同则**最小化差值总和**（``sum |observed - expected|``）；
  3. 仍相同时把配对按（观测索引, 预期索引）排序得到配对序列，取
     **字典序最小**者——更小编号的观测片段优先被配出（未配对因此集中在
     编号最大的观测上），必须配对时优先占用预期索引更小的片段；
* 重复长度不做折叠：长度相同的两个片段是不同索引，可分别配对；
* 未配对片段总数 = 观测未配对数 + 预期未配对数；候选按
  （未配对总数, 差值总和, 请求中的 id 顺序）排序。未配对总数不超过
  ``max_unpaired`` 的候选才算命中；最优（命中候选中分数最高）分数
  并列时标记歧义。

算法
----
带记忆化的序列 DP：状态 ``(i, used_mask)`` 表示依次决策观测片段
``i..n-1``、``used_mask`` 标记已占用的预期索引，返回
``(最大配对数, 最小差值总和)``，每个观测可跳过或与任一未占用且容差内的
预期片段配对。状态数有界（每侧至多 :data:`MAX_IDENTIFICATION_FRAGMENTS`
个片段，``(n+1) * 2^m`` 有界）。

字典序择优用贪心完成：按观测索引升序，在每个索引上按预期索引升序探测
第一个可扩展为全局最优（配对数与差值总和均达到后缀最优）的配对；不存在
这样的配对时该观测保持未配对。于是更小编号的观测优先被配出、未配对集中
在编号最大的观测上，被配出的观测占用编号最小的可行预期。可行性判定全部
基于后缀最优值，结果与枚举顺序无关，相同请求必然得到相同方案。

整个识别只在内存中进行，不读写数据库。
"""

from __future__ import annotations

from functools import lru_cache

# 每侧（观测列表 / 单条已保存记录）参与识别的片段数上限：
# 记忆化状态数 = (n+1) * 2^m，超出该上限无法在保证耗时的前提下完成全枚举，
# 观测侧由请求模型校验、预期侧由路由在取到记录后校验，均以 422 拒绝。
MAX_IDENTIFICATION_FRAGMENTS = 12


def best_assignment(
    observed_lengths: list[int],
    expected_lengths: list[int],
    tolerance: int,
) -> dict:
    """求一份观测清单与一份预期清单之间的最优配对。

    返回 ``{"pairs", "unpaired_observed_indices", "unpaired_expected_indices",
    "pair_count", "diff_sum"}``；``pairs`` 为
    ``{"observed_index", "expected_index", "observed_length",
    "expected_length", "length_diff"}`` 列表，按观测索引升序；
    未配对索引列表按索引升序。调用方须已完成全部校验（每侧 1~12 个
    片段、长度为正整数、tolerance 为非负整数）。
    """
    n = len(observed_lengths)
    m = len(expected_lengths)

    @lru_cache(maxsize=None)
    def solve(idx: int, used_mask: int) -> tuple[int, int]:
        """决策观测 idx..n-1（used_mask 标记已占用预期）的（最大配对数, 最小差值和）。"""
        if idx == n:
            return 0, 0
        # 观测 idx 不配对
        best_pairs, best_cost = solve(idx + 1, used_mask)
        observed = observed_lengths[idx]
        for j in range(m):
            bit = 1 << j
            if used_mask & bit:
                continue
            diff = abs(observed - expected_lengths[j])
            if diff > tolerance:
                continue
            sub_pairs, sub_cost = solve(idx + 1, used_mask | bit)
            pairs = sub_pairs + 1
            cost = sub_cost + diff
            # 先最大化配对数，再最小化差值总和
            if pairs > best_pairs or (pairs == best_pairs and cost < best_cost):
                best_pairs, best_cost = pairs, cost
        return best_pairs, best_cost

    total_pairs, total_cost = solve(0, 0)

    # 字典序贪心：观测索引升序；每个观测按预期索引升序取第一个仍可扩展为
    # 全局最优（配对数、差值总和双达标）的配对；不存在则该观测未配对。
    pairs: list[dict] = []
    unpaired_observed: list[int] = []
    used_mask = 0
    fixed_pairs = 0
    fixed_cost = 0
    for i in range(n):
        matched = False
        for j in range(m):
            bit = 1 << j
            if used_mask & bit:
                continue
            diff = abs(observed_lengths[i] - expected_lengths[j])
            if diff > tolerance:
                continue
            sub_pairs, sub_cost = solve(i + 1, used_mask | bit)
            if (
                fixed_pairs + 1 + sub_pairs == total_pairs
                and fixed_cost + diff + sub_cost == total_cost
            ):
                pairs.append(
                    {
                        "observed_index": i,
                        "expected_index": j,
                        "observed_length": observed_lengths[i],
                        "expected_length": expected_lengths[j],
                        "length_diff": diff,
                    }
                )
                used_mask |= bit
                fixed_pairs += 1
                fixed_cost += diff
                matched = True
                break
        if not matched:
            unpaired_observed.append(i)

    unpaired_expected = [j for j in range(m) if not (used_mask & (1 << j))]
    return {
        "pairs": pairs,
        "unpaired_observed_indices": unpaired_observed,
        "unpaired_expected_indices": unpaired_expected,
        "pair_count": total_pairs,
        "diff_sum": total_cost,
    }


def identify_sample(
    candidates: list[dict],
    observed_lengths: list[int],
    tolerance: int,
    max_unpaired: int,
) -> dict:
    """对若干候选记录逐一最优配对，排序并判定命中/歧义。

    ``candidates`` 为 ``{"digest_id": int, "fragment_lengths": [int, ...]}``
    序列，顺序即请求中 digest_ids 的顺序。调用方须已完成全部校验
    （候选数 1~8、每侧片段数 1~12、记录均存在）。
    """
    n = len(observed_lengths)
    results: list[dict] = []
    for position, candidate in enumerate(candidates):
        expected_lengths = candidate["fragment_lengths"]
        assignment = best_assignment(observed_lengths, expected_lengths, tolerance)
        results.append(
            {
                "digest_id": candidate["digest_id"],
                "position": position,
                "fragment_count": len(expected_lengths),
                "pair_count": assignment["pair_count"],
                "unpaired_count": n
                + len(expected_lengths)
                - 2 * assignment["pair_count"],
                "diff_sum": assignment["diff_sum"],
                "pairs": assignment["pairs"],
                "unpaired_observed_indices": assignment[
                    "unpaired_observed_indices"
                ],
                "unpaired_expected_indices": assignment[
                    "unpaired_expected_indices"
                ],
            }
        )

    # 命中门槛：未配对总数不超过 max_unpaired
    eligible = [r for r in results if r["unpaired_count"] <= max_unpaired]
    best_positions: set[int] = set()
    ambiguous = False
    best_digest_id: int | None = None
    if eligible:
        # 最优分数 = （未配对总数, 差值总和）最小；请求位置仅决定并列时的
        # 代表 id，不影响"并列"本身——同分数候选超过一个即歧义。
        best_unpaired = min(r["unpaired_count"] for r in eligible)
        best_diff = min(
            r["diff_sum"] for r in eligible if r["unpaired_count"] == best_unpaired
        )
        best_positions = {
            r["position"]
            for r in eligible
            if r["unpaired_count"] == best_unpaired and r["diff_sum"] == best_diff
        }
        ambiguous = len(best_positions) > 1
        representative = min(best_positions)
        best_digest_id = results[representative]["digest_id"]

    # 排序：未配对总数、差值总和、请求 ID 顺序
    ordered = sorted(
        results, key=lambda r: (r["unpaired_count"], r["diff_sum"], r["position"])
    )
    output_candidates = []
    for r in ordered:
        output_candidates.append(
            {
                "digest_id": r["digest_id"],
                "fragment_count": r["fragment_count"],
                "pair_count": r["pair_count"],
                "unpaired_count": r["unpaired_count"],
                "diff_sum": r["diff_sum"],
                "eligible": r["unpaired_count"] <= max_unpaired,
                "best": r["position"] in best_positions,
                "pairs": r["pairs"],
                "unpaired_observed_indices": r["unpaired_observed_indices"],
                "unpaired_expected_indices": r["unpaired_expected_indices"],
            }
        )

    return {
        "observed_fragment_lengths": list(observed_lengths),
        "tolerance": tolerance,
        "max_unpaired": max_unpaired,
        "candidate_count": len(results),
        "hit": bool(eligible),
        "ambiguous": ambiguous,
        "best_digest_id": best_digest_id,
        "candidates": output_candidates,
    }
