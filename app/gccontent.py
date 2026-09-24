"""固定窗口 GC 含量范围计算（只读）。

上下界口径
----------
序列仅含 IUPAC 碱基（调用方负责校验）。对每个碱基考察其可能碱基集合：

* 可能集合是 ``{G, C}`` 的子集时，该位置**必然**为 GC，计入下界——
  确定碱基 C、G 必然计入，简并码中只有 ``S={G,C}`` 必然计入；
* 可能集合与 ``{G, C}`` 有交集时，该位置**可能**为 GC，计入上界——
  简并码 R/Y/K/M/B/D/H/V/N 的可能碱基中均含 G 或 C（W={A,T} 不含），
  确定碱基 C、G 同样计入、A、T 同样不计入。

因此确定碱基 A/C/G/T 对上下界贡献相同（C/G 恒计入、A/T 恒不计入），
上下界的差异完全来自简并符号。

窗口口径
--------
* 窗口长度等于序列长度时，无论步长为何都只返回唯一一个整长窗口；
* 窗口长度小于序列长度时按步长滑动：
  - 线性：只保留完整落在 ``[0, L)`` 内的窗口，起点为
    ``0, step, 2*step, ...`` 中所有满足 ``start + W <= L`` 者，
    末端可能剩余不足一个步长的尾巴，不另补窗口；
  - 环形：起点为 ``0, step, 2*step, ...`` 中所有 ``< L`` 者，
    共 ``ceil(L/step)`` 个，最后一个窗口跨越原点回到开头。
* 窗口长度恒不超过序列长度，因此环形窗口最多覆盖一圈，
  同一物理位置不会在单个窗口内被重复计入。
"""

from __future__ import annotations

from .iupac import IUPAC_SETS

GC_BASES = frozenset("GC")

# 百分比保留的小数位；整数计数同时返回，精确值可由 计数 / 窗口长度 复原
PERCENT_DECIMALS = 6


def _gc_flags(sequence: str) -> tuple[list[bool], list[bool]]:
    """逐碱基返回 ``(必然为 GC, 可能为 GC)`` 两个标记列表。"""
    certain: list[bool] = []
    possible: list[bool] = []
    for base in sequence:
        choices = IUPAC_SETS[base]
        certain.append(choices <= GC_BASES)
        possible.append(bool(choices & GC_BASES))
    return certain, possible


def _prefix_sums(flags: list[bool]) -> list[int]:
    """前缀和：``prefix[i]`` 为 ``flags[0:i]`` 中真标记的个数。"""
    prefix = [0]
    total = 0
    for flag in flags:
        total += int(flag)
        prefix.append(total)
    return prefix


def window_starts(
    length: int, window_length: int, step: int, circular: bool
) -> list[int]:
    """按模块文档中的窗口口径返回全部窗口起点。调用方须已保证
    ``1 <= window_length <= length`` 且 ``step >= 1``。"""
    if window_length == length:
        return [0]
    if circular:
        # 起点 0, step, ... 中所有 < L 者：最后一个窗口跨过原点闭合环面
        return list(range(0, length, step))
    # 线性：只保留完整窗口，最后一个起点 <= L - W，末端可能留有短尾
    return list(range(0, length - window_length + 1, step))


def gc_content_windows(
    sequence: str, topology: str, window_length: int, step: int
) -> dict:
    """计算固定窗口的 GC 含量上下界，返回 ``{"window_count", "windows"}``。

    每个窗口含 index / start / end / length / wraps_origin、必然与可能
    GC 计数，以及由计数换算的最低/最高 GC 百分比。环形窗口的 end 为
    ``(start + W) % L``（可能小于 start），与酶切片段坐标口径一致。
    """
    length = len(sequence)
    circular = topology == "circular"

    certain_flags, possible_flags = _gc_flags(sequence)
    if circular:
        # 标志串复制一份后，任意 start < L、W <= L 的窗口都可直接切片，
        # 且切片长度恒 <= L —— 每个窗口最多覆盖一圈，同一位置不会计两次
        certain_flags = certain_flags + certain_flags
        possible_flags = possible_flags + possible_flags
    certain_prefix = _prefix_sums(certain_flags)
    possible_prefix = _prefix_sums(possible_flags)

    windows = []
    for index, start in enumerate(
        window_starts(length, window_length, step, circular)
    ):
        end_raw = start + window_length
        min_count = certain_prefix[end_raw] - certain_prefix[start]
        max_count = possible_prefix[end_raw] - possible_prefix[start]
        windows.append(
            {
                "index": index,
                "start": start,
                "end": end_raw % length if circular else end_raw,
                "length": window_length,
                "wraps_origin": circular and end_raw > length,
                "min_gc_count": min_count,
                "max_gc_count": max_count,
                "min_gc_percent": round(
                    100.0 * min_count / window_length, PERCENT_DECIMALS
                ),
                "max_gc_percent": round(
                    100.0 * max_count / window_length, PERCENT_DECIMALS
                ),
            }
        )
    return {"window_count": len(windows), "windows": windows}
