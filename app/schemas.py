"""请求/响应模型与输入校验。

所有校验在数据库写入之前完成：任何一个酶的偏移非法都会使
整个请求以 422 被拒绝，不会有任何数据落库。
"""

from __future__ import annotations

import math
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from .iupac import VALID_BASES


def _check_iupac(value: str, what: str) -> str:
    seq = value.strip().upper()
    if not seq:
        raise ValueError(f"{what} must not be empty")
    bad = sorted(set(seq) - VALID_BASES)
    if bad:
        raise ValueError(
            f"{what} contains non-IUPAC bases: {''.join(bad)}; "
            "allowed: ACGTRYSWKMBDHVN"
        )
    return seq


class EnzymeIn(BaseModel):
    """一种限制性内切酶：识别模式 + 正/反链切口偏移。

    偏移以识别位点 5' 端（top 链坐标）为基准，合法范围为
    [0, len(pattern)]；超出范围（含负数）即非法，整个请求被拒绝。
    """

    name: str = Field(min_length=1, max_length=100, examples=["EcoRI"])
    pattern: str = Field(examples=["GAATTC"])
    cut_top: int = Field(description="top 链切口偏移，0..len(pattern)", examples=[1])
    cut_bottom: int = Field(description="bottom 链切口偏移，0..len(pattern)", examples=[5])

    @field_validator("name")
    @classmethod
    def _name_ok(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("enzyme name must not be blank")
        return v

    @field_validator("pattern")
    @classmethod
    def _pattern_ok(cls, v: str) -> str:
        return _check_iupac(v, "pattern")

    @model_validator(mode="after")
    def _offsets_ok(self) -> "EnzymeIn":
        k = len(self.pattern)
        for field_name in ("cut_top", "cut_bottom"):
            value = getattr(self, field_name)
            if not 0 <= value <= k:
                raise ValueError(
                    f"enzyme '{self.name}': {field_name}={value} is out of "
                    f"range [0, {k}] for pattern '{self.pattern}'"
                )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(examples=["AAAAGAATTCTTTT"])
    topology: Literal["linear", "circular"] = "linear"
    enzymes: list[EnzymeIn] = Field(min_length=1)

    @field_validator("sequence")
    @classmethod
    def _sequence_ok(cls, v: str) -> str:
        return _check_iupac(v, "sequence")


class VariantIn(BaseModel):
    """一个 0-based 单碱基替换：将 position 处碱基替换为 alt。

    position 上界依赖已保存记录的序列长度，因此越界检查在路由中
    （取到记录之后）完成，以 422 拒绝；此处只做非负与碱基合法性校验。
    """

    position: int = Field(ge=0, examples=[4])
    alt: str = Field(min_length=1, max_length=1, examples=["C"])

    @field_validator("alt")
    @classmethod
    def _alt_ok(cls, v: str) -> str:
        v = v.strip().upper()
        if v not in VALID_BASES:
            raise ValueError(
                f"alt '{v}' is not an IUPAC base; allowed: ACGTRYSWKMBDHVN"
            )
        return v


class VariantAnalysisRequest(BaseModel):
    """变异影响分析请求：至少一个替换；同一位置重复出现即 422 拒绝。"""

    variants: list[VariantIn] = Field(min_length=1)

    @model_validator(mode="after")
    def _no_duplicate_positions(self) -> "VariantAnalysisRequest":
        seen: set[int] = set()
        duplicates: set[int] = set()
        for v in self.variants:
            if v.position in seen:
                duplicates.add(v.position)
            seen.add(v.position)
        if duplicates:
            raise ValueError(
                "duplicate variant positions: "
                + ", ".join(str(p) for p in sorted(duplicates))
            )
        return self



class SiteOut(BaseModel):
    enzyme: str
    start: int
    strand: Literal["+", "-"]
    matched: str


class CutOut(BaseModel):
    top_cut: int
    bottom_cut: int
    end_type: Literal["blunt", "5prime", "3prime"]
    overhang: str
    enzymes: list[str]
    sites: list[SiteOut]


class EndOut(BaseModel):
    type: Literal["blunt", "5prime", "3prime", "none"]
    overhang: str


class FragmentOut(BaseModel):
    index: int
    start: int
    end: int
    length: int
    sequence: str
    wraps_origin: bool
    left_end: EndOut
    right_end: EndOut


class DigestResult(BaseModel):
    cuts: list[CutOut]
    fragments: list[FragmentOut]


class DigestResponse(DigestResult):
    id: int
    topology: Literal["linear", "circular"]
    sequence: str
    sequence_length: int
    enzymes: list[EnzymeIn]
    created_at: str


class DigestSummary(BaseModel):
    id: int
    topology: str
    sequence_length: int
    enzyme_count: int
    cut_count: int
    fragment_count: int
    created_at: str


# ---------- 变异影响分析 ----------


class AppliedVariantOut(BaseModel):
    position: int
    ref: str
    alt: str


class ChangedCutOut(BaseModel):
    """坐标相同、但末端或归属酶/位点发生变化的切口。"""

    before: CutOut
    after: CutOut
    changes: list[str]


class CutDiffOut(BaseModel):
    added: list[CutOut]
    removed: list[CutOut]
    changed: list[ChangedCutOut]


class ChangedFragmentOut(BaseModel):
    start: int
    before_length: int
    after_length: int
    length_delta: int


class FragmentDiffOut(BaseModel):
    """片段长度差异：before/after_lengths 为物理顺序的片段长度清单。"""

    before_lengths: list[int]
    after_lengths: list[int]
    added: list[FragmentOut]
    removed: list[FragmentOut]
    changed: list[ChangedFragmentOut]


class VariantAnalysisResponse(BaseModel):
    id: int
    sequence_length: int
    variant_sequence: str
    variants: list[AppliedVariantOut]
    cuts: CutDiffOut
    fragments: FragmentDiffOut


# ---------- 酶组合面板规划 ----------


class EnzymePanelPlanRequest(BaseModel):
    """酶组合面板规划请求：目标片段数 + 1~8 种候选酶。

    候选酶沿用 EnzymeIn 的全部校验（IUPAC 模式、切口偏移范围），
    任一候选非法即整个请求以 422 拒绝。
    """

    target_fragment_count: int = Field(ge=1, examples=[3])
    candidates: list[EnzymeIn] = Field(min_length=1, max_length=8)


class EnzymePanelPlanResponse(BaseModel):
    """规划结果：选中的最少酶集合、其切口与达到的片段数（恒等于目标值）。"""

    id: int
    target_fragment_count: int
    selected_enzymes: list[EnzymeIn]
    fragment_count: int
    cuts: list[CutOut]


# ---------- 部分酶切概率分布 ----------


class CutProbabilityIn(BaseModel):
    """一个已保存物理切口 (top_cut, bottom_cut) 的独立切割概率。

    坐标身份与完整酶切结果一致；概率取闭区间 [0, 1]（0 表示必不切、
    1 表示必切），NaN/无穷大非法。坐标与已保存切口是否一一匹配
    （缺失/多余）依赖记录内容，在路由取到记录后校验 -> 422。
    """

    top_cut: int = Field(examples=[5])
    bottom_cut: int = Field(examples=[9])
    probability: float = Field(ge=0.0, le=1.0, examples=[0.8])

    @field_validator("probability")
    @classmethod
    def _probability_finite(cls, v: float) -> float:
        if not math.isfinite(v):
            raise ValueError("probability must be a finite number")
        return v


class PartialDigestRequest(BaseModel):
    """部分酶切概率分布请求：为每个已保存物理切口指定独立切割概率。

    请求内同一物理坐标重复出现即 422 拒绝；允许空列表（记录本身
    无切口时，唯一的结果是整条分子不被切割）。
    """

    cuts: list[CutProbabilityIn] = Field(default_factory=list)

    @model_validator(mode="after")
    def _no_duplicate_coordinates(self) -> "PartialDigestRequest":
        seen: set[tuple[int, int]] = set()
        duplicates: set[tuple[int, int]] = set()
        for item in self.cuts:
            key = (item.top_cut, item.bottom_cut)
            if key in seen:
                duplicates.add(key)
            seen.add(key)
        if duplicates:
            raise ValueError(
                "duplicate cut coordinates: "
                + ", ".join(f"({t}, {b})" for t, b in sorted(duplicates))
            )
        return self


class PartialDigestOutcome(BaseModel):
    """一种片段长度多重集结果及其合并后的发生概率。"""

    fragment_lengths: list[int]
    fragment_count: int
    probability: float


class PartialDigestDistributionResponse(BaseModel):
    """部分酶切分布：合并后的结果按概率降序稳定排列。"""

    id: int
    topology: Literal["linear", "circular"]
    sequence_length: int
    cut_count: int
    outcome_count: int
    total_probability: float
    outcomes: list[PartialDigestOutcome]


# ---------- 甲基化阻断分析 ----------


class MethylationRuleIn(BaseModel):
    """一种酶的甲基化敏感规则：识别模式内的敏感偏移（0-based）。

    偏移以识别模式 5' 端为基准，合法范围为 [0, len(pattern) - 1]；
    上界依赖该酶在已保存记录中的模式长度，因此非负由模型校验、
    上界在路由取到记录后校验（未知酶名同样依赖记录），越界 -> 422。
    同一酶在请求中重复出现即 422。
    """

    enzyme: str = Field(min_length=1, max_length=100, examples=["EcoRI"])
    sensitive_offsets: list[Annotated[int, Field(ge=0)]] = Field(
        min_length=1, examples=[[0, 1]]
    )

    @field_validator("enzyme")
    @classmethod
    def _enzyme_name_ok(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("enzyme name must not be blank")
        return v


class MethylationAnalysisRequest(BaseModel):
    """甲基化阻断分析请求：甲基化坐标（top 链 0-based）+ 酶敏感规则。

    甲基化坐标与酶规则各自不得重复；坐标上界依赖已保存记录的序列长度，
    越界检查在路由中（取到记录之后）完成，以 422 拒绝。
    """

    methylated_positions: list[Annotated[int, Field(ge=0)]] = Field(
        min_length=1, examples=[[4, 5]]
    )
    enzyme_rules: list[MethylationRuleIn] = Field(min_length=1)

    @model_validator(mode="after")
    def _no_duplicates(self) -> "MethylationAnalysisRequest":
        problems = []

        seen: set[int] = set()
        duplicates: set[int] = set()
        for position in self.methylated_positions:
            if position in seen:
                duplicates.add(position)
            seen.add(position)
        if duplicates:
            problems.append(
                "duplicate methylated positions: "
                + ", ".join(str(p) for p in sorted(duplicates))
            )

        seen_enzymes: set[str] = set()
        duplicate_enzymes: set[str] = set()
        for rule in self.enzyme_rules:
            if rule.enzyme in seen_enzymes:
                duplicate_enzymes.add(rule.enzyme)
            seen_enzymes.add(rule.enzyme)
        if duplicate_enzymes:
            problems.append(
                "duplicate enzyme rules: "
                + ", ".join(sorted(duplicate_enzymes))
            )

        if problems:
            raise ValueError("; ".join(problems))
        return self


class BlockedSiteOut(BaseModel):
    """一个被甲基化阻断的识别位点及其命中的甲基化物理坐标。"""

    top_cut: int
    bottom_cut: int
    site: SiteOut
    methylated_positions: list[int]


class MethylationAnalysisResponse(BaseModel):
    """甲基化分析结果：被阻断位点、保留/移除切口与重建片段均稳定排序。"""

    id: int
    topology: Literal["linear", "circular"]
    sequence_length: int
    methylated_positions: list[int]
    blocked_sites: list[BlockedSiteOut]
    retained_cuts: list[CutOut]
    removed_cuts: list[CutOut]
    fragments: list[FragmentOut]


# ---------- 凝胶条带比较 ----------


class GelComparisonRequest(BaseModel):
    """凝胶条带比较请求：2~8 条已保存记录 + 分辨率比例。

    digest_ids 不可重复（重复即 422）；resolution_ratio 取闭区间
    [0, 0.2]，NaN / 无穷大非法。记录是否存在依赖数据库，在路由中
    逐个确认，任一不存在 -> 404。
    """

    digest_ids: list[int] = Field(min_length=2, max_length=8, examples=[[1, 2]])
    resolution_ratio: float = Field(ge=0.0, le=0.2, examples=[0.05])

    @field_validator("resolution_ratio")
    @classmethod
    def _ratio_finite(cls, v: float) -> float:
        if not math.isfinite(v):
            raise ValueError("resolution_ratio must be a finite number")
        return v

    @model_validator(mode="after")
    def _no_duplicate_ids(self) -> "GelComparisonRequest":
        seen: set[int] = set()
        duplicates: set[int] = set()
        for digest_id in self.digest_ids:
            if digest_id in seen:
                duplicates.add(digest_id)
            seen.add(digest_id)
        if duplicates:
            raise ValueError(
                "duplicate digest ids: "
                + ", ".join(str(i) for i in sorted(duplicates))
            )
        return self


class GelBandOut(BaseModel):
    """一个聚类条带：长度范围、各泳道片段数、总碱基数与诊断标记。

    lane_fragment_counts 与请求 digest_ids 的顺序一一对应；
    diagnostic 为 True 表示该条带仅出现在一个泳道。
    """

    min_length: int
    max_length: int
    lane_fragment_counts: list[int]
    fragment_count: int
    total_bases: int
    diagnostic: bool


class GelComparisonResponse(BaseModel):
    """凝胶比较结果：条带按 max_length 降序稳定排列。"""

    digest_ids: list[int]
    resolution_ratio: float
    lane_count: int
    band_count: int
    bands: list[GelBandOut]


# ---------- 末端连接兼容性分析 ----------


class LigationCompatibilityRequest(BaseModel):
    """末端连接兼容性分析请求：2~8 条已保存记录 + 可选 include_blunt。

    digest_ids 不可重复（重复即 422）；include_blunt 默认为 True，
    为 false 时平端-平端配对不纳入结果。记录是否存在依赖数据库，
    在路由中逐个确认，任一不存在 -> 404。
    """

    digest_ids: list[int] = Field(min_length=2, max_length=8, examples=[[1, 2]])
    include_blunt: bool = Field(
        default=True,
        description="是否纳入平端-平端配对；false 时只返回黏端配对",
        examples=[True],
    )

    @model_validator(mode="after")
    def _no_duplicate_ids(self) -> "LigationCompatibilityRequest":
        seen: set[int] = set()
        duplicates: set[int] = set()
        for digest_id in self.digest_ids:
            if digest_id in seen:
                duplicates.add(digest_id)
            seen.add(digest_id)
        if duplicates:
            raise ValueError(
                "duplicate digest ids: "
                + ", ".join(str(i) for i in sorted(duplicates))
            )
        return self


class LigationEndOut(BaseModel):
    """一个参与配对的片段末端：记录 id、片段索引、端方向与该端类型/突出端。

    天然末端（type == "none"）不参与任何配对，因此不会出现在结果中。
    """

    digest_id: int
    fragment_index: int
    end: Literal["left", "right"]
    end_type: Literal["blunt", "5prime", "3prime"]
    overhang: str


class LigationPairOut(BaseModel):
    """一对可连接的末端：first 恒为请求顺序靠前记录的一端。

    compatibility 为 "sticky"（黏端）或 "blunt"（平端）；
    ligated_overhang 为接头处退火的 5' 突出端序列（取该对中 5' 端的
    突出端），平端配对时为空字符串。
    """

    first: LigationEndOut
    second: LigationEndOut
    compatibility: Literal["sticky", "blunt"]
    ligated_overhang: str


class LigationCompatibilityResponse(BaseModel):
    """连接兼容性结果：配对按（记录顺序、片段索引、端方向）稳定排列。"""

    digest_ids: list[int]
    include_blunt: bool
    record_count: int
    pair_count: int
    pairs: list[LigationPairOut]


# ---------- 片段拼接链规划 ----------


class LigationAssemblyRequest(BaseModel):
    """片段拼接链规划请求：2~8 条已保存记录 + 平端开关 + 最长链长。

    digest_ids 不可重复（重复即 422）；include_blunt 默认为 True，为 false
    时平端不能作为连接末端；max_chain_length 为拼接链允许包含的最大节点
    （片段）数，闭区间 [2, 8]。记录是否存在依赖数据库，在路由中逐个确认，
    任一不存在 -> 404；参与规划的片段总数过多（超过固定枚举上限）-> 422。
    """

    digest_ids: list[int] = Field(min_length=2, max_length=8, examples=[[1, 2]])
    include_blunt: bool = Field(
        default=True,
        description="是否允许平端连接；false 时只允许黏端连接",
        examples=[True],
    )
    max_chain_length: int = Field(
        ge=2,
        le=8,
        description="拼接链允许包含的最大节点（片段）数；每条记录最多贡献一个片段",
        examples=[8],
    )

    @model_validator(mode="after")
    def _no_duplicate_ids(self) -> "LigationAssemblyRequest":
        seen: set[int] = set()
        duplicates: set[int] = set()
        for digest_id in self.digest_ids:
            if digest_id in seen:
                duplicates.add(digest_id)
            seen.add(digest_id)
        if duplicates:
            raise ValueError(
                "duplicate digest ids: "
                + ", ".join(str(i) for i in sorted(duplicates))
            )
        return self


class AssemblyEndOut(BaseModel):
    """拼接链中一个片段端点：记录 id、片段索引、端方向与该端类型/突出端。

    end 为该端在**原片段**上的方向（``left`` / ``right``）；片段翻转摆入
    链中时，其链左端可能对应原片段的 right 端。天然末端 ``none`` 只可能
    出现在整条拼接链的最左/最右端（链内连接恒为兼容的切割末端）。
    """

    digest_id: int
    fragment_index: int
    end: Literal["left", "right"]
    end_type: Literal["blunt", "5prime", "3prime", "none"]
    overhang: str


class AssemblyNodeOut(BaseModel):
    """拼接链上的一个节点：一条已保存记录贡献的一个片段及其摆放方向。

    orientation 为 ``forward``（按原片段 left→right 摆入）或 ``reverse``
    （翻转摆入）；left_end / right_end 为该节点在拼接链上的左/右端点。
    """

    digest_id: int
    fragment_index: int
    orientation: Literal["forward", "reverse"]
    left_end: AssemblyEndOut
    right_end: AssemblyEndOut


class AssemblyConnectionOut(BaseModel):
    """拼接链中相邻两节点间的一处连接。

    from_end 恒为左侧节点的出端（链右端），to_end 恒为右侧节点的入端
    （链左端），二者来自不同记录且通过末端兼容判定；compatibility 为
    ``sticky``（黏端）或 ``blunt``（平端），ligated_overhang 为接头处
    退火的 5' 突出端序列，平端连接时为空字符串。
    """

    index: int
    from_end: AssemblyEndOut
    to_end: AssemblyEndOut
    compatibility: Literal["sticky", "blunt"]
    ligated_overhang: str


class LigationAssemblyResponse(BaseModel):
    """拼接链规划结果：节点数最多的线性链；同长按（记录请求顺序、片段索引、
    端方向）稳定择优。

    所有提交记录均不存在任何可连接末端时 nodes / connections 为空，
    stop_reason 为 ``no_compatible_ends``。used / unused_digest_ids 分别为
    入选 / 未入选链的记录 id，均按请求顺序排列。
    """

    digest_ids: list[int]
    include_blunt: bool
    max_chain_length: int
    record_count: int
    node_count: int
    connection_count: int
    nodes: list[AssemblyNodeOut]
    connections: list[AssemblyConnectionOut]
    used_digest_ids: list[int]
    unused_digest_ids: list[int]
    stop_reason: Literal[
        "no_compatible_ends",
        "all_records_used",
        "max_chain_length_reached",
        "no_compatible_extension",
    ]
    stop_reason_detail: str


# ---------- PCR 扩增产物预测 ----------


class PcrProductsRequest(BaseModel):
    """PCR 扩增产物预测请求：两条 IUPAC 引物（5'->3'）+ 最大产物长度。

    纯只读计算，读取已保存记录的原始序列与拓扑，不引用酶切片段、不写库。
    两条引物至少各 2 个碱基（3' 端最后两位强制匹配，单碱基引物无从谈起）；
    ``max_product_length`` 为正整数且不得超过模板长度——上界依赖已保存
    记录的序列长度，越界检查在路由中（取到记录之后）完成，以 422 拒绝。
    引物长于模板时不产生位点，结果为空数组（并非非法请求）。
    """

    forward_primer: str = Field(
        description="正引物序列（5'->3'），沿模板正向匹配", examples=["ATGCGAT"]
    )
    reverse_primer: str = Field(
        description="反引物序列（5'->3'），以反向互补序列在模板上匹配",
        examples=["TTCGCGC"],
    )
    max_product_length: int = Field(
        ge=1, description="最大产物长度（含），不得超过模板长度", examples=[1000]
    )

    @field_validator("forward_primer", "reverse_primer")
    @classmethod
    def _primer_ok(cls, v: str) -> str:
        seq = _check_iupac(v, "primer")
        # 在去除空白后的规范序列上检查长度，使 " A" 这类输入同样被拒绝
        if len(seq) < 2:
            raise ValueError("primer must contain at least 2 IUPAC bases")
        return seq


class PcrBindingOut(BaseModel):
    """一条引物的结合信息：top 链 0-based 半开区间 ``[start, end)`` 与错配数。

    反引物的坐标同样给出其反向互补序列在 top 链上的结合区；环形结合区
    跨原点时 ``end`` 为取模后的坐标（可能小于 ``start``）。``mismatches``
    恒为 0 或 1（3' 端最后两位不计入错配——该位置不兼容的位点不会被枚举）。
    """

    start: int
    end: int
    mismatches: int


class PcrProductOut(BaseModel):
    """一个 PCR 产物：top 链 5'->3' 序列、两端结合坐标、长度、各自错配数。

    ``start`` / ``end`` 为产物在 top 链上的 0-based 半开区间端点；环形产物
    跨原点时 ``end`` 为取模后的坐标（可能小于 ``start``）且
    ``wraps_origin`` 为 true。
    """

    sequence: str
    start: int
    end: int
    length: int
    wraps_origin: bool
    forward_binding: PcrBindingOut
    reverse_binding: PcrBindingOut


class PcrProductsResponse(BaseModel):
    """PCR 产物预测结果：产物按（起点、长度）升序排列，无产物时为空数组。"""

    id: int
    topology: Literal["linear", "circular"]
    sequence_length: int
    forward_primer: str
    reverse_primer: str
    max_product_length: int
    product_count: int
    products: list[PcrProductOut]


# ---------- 开放阅读框（ORF）检索 ----------


class OrfSearchRequest(BaseModel):
    """开放阅读框检索请求：最短氨基酸数（不含终止密码子）。

    纯只读计算，读取已保存记录的原始序列与拓扑，不引用酶切片段、不写库。
    ``min_aa_length`` 为正整数；无任何依赖记录内容的校验，记录不存在
    时先返回 404。
    """

    min_aa_length: int = Field(
        ge=1, description="最短氨基酸数（不含终止密码子，含起始甲硫氨酸）",
        examples=[20],
    )


class OrfOut(BaseModel):
    """一个开放阅读框：链方向、top 链坐标、定向序列与长度。

    ``start`` / ``end`` 为该 ORF 在 top 链（原始正链）上覆盖区间的
    0-based 半开端点；环形跨原点时 ``end`` 为取模后的坐标（可能小于
    ``start``）且 ``wraps_origin`` 为 true。``sequence`` 为 ORF 自身
    5'->3' 的定向核苷酸序列（含终止密码子）；``length`` 为其核苷酸
    长度；``aa_length`` 为氨基酸数（不含终止）。
    """

    strand: Literal["+", "-"]
    start: int
    end: int
    wraps_origin: bool
    sequence: str
    length: int
    aa_length: int


class OrfSearchResponse(BaseModel):
    """ORF 检索结果：按（链方向、起点、长度）稳定排序，无结果时为空数组。"""

    id: int
    topology: Literal["linear", "circular"]
    sequence_length: int
    min_aa_length: int
    orf_count: int
    orfs: list[OrfOut]


# ---------- 未知样本识别（片段长度配对） ----------


class DigestIdentificationRequest(BaseModel):
    """未知样本识别请求：1~8 条已保存记录 + 观测片段长度 + 容差 + 最大未配对数。

    纯只读计算，读取各记录酶切片段的长度，不引用原始序列、不写库。
    ``digest_ids`` 不可重复（重复即 422）；``observed_fragment_lengths``
    为 1~12 个正整数（凝胶上读得的条带长度，重复长度按不同片段计数）；
    ``tolerance`` 为非负整数（bp）；``max_unpaired`` 为非负整数。
    记录是否存在、其片段数是否超过 12，均依赖数据库内容，在路由中
    逐个确认，任一不存在 -> 404。
    """

    digest_ids: list[int] = Field(min_length=1, max_length=8, examples=[[1, 2]])
    observed_fragment_lengths: list[Annotated[int, Field(ge=1)]] = Field(
        min_length=1,
        max_length=12,
        description="未知样本的观测片段长度（bp），1~12 个正整数，重复长度按不同片段计数",
        examples=[[100, 200, 305]],
    )
    tolerance: int = Field(
        ge=0,
        description="允许配对的最大长度差（bp，含），非负整数",
        examples=[10],
    )
    max_unpaired: int = Field(
        ge=0,
        description="命中所允许的两侧未配对片段总数上限（含），非负整数",
        examples=[2],
    )

    @model_validator(mode="after")
    def _no_duplicate_ids(self) -> "DigestIdentificationRequest":
        seen: set[int] = set()
        duplicates: set[int] = set()
        for digest_id in self.digest_ids:
            if digest_id in seen:
                duplicates.add(digest_id)
            seen.add(digest_id)
        if duplicates:
            raise ValueError(
                "duplicate digest ids: "
                + ", ".join(str(i) for i in sorted(duplicates))
            )
        return self


class IdentificationPairOut(BaseModel):
    """一对配对成功的片段：观测/预期索引、长度与差值。

    索引分别对应请求的观测长度列表与 ``GET /digests/{id}`` 结果中的
    片段索引（0-based，按物理顺序）；``length_diff`` 为
    ``|observed_length - expected_length|``，恒不超过容差。
    """

    observed_index: int
    expected_index: int
    observed_length: int
    expected_length: int
    length_diff: int


class IdentificationCandidateOut(BaseModel):
    """一个候选记录的最优配对结果。

    ``unpaired_count`` 为两侧未配对片段总数
    （观测数 + 预期数 - 2 × 配对数）；``eligible`` 表示其未超过
    ``max_unpaired``（计入命中）；``best`` 表示其分数（未配对总数,
    差值总和）在全部命中候选中并列最优。配对按观测索引升序，
    未配对索引按各自列表的索引升序。
    """

    digest_id: int
    fragment_count: int
    pair_count: int
    unpaired_count: int
    diff_sum: int
    eligible: bool
    best: bool
    pairs: list[IdentificationPairOut]
    unpaired_observed_indices: list[int]
    unpaired_expected_indices: list[int]


class DigestIdentificationResponse(BaseModel):
    """识别结果：候选按（未配对总数, 差值总和, 请求 ID 顺序）稳定排序。

    无候选满足 ``max_unpaired`` 时 ``hit`` 为 false、``ambiguous`` 为
    false、``best_digest_id`` 为 null（无命中）；有命中且最优分数被
    多个候选并列时 ``ambiguous`` 为 true，``best_digest_id`` 取请求
    顺序最靠前的并列候选作为代表。
    """

    digest_ids: list[int]
    observed_fragment_lengths: list[int]
    tolerance: int
    max_unpaired: int
    candidate_count: int
    hit: bool
    ambiguous: bool
    best_digest_id: int | None
    candidates: list[IdentificationCandidateOut]


# ---------- 固定窗口 GC 含量范围 ----------


class GcContentRequest(BaseModel):
    """固定窗口 GC 含量范围请求：IUPAC 序列 + 窗口长度 + 步长。

    纯只读计算，不引用任何已保存记录、不写数据库。``window_length`` /
    ``step`` 为正整数，且窗口长度不得超过序列长度（越界 -> 422）；
    窗口长度等于序列长度时只返回一个整长窗口（step 不参与取窗）。
    """

    sequence: str = Field(examples=["ATGCSWNATG"])
    topology: Literal["linear", "circular"] = "linear"
    window_length: int = Field(ge=1, description="固定窗口长度（碱基数）", examples=[4])
    step: int = Field(ge=1, description="窗口小于序列长度时的滑动步长", examples=[2])

    @field_validator("sequence")
    @classmethod
    def _sequence_ok(cls, v: str) -> str:
        return _check_iupac(v, "sequence")

    @model_validator(mode="after")
    def _window_fits_sequence(self) -> "GcContentRequest":
        if self.window_length > len(self.sequence):
            raise ValueError(
                f"window_length {self.window_length} exceeds sequence length "
                f"{len(self.sequence)}"
            )
        return self


class GcWindowOut(BaseModel):
    """单个窗口的 GC 含量范围。

    min/max_gc_count 为窗口内"必然为 GC"/"可能为 GC"的碱基数：
    确定碱基 C/G 两者都计入、A/T 两者都不计；简并符号按其可能碱基
    分别计入上/下界（仅 S={G,C} 计入下界）。min/max_gc_percent 为
    相应计数除以窗口长度的百分比，保留 6 位小数。环形窗口的 end 为
    取模后的坐标（可能小于 start），``wraps_origin`` 标明是否跨原点。
    """

    index: int
    start: int
    end: int
    length: int
    wraps_origin: bool
    min_gc_count: int
    max_gc_count: int
    min_gc_percent: float
    max_gc_percent: float


class GcContentResponse(BaseModel):
    """GC 含量窗口结果：窗口按起点升序排列（环形最后一个窗口可能跨原点）。"""

    topology: Literal["linear", "circular"]
    sequence_length: int
    window_length: int
    step: int
    window_count: int
    windows: list[GcWindowOut]

