# DNA 限制性酶切模拟 API

基于 **Python + FastAPI + SQLite** 的限制性内切酶消化模拟服务。提交线性或环形 DNA
序列（仅 IUPAC 碱基）以及多种酶的识别模式与正/反链切口偏移，服务在**两个方向**
查找识别位点（允许重叠、含简并碱基、可跨越环形序列首尾），对同一物理切口去重后
按拓扑结构生成片段，返回坐标、长度、两端的平端/黏端类型及突出端序列。

## 快速启动

一条指令完成构建、依赖安装与数据库初始化：

```bash
docker compose up --build
```

依赖在镜像构建时安装（`pip install -r requirements.txt`），SQLite 数据库在容器
启动时自动建库建表，无需任何手工步骤。服务监听容器内 **8000** 端口，compose 通过
`ports: ["8000:8000"]` 映射到宿主机同一端口——启动后直接访问
**http://localhost:8000/docs** 即可看到包含全部接口的交互式文档（如宿主机 8000
端口被占用，可把 compose 中映射改为 `8001:8000` 后从 8001 访问）。

其它常用 Docker 操作：

```bash
# 后台运行
docker compose up --build -d

# 查看日志 / 停止
docker compose logs -f
docker compose down

# 仅构建镜像
docker build -t digest-api .
docker run --rm -p 8000:8000 -v "$PWD/data:/app/data" digest-api

# 在容器内运行测试套件
docker build -t digest-api .
docker run --rm digest-api sh -c "pip install --no-cache-dir pytest httpx && pytest"
```

## 配置

| 配置项 | 默认值 | 说明 |
| --- | --- | --- |
| 环境变量 `DATABASE_PATH` | `/app/data/digest.db`（容器内） | SQLite 数据库文件路径，目录不存在时自动创建 |
| 端口映射 | `8000:8000`（compose 中 `ports`） | 宿主机端口：容器端口 |
| 数据卷 | `digest-data` → `/app/data` | 持久化数据库文件 |

不使用 Docker 时本地运行：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## API 地址

- 服务地址：**http://localhost:8000**
- 交互式文档（Swagger UI）：**http://localhost:8000/docs**
- 健康检查：`GET /health`

### 接口一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/digests` | 提交序列与酶列表，执行酶切、保存并返回完整结果（201） |
| `GET` | `/digests` | 列出全部提交的摘要（按 id 升序） |
| `GET` | `/digests/{id}` | 按 id 查询保存的完整酶切结果（不存在返回 404） |
| `POST` | `/digests/{id}/variant-analysis` | 对已保存记录施加一组单碱基替换并重跑酶切，返回切口与片段长度差异（不落库，记录不存在返回 404） |
| `POST` | `/digests/{id}/enzyme-panel-plan` | 枚举候选酶组合，寻找达到目标片段数的最少酶集合（只读不落库，记录不存在返回 404，无解返回 422） |
| `POST` | `/digests/{id}/partial-digest-distribution` | 为每个已保存物理切口指定独立切割概率，枚举成功/失败状态并返回片段长度多重集的概率分布（只读不落库，记录不存在返回 404） |
| `POST` | `/digests/{id}/methylation-analysis` | 只读分析甲基化对已保存酶切位点的阻断影响：按酶的敏感偏移定位甲基化碱基，阻断/保留/移除切口并重建片段（不落库，记录不存在返回 404） |
| `POST` | `/digests/{id}/pcr-products` | 只读预测 PCR 扩增产物：读取已保存记录的原始序列与拓扑（不依赖酶切片段），按 IUPAC 集合交集匹配正/反引物（各允许一处错配、3' 端最后两位强制匹配），枚举相向且结合区不重叠的位点对，线性不越界、环形可跨原点但不超过一圈，返回产物序列、结合坐标、长度、错配数与跨原点标记（不落库，记录不存在返回 404，参数非法返回 422，无产物返回空数组） |
| `POST` | `/digests/{id}/orfs` | 只读检索开放阅读框：基于已保存记录的原始 IUPAC 序列与拓扑（不依赖酶切片段），正链与反向互补链分别逐碱基寻找字面量 ATG，从各起点按三联体读取到同框第一个字面量 TAA/TAG/TGA 形成 ORF，简并密码子不充当起止信号，线性不越界、环形可跨原点但完整 ORF（含终止）不超过一圈，保留满足最短氨基酸数（不含终止）的全部嵌套 ORF，返回链方向、原始正链坐标、跨原点标记、定向核苷酸序列与氨基酸长度（不落库，记录不存在返回 404，参数非法返回 422，无结果返回空数组） |
| `POST` | `/gel-comparisons` | 只读比较 2~8 条已保存记录的凝胶条带：汇总片段长度按降序聚类（条带内 `(max-min)/max <= resolution_ratio`），返回条带长度范围、各泳道片段数、总碱基数与诊断条带（不落库，记录不存在返回 404） |
| `POST` | `/ligation-compatibility` | 只读枚举 2~8 条已保存记录之间可连接的片段末端对：黏端要求一方 5' 突出端序列等于另一方 3' 突出端序列的反向互补，平端仅与平端兼容（可通过 `include_blunt` 关闭），返回按记录/片段/端方向稳定排序的配对、兼容类型、连接后突出端与总数（不落库，记录不存在返回 404） |
| `POST` | `/ligation-assembly` | 只读规划 2~8 条已保存记录之间的最长片段拼接链：复用末端反向互补兼容判定构建有向图，搜索不重复使用记录的最长线性链（每条记录最多贡献一个片段，可正向/反向摆入，`max_chain_length` 限制链上节点数），同长度方案按记录请求顺序、片段索引与端方向稳定择优，返回节点端点、连接类型/突出端、总连接数与未能延伸的原因（不落库，记录不存在返回 404，候选片段过多返回 422） |
| `POST` | `/gc-content` | 只读对提交的 IUPAC 序列按固定窗口计算 GC 含量范围：每窗口返回确定碱基贡献的最低/最高 GC 计数与百分比，简并符号按其可能碱基分别计入上/下界；窗口短于序列时按指定步长滑动，环形窗口可跨越末尾回到开头且每个窗口最多覆盖一圈（不落库，不依赖已保存记录） |
| `POST` | `/digest-identification` | 只读从未知样本的观测片段长度识别 1~8 条已保存记录：观测与预期片段一对一配对、长度差不超过整数 bp 容差，先最大化配对数再最小化差值总和，同分配对按索引字典序固定；候选按未配对总数、差值总和、请求 ID 顺序排序，最优分数并列标记歧义，无候选满足最大未配对数时返回无命中（不落库） |
| `GET` | `/health` | 健康检查 |

### 提交示例

```bash
curl -X POST http://localhost:8000/digests \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAAAGAATTCTTTT",
    "topology": "linear",
    "enzymes": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    ]
  }'
```

- `sequence`：仅允许 IUPAC 碱基 `ACGTRYSWKMBDHVN`（大小写均可，存储为大写）。
- `topology`：`linear`（线性）或 `circular`（环形）。
- `enzymes`：至少一种酶。
  - `pattern`：识别模式，同样仅允许 IUPAC 碱基；
  - `cut_top` / `cut_bottom`：正/反链切口偏移，以识别位点 5' 端为基准，
    合法范围为 **`0 ~ len(pattern)`**。

响应（节选）：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 14,
  "cuts": [
    {
      "top_cut": 5,
      "bottom_cut": 9,
      "end_type": "5prime",
      "overhang": "AATT",
      "enzymes": ["EcoRI"],
      "sites": [{"enzyme": "EcoRI", "start": 4, "strand": "+", "matched": "GAATTC"}]
    }
  ],
  "fragments": [
    {"index": 0, "start": 0, "end": 5, "length": 5, "sequence": "AAAAG",
     "wraps_origin": false,
     "left_end": {"type": "none", "overhang": ""},
     "right_end": {"type": "5prime", "overhang": "AATT"}},
    {"index": 1, "start": 5, "end": 14, "length": 9, "sequence": "AATTCTTTT",
     "wraps_origin": false,
     "left_end": {"type": "5prime", "overhang": "AATT"},
     "right_end": {"type": "none", "overhang": ""}}
  ]
}
```

## 变异影响分析

对已保存的酶切记录施加一组 **0-based 单碱基替换**，在**内存中**用变异后序列
重跑酶切（拓扑与酶列表沿用原记录），并与原结果逐项对比。**分析过程不写库**：
原记录保持不变，也不会产生新的提交。

### 请求

```bash
curl -X POST http://localhost:8000/digests/1/variant-analysis \
  -H 'Content-Type: application/json' \
  -d '{"variants": [{"position": 4, "alt": "C"}]}'
```

- `variants`：至少一个替换；每项 `position` 为 0-based 位置，`alt` 为单个
  IUPAC 碱基（大小写均可，存储为大写）。
- 允许同一碱基的"同义替换"（`alt` 与原碱基相同），此时差异全部为空。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| `position` 越界（`>= 序列长度`，含上界） | **422** |
| `position` 为负数 / `alt` 为非 IUPAC 碱基 | **422** |
| 同一 `position` 重复出现 | **422** |
| `variants` 为空列表 | **422** |

### 响应与差异口径

```json
{
  "id": 1,
  "sequence_length": 14,
  "variant_sequence": "AAAACAATTCTTTT",
  "variants": [{"position": 4, "ref": "G", "alt": "C"}],
  "cuts": {
    "added": [],
    "removed": [ {"top_cut": 5, "bottom_cut": 9, "...": "（完整切口对象）"} ],
    "changed": []
  },
  "fragments": {
    "before_lengths": [5, 9],
    "after_lengths": [14],
    "added": [],
    "removed": [ {"start": 5, "...": "（完整片段对象）"} ],
    "changed": [
      {"start": 0, "before_length": 5, "after_length": 14, "length_delta": 9}
    ]
  }
}
```

- **切口身份 = 物理坐标 `(top_cut, bottom_cut)`**：
  - `added`：仅变异后存在的切口（**新增**）；
  - `removed`：仅原结果存在的切口（**消失**）；
  - `changed`：坐标相同、但末端或归属发生变化，`before` / `after` 为两个完整
    切口对象，`changes` 列出变化字段（`end_type`、`overhang`、`enzymes`、
    `sites` 的子集）——例如简并位点存活但突出端序列被替换改变；
  - 切口"移动位置"（坐标/末端位置改变）表现为一条 `removed` 加一条 `added`。
- **片段差异**：替换不改变序列长度，片段按骨架起点 `start` 配对。`changed`
  给出起点相同片段的长度差异（`length_delta = after_length - before_length`）；
  `added` / `removed` 为仅一侧存在的完整片段；`before_lengths` /
  `after_lengths` 为物理顺序的片段长度清单，便于快速对比。
- **稳定排序**：回显变异按 `position` 升序（与传入顺序无关）；切口差异按
  `(top_cut, bottom_cut)` 排序；片段差异按 `start` 排序，相同请求结果必然一致。

## 酶组合面板规划

基于已保存记录的**序列与拓扑**，枚举请求中给出的候选酶（**至多 8 种**）的所有
组合并复用酶切计算，寻找**恰好达到目标片段数的最少酶集合**。**规划过程不写库**：
原记录保持不变，也不会产生新的提交。

### 请求

以下示例假设 `id=1` 的记录序列为 `AAAAGAATTCTTTTGGATCCGGGG`（线性，含 EcoRI、
BamHI 各一个位点；可先通过 `POST /digests` 创建）：

```bash
curl -X POST http://localhost:8000/digests/1/enzyme-panel-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "target_fragment_count": 3,
    "candidates": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "BamHI", "pattern": "GGATCC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    ]
  }'
```

- `target_fragment_count`：目标片段数，整数且 `>= 1`。
- `candidates`：1~8 种候选酶，字段与校验规则和提交酶切时的 `enzymes` 完全一致；
  完全相同（名称/模式/偏移均相同）的候选会自动去重。

### 响应

```json
{
  "id": 1,
  "target_fragment_count": 3,
  "selected_enzymes": [
    {"name": "BamHI", "pattern": "GGATCC", "cut_top": 1, "cut_bottom": 5},
    {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
  ],
  "fragment_count": 3,
  "cuts": [
    {"top_cut": 5, "bottom_cut": 9, "end_type": "5prime", "overhang": "AATT",
     "enzymes": ["EcoRI"], "sites": ["..."]},
    {"top_cut": 15, "bottom_cut": 19, "end_type": "5prime", "overhang": "GATC",
     "enzymes": ["BamHI"], "sites": ["..."]}
  ]
}
```

- `selected_enzymes`：选中的最少酶集合，按酶名排序返回；`fragment_count` 为其实际
  达到的片段数（恒等于目标值）；`cuts` 为该组合切口的完整对象（含坐标、末端类型、
  归属酶与位点），与 `POST /digests` 结果中的切口结构一致。

### 择优规则与状态码

- **最少酶优先**：按组合大小从小到大枚举，第一个达到目标片段数的规模即为答案；
  空集也参与枚举——不切割时线性/环形分子均为 1 条片段，因此目标为 1 时可能返回
  空的 `selected_enzymes`。
- **酶数相同时按酶名序列择优**：比较排序后的酶名元组，字典序最小者胜出；候选在
  请求中的传入顺序不影响结果，相同请求必然得到相同方案。
- **无解返回 422**：所有候选组合（含空集）都无法恰好达到目标片段数时，返回
  `422` 并说明原因。

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| 所有候选组合均无法达到目标片段数 | **422** |
| `target_fragment_count` < 1 | **422** |
| `candidates` 为空或超过 8 种 | **422** |
| 候选酶非法（非 IUPAC 模式 / 切口偏移越界） | **422** |

## 部分酶切概率分布

对一次**已保存的完整酶切结果**，为其中每个**物理切口**指定一个 0～1 的**独立
切割概率**（partial digest：该切口以给定概率切成功，否则保持未切）。服务枚举
全部成功/失败状态（**至多 12 个切口，最多 4096 种状态**），复用线性/环形片段
重建逻辑生成产物，把**片段长度多重集相同**的状态合并并累加概率，按概率**降序
稳定返回**。**分析过程不写库**：原记录保持不变，也不会产生新的提交。

### 请求

```bash
curl -X POST http://localhost:8000/digests/1/partial-digest-distribution \
  -H 'Content-Type: application/json' \
  -d '{
    "cuts": [
      {"top_cut": 5, "bottom_cut": 9, "probability": 0.8}
    ]
  }'
```

- `cuts`：每个**已保存物理切口**恰好一项（坐标身份为 `(top_cut, bottom_cut)`，
  与 `GET /digests/{id}` 结果一致）：
  - `probability`：该切口独立切成功的概率，闭区间 **`[0, 1]`**；`0` 表示必不切、
    `1` 表示必切，`NaN` / 无穷大 / 越界均为 **422**；
  - 请求条目的顺序不影响结果（按已保存切口 `(top_cut, bottom_cut)` 升序枚举）。
- 记录本身没有切口时，允许传空列表 `[]`，得到概率为 1 的唯一全长结果。

### 响应

以上例（记录序列 `AAAAGAATTCTTTT`，EcoRI 单切口 `(5,9)`）为例：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 14,
  "cut_count": 1,
  "outcome_count": 2,
  "total_probability": 1.0,
  "outcomes": [
    {"fragment_lengths": [5, 9], "fragment_count": 2, "probability": 0.8},
    {"fragment_lengths": [14], "fragment_count": 1, "probability": 0.2}
  ]
}
```

- `outcomes[]`：
  - `fragment_lengths`：**排序后的片段长度多重集**——不同切割状态只要产物长度
    分布相同（如在线性分子上只切对称位置切口之一）就合并为一项；
  - `fragment_count`：片段数（= 长度多重集大小）；
  - `probability`：该长度分布的合并概率。
- `total_probability`：所有返回结果概率之和（独立事件全枚举，数学上恒为 1），
  服务以 `math.fsum` 精确累加，**与 1 的误差远小于 1e-9**。
- **稳定排序**：先按 `probability` 降序；概率相同的按长度多重集字典序排序，
  相同请求（无论 `cuts` 条目顺序）结果必然一致。概率为 0 的不可能结果不返回。
- **拓扑规则与完整酶切完全一致**：线性分子的端点切口不产生新边界；环形分子
  单切口开环为 1 条全长片段、无切口保留整条环（因此环形分子的"无切口"与
  "恰好一个切口成功"会合并到同一长度多重集 `[L]`）。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| 请求缺少某个已保存切口的概率 | **422** |
| 请求给出了已保存结果中不存在的切口坐标 | **422** |
| 请求内同一坐标重复出现 | **422** |
| 概率不在 `[0, 1]`（含 `NaN` / 无穷大） | **422** |
| 已保存切口数超过 12（枚举上限 4096 状态） | **422** |

## 甲基化阻断分析

对一次**已保存的完整酶切结果**，给出一组 **top 链 0-based 甲基化坐标**以及每种酶的
**甲基化敏感偏移规则**，服务在**内存中**评估甲基化对已保存识别位点的阻断效果：
任一敏感位置被甲基化即阻断该位点；同一物理切口只要还有未阻断的位点就保留（并重算
其 `enzymes` / `sites`），否则移除；最后由保留切口重建片段。**分析过程不写库**：
原记录保持不变，也不会产生新的提交。

### 坐标映射

敏感偏移是相对于识别模式 5' 端的 0-based 下标，按位点所在链换算为 top 链物理坐标：

- **正向位点**（`strand: "+"`）：`start + offset`；
- **反向位点**（`strand: "-"`）：`start + 模式长度 - 1 - offset`；
- **环形序列**：上述坐标对序列长度取模，因此跨原点位点的偏移也能定位到真实碱基；
  线性序列不做取模（位点本身不出界）。

位点的任一敏感物理坐标出现在 `methylated_positions` 中，即视为该位点**被阻断**。

### 请求

以下示例假设 `id=1` 的记录序列为 `GAATTCAAAAGAATTC`（线性，两个 EcoRI 位点，
分别起于 0 和 10，切口 `(1,5)` 与 `(11,15)`；可先通过 `POST /digests` 创建）：

```bash
curl -X POST http://localhost:8000/digests/1/methylation-analysis \
  -H 'Content-Type: application/json' \
  -d '{
    "methylated_positions": [0],
    "enzyme_rules": [
      {"enzyme": "EcoRI", "sensitive_offsets": [0]}
    ]
  }'
```

- `methylated_positions`：**top 链 0-based 物理坐标**列表，至少一个，不可重复；
  坐标上界为已保存记录的序列长度（越界返回 422）。
- `enzyme_rules`：至少一条规则，同一酶名不可重复；
  - `enzyme`：酶名，必须是该记录已保存的酶（未知酶名返回 422）；
  - `sensitive_offsets`：识别模式内的敏感偏移（0-based），至少一个，合法范围
    **`0 ~ len(pattern)-1`**（越界返回 422）；同一规则内重复的偏移允许出现，
    会自动按集合处理。
- 请求条目的顺序不影响结果（甲基化坐标回显时排序去重）。

### 响应

上例中起于 0 的位点在其偏移 0（物理坐标 0）被甲基化而阻断，切口 `(1,5)` 移除；
起于 10 的位点不受影响，切口 `(11,15)` 保留，片段由 `[1,10,5]` 重建为 `[11,5]`：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 16,
  "methylated_positions": [0],
  "blocked_sites": [
    {
      "top_cut": 1,
      "bottom_cut": 5,
      "site": {"enzyme": "EcoRI", "start": 0, "strand": "+", "matched": "GAATTC"},
      "methylated_positions": [0]
    }
  ],
  "retained_cuts": [
    {
      "top_cut": 11,
      "bottom_cut": 15,
      "end_type": "5prime",
      "overhang": "AATT",
      "enzymes": ["EcoRI"],
      "sites": [
        {"enzyme": "EcoRI", "start": 10, "strand": "+", "matched": "GAATTC"}
      ]
    }
  ],
  "removed_cuts": [
    {
      "top_cut": 1,
      "bottom_cut": 5,
      "end_type": "5prime",
      "overhang": "AATT",
      "enzymes": ["EcoRI"],
      "sites": [
        {"enzyme": "EcoRI", "start": 0, "strand": "+", "matched": "GAATTC"}
      ]
    }
  ],
  "fragments": [
    {"index": 0, "start": 0, "end": 11, "length": 11, "sequence": "GAATTCAAAAG",
     "wraps_origin": false,
     "left_end": {"type": "none", "overhang": ""},
     "right_end": {"type": "5prime", "overhang": "AATT"}},
    {"index": 1, "start": 11, "end": 16, "length": 5, "sequence": "AATTC",
     "wraps_origin": false,
     "left_end": {"type": "5prime", "overhang": "AATT"},
     "right_end": {"type": "none", "overhang": ""}}
  ]
}
```

- `blocked_sites[]`：每个被阻断的识别位点一项，含其所属物理切口坐标、完整位点对象
  与实际命中的 `methylated_positions`（一个位点可能同时命中多个甲基化坐标）。
- `retained_cuts[]`：仍存在未阻断位点的物理切口，结构与完整酶切的切口一致，
  其中 `enzymes` / `sites` 已按存活位点**重新合并去重**。
- `removed_cuts[]`：全部位点均被阻断而消失的物理切口，回显阻断前的完整切口对象。
- `fragments[]`：由保留切口复用完整酶切的片段逻辑重建，线性端点切口、环形单切口
  开环、无切口保留整条环等规则与 `POST /digests` 完全一致。
- **稳定排序**：`retained_cuts` / `removed_cuts` 按 `(top_cut, bottom_cut)` 排序；
  `blocked_sites` 按 `(top_cut, bottom_cut, 酶名, 起点, 链)` 排序；片段按起点排序。
  相同请求（无论列表条目顺序）结果必然一致。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| 甲基化坐标越界（`>= 序列长度`，含上界） | **422** |
| 甲基化坐标为负数 / 重复出现 | **422** |
| `enzyme_rules` 中出现记录未保存的酶名 | **422** |
| 同一酶名在 `enzyme_rules` 中重复出现 | **422** |
| 敏感偏移超出 `[0, len(pattern)-1]`（含负数） | **422** |
| `methylated_positions` 或 `enzyme_rules` 为空列表 | **422** |



- **坐标**：均为 0-based，作用于 top 链（提交序列方向）；切口位置 `c` 表示断裂
  位于第 `c-1` 与第 `c` 个碱基之间。片段区间为半开区间 `[start, end)`。
- **末端类型 `end_type` / `left_end.type` / `right_end.type`**：
  - `blunt`：平端（正/反链切口位置相同）；
  - `5prime`：5' 黏端；`3prime`：3' 黏端；
  - `none`：线性分子的天然末端，或未切割的完整分子。
- **`overhang`**：突出端序列，按 top 链 5'→3' 方向取自两切口之间的区段；
  平端与 `none` 时为空字符串。
- **环形序列**：切口位置对序列长度取模；`wraps_origin: true` 表示该片段跨越
  首尾原点，其 `end` 为取模后的坐标（可能小于 `start`）。单个切口开环时返回
  一条全长片段；无切口时保留整条环状分子。
- **结果顺序**：切口按 `(top_cut, bottom_cut)` 排序，片段按起点排序，酶名与
  位点列表均排序——**相同输入的结果顺序必然一致**。

## PCR 扩增产物预测

基于已保存记录的**原始序列与线性/环形拓扑**（**不依赖酶切片段**，即使记录没有
任何切口也可使用），按请求给出的**正、反两条 IUPAC 引物**（均按 5'→3' 书写）
预测 PCR 扩增产物。**预测过程不写库**：原记录保持不变，也不会产生新的提交。

### 匹配与产物口径

- **正引物**结合 top 链、沿模板**正向**（向右）延伸；**反引物**先取**反向互补**
  再在 top 链上匹配、沿模板**反向**（向左）延伸。
- 引物碱基与模板碱基的 **IUPAC 可能集合有交集**即匹配（如引物 `A` 可匹配模板
  `R={A,G}`）；每条引物整个结合区允许**至多 1 处错配**，但各自 **3' 端最后两位
  必须匹配**（该两位任意一处不兼容即整条结合不成立）。
- 枚举正、反结合位点中**相向且结合区互不重叠**的位点对：
  - **线性**：正结合区完整位于反结合区之前（`forward.start + 引物长 <= reverse.start`），
    产物区间不得越出模板边界；
  - **环形**：以正位点起点沿环向前到反位点起点的有向距离判定，允许产物**跨原点**
    （`wraps_origin: true`），但长度**不得超过一圈**（即不超过模板长度）。
- 长度超过 `max_product_length` 的产物不返回；产物按**起点、长度**去重并升序排列。

### 请求

以下示例假设 `id=1` 的记录为 28 bp 线性序列
`ATGCGATCGCTTTTTTTTTTTTAGCTGC`（可先通过 `POST /digests` 创建；模板上无需任何
酶切位点）：

```bash
curl -X POST http://localhost:8000/digests/1/pcr-products \
  -H 'Content-Type: application/json' \
  -d '{
    "forward_primer": "ATGCGATCGC",
    "reverse_primer": "GCAGCTAA",
    "max_product_length": 28
  }'
```

- `forward_primer` / `reverse_primer`：5'→3' 书写的 IUPAC 引物，至少 **2** 个
  碱基（大小写均可，按大写处理）；引物长于模板时不产生结合位点，结果为空数组。
- `max_product_length`：最大产物长度（含），正整数，**不得超过模板长度**。

### 响应

反引物 `GCAGCTAA` 的反向互补为 `TTAGCTGC`，在位置 20 精确结合；正引物在位置 0
精确结合，唯一产物覆盖 `[0, 28)`：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 28,
  "forward_primer": "ATGCGATCGC",
  "reverse_primer": "GCAGCTAA",
  "max_product_length": 28,
  "product_count": 1,
  "products": [
    {
      "sequence": "ATGCGATCGCTTTTTTTTTTTTAGCTGC",
      "start": 0,
      "end": 28,
      "length": 28,
      "wraps_origin": false,
      "forward_binding": {"start": 0, "end": 10, "mismatches": 0},
      "reverse_binding": {"start": 20, "end": 28, "mismatches": 0}
    }
  ]
}
```

- `products[]` 按（`start`, `length`）**升序、去重**返回：
  - `sequence`：top 链 5'→3' 方向的产物序列（环形跨原点时拼接环首尾，长度恒为
    一圈以内）；
  - `start` / `end`：产物在 top 链上的 0-based 半开区间端点；环形产物跨原点时
    `end` 为取模后的坐标（可能小于 `start`），`wraps_origin` 为 `true`；
  - `length`：产物长度（含两端引物结合区）；
  - `forward_binding` / `reverse_binding`：两端引物结合区（反引物坐标为其反向
    互补序列在 top 链上的结合区）与各自错配数（恒为 0 或 1）。
- **无产物时 `products` 为空数组 `[]`**（背向位点、越界、结合区重叠、引物过长、
  全部位点均被 `max_product_length` 过滤等）。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| `max_product_length` 超过模板长度 | **422** |
| `max_product_length` 小于 1 / 非整数 | **422** |
| 引物为空或短于 2 个碱基 / 含非 IUPAC 字符 | **422** |
| 缺少必填字段 | **422** |
| 无符合条件的产物 | **200**，`products` 为 `[]` |

## 开放阅读框（ORF）检索

基于已保存记录的**原始 IUPAC 序列与线性/环形拓扑**（**不依赖酶切片段**，
即使记录没有任何切口也可使用），在**正链**与**反向互补链**上检索开放阅读框。
**检索过程不写库**：原记录保持不变，也不会产生新的提交。

服务地址与端口同其它接口（默认 **http://localhost:8000**，容器内 8000 端口
经 compose 映射到宿主机 8000）：

### 检索口径

- 在两条链上分别**逐碱基**寻找**字面量** `ATG` 作为起始密码子；含简并符号
  的密码子（如 `NTG`、`ATN`）**不充当起始信号**。
- 从每个起点按三联体逐个密码子读取，只有遇到**同框**的第一个**字面量**
  `TAA` / `TAG` / `TGA` 才形成 ORF；不同阅读框上的终止密码子忽略；含简并
  符号的密码子（如 `TAR`、`TGA` 中的简并变体）**不充当终止信号**，读取继续。
- **线性**：读取不得越出模板边界；边界前没有同框终止密码子的起点不产生 ORF。
- **环形**：读取可**跨原点**，但**完整 ORF（含终止密码子）不得超过模板一圈**；
  读取超过一圈仍未见终止密码子的起点不产生 ORF。
- 保留氨基酸数（密码子总数减 1，**不含终止**）不少于 `min_aa_length` 的
  **全部嵌套 ORF**——每个 ATG 各自独立寻找同框第一个终止密码子，同框内部
  ATG 给出的更短嵌套 ORF 全部保留。

### 请求

以下示例假设 `id=1` 的记录为线性序列 `ATGAAATTTTAG`（可先通过
`POST /digests` 创建；模板上无需任何酶切位点）：

```bash
curl -X POST http://localhost:8000/digests/1/orfs \
  -H 'Content-Type: application/json' \
  -d '{"min_aa_length": 2}'
```

- `min_aa_length`：**必填**正整数（`>= 1`），ORF 至少应含的氨基酸数
  （不含终止密码子，含起始甲硫氨酸）。

### 响应

上例唯一 ORF 为 `ATG AAA TTT TAG`（3 个氨基酸，核苷酸长 12）：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 12,
  "min_aa_length": 2,
  "orf_count": 1,
  "orfs": [
    {
      "strand": "+",
      "start": 0,
      "end": 12,
      "wraps_origin": false,
      "sequence": "ATGAAATTTTAG",
      "length": 12,
      "aa_length": 3
    }
  ]
}
```

- `orfs[]` 按（`strand`, `start`, `length`）**稳定升序**排列，正链 `+`
  排在反链 `-` 之前：
  - `strand`：链方向，`+` 为原始正链（top 链），`-` 为反向互补链；
  - `start` / `end`：ORF 在**原始正链**上覆盖区间的 0-based 半开端点；
    环形 ORF 跨原点时 `end` 为取模后的坐标（可能小于 `start`），
    `wraps_origin` 为 `true`；反链坐标同样给出其在正链上的覆盖区间
    （与 ORF 自身 5'->3' 方向相反）；
  - `sequence`：ORF 自身 **5'->3' 方向的定向核苷酸序列**（反链 ORF 为
    反向互补链上的序列，环形跨原点时拼接环首尾），含终止密码子；
  - `length`：核苷酸长度（恒为 3 的倍数）；`aa_length`：氨基酸数
    （`= length / 3 - 1`，不含终止）。
- **无结果时 `orfs` 为空数组 `[]`**（无字面量 ATG、线性越界、环形超过
  一圈、全部 ORF 短于 `min_aa_length` 等）。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| `min_aa_length` 小于 1 / 非整数 / 缺失 | **422** |
| 无符合条件的 ORF | **200**，`orfs` 为 `[]` |

## 凝胶条带比较

只读比较 **2~8 条已保存记录**的凝胶条带：把各记录（泳道）完整酶切结果的全部
片段长度汇总，按**长度降序**聚类为条带——新长度并入当前条带当且仅当它与条带
**最大长度**的相对跨度 `(max - min) / max <= resolution_ratio`。判定始终相对
条带最大值而非相邻元素，因此不会出现"相邻两两可合并、首尾却超差"的链式归并。
**分析过程不写库**：原记录保持不变，也不会产生新的提交。

### 请求

以下示例假设 `id=1`、`id=2` 为两条已保存记录（可先通过 `POST /digests` 创建）：

```bash
curl -X POST http://localhost:8000/gel-comparisons \
  -H 'Content-Type: application/json' \
  -d '{"digest_ids": [1, 2], "resolution_ratio": 0.05}'
```

- `digest_ids`：2~8 个已保存记录 id，不可重复；顺序即泳道顺序。
- `resolution_ratio`：分辨率比例，闭区间 **`[0, 0.2]`**；`0` 表示仅完全相同的
  长度才并入同一条带。

### 响应

```json
{
  "digest_ids": [1, 2],
  "resolution_ratio": 0.05,
  "lane_count": 2,
  "band_count": 3,
  "bands": [
    {"min_length": 60, "max_length": 60, "lane_fragment_counts": [1, 0],
     "fragment_count": 1, "total_bases": 60, "diagnostic": true},
    {"min_length": 40, "max_length": 41, "lane_fragment_counts": [1, 1],
     "fragment_count": 2, "total_bases": 81, "diagnostic": false},
    {"min_length": 20, "max_length": 21, "lane_fragment_counts": [0, 2],
     "fragment_count": 2, "total_bases": 41, "diagnostic": true}
  ]
}
```

- `bands[]` 按 `max_length` **降序**稳定排列（各条带的 max 必然互不相同，
  次序完全确定）：
  - `min_length` / `max_length`：条带覆盖的片段长度范围；
  - `lane_fragment_counts`：各泳道落入该条带的片段数，与请求 `digest_ids`
    的顺序一一对应；
  - `fragment_count`：条带内片段总数；`total_bases`：条带内片段长度之和；
  - `diagnostic`：仅出现在一个泳道的**诊断条带**为 `true`。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 任一记录不存在（按请求顺序报告第一个缺失 id） | **404** |
| `digest_ids` 重复 / 少于 2 条 / 多于 8 条 | **422** |
| `resolution_ratio` 不在 `[0, 0.2]`（含 NaN / 无穷大） | **422** |

## 末端连接兼容性分析

只读枚举 **2~8 条已保存记录**之间可连接的片段末端对：取每条记录完整酶切结果中
每个片段的左/右末端，枚举**跨记录**的末端组合（同一记录内的末端互不配对，同一
片段端更不会与自身配对），按兼容规则筛选后稳定返回。**分析过程不写库**：原记录
保持不变，也不会产生新的提交。

### 兼容规则

- **黏端 `sticky`**：一端为 5' 突出端、另一端为 3' 突出端，且 5' 突出端序列
  **恰好等于** 3' 突出端序列的反向互补（IUPAC 逐字符反向互补后的字符串相等，
  不做简并匹配）；5'-5'、3'-3' 即使突出端相同也不兼容。
- **平端 `blunt`**：仅与平端兼容，且仅当 `include_blunt` 为真（默认）时纳入结果；
  设为 `false` 时只返回黏端配对。
- 其余组合（含 `none` 天然末端、黏端-平端混合）一律不兼容。
- `ligated_overhang`（连接后突出端）：黏端对取该对中 **5' 突出端**的序列
  （即接头处退火的突出端），平端对为空字符串。

### 请求

以下示例创建两条记录：`id=1` 用 EcoRI（`GAATTC`，切口 1/5，5' 突出端 `AATT`），
`id=2` 用同一识别模式但切口 5/1 的酶（3' 突出端 `AATT`，其反向互补仍为 `AATT`）：

```bash
curl -X POST http://localhost:8000/ligation-compatibility \
  -H 'Content-Type: application/json' \
  -d '{"digest_ids": [1, 2], "include_blunt": true}'
```

- `digest_ids`：2~8 个已保存记录 id，不可重复；记录顺序决定配对的稳定排序，
  每对的 `first` 恒为请求顺序靠前记录的一端。
- `include_blunt`：**可选**，默认 `true`；为 `false` 时不返回平端-平端配对。

### 响应

两条记录各有 2 个黏端（片段 0 的右端、片段 1 的左端），产生 2 × 2 = 4 个配对，
按（first 记录、片段索引、端方向，second 同序）排列：

```json
{
  "digest_ids": [1, 2],
  "include_blunt": true,
  "record_count": 2,
  "pair_count": 4,
  "pairs": [
    {"first":  {"digest_id": 1, "fragment_index": 0, "end": "right",
                "end_type": "5prime", "overhang": "AATT"},
     "second": {"digest_id": 2, "fragment_index": 0, "end": "right",
                "end_type": "3prime", "overhang": "AATT"},
     "compatibility": "sticky", "ligated_overhang": "AATT"},
    {"first":  {"digest_id": 1, "fragment_index": 0, "end": "right",
                "end_type": "5prime", "overhang": "AATT"},
     "second": {"digest_id": 2, "fragment_index": 1, "end": "left",
                "end_type": "3prime", "overhang": "AATT"},
     "compatibility": "sticky", "ligated_overhang": "AATT"},
    {"first":  {"digest_id": 1, "fragment_index": 1, "end": "left",
                "end_type": "5prime", "overhang": "AATT"},
     "second": {"digest_id": 2, "fragment_index": 0, "end": "right",
                "end_type": "3prime", "overhang": "AATT"},
     "compatibility": "sticky", "ligated_overhang": "AATT"},
    {"first":  {"digest_id": 1, "fragment_index": 1, "end": "left",
                "end_type": "5prime", "overhang": "AATT"},
     "second": {"digest_id": 2, "fragment_index": 1, "end": "left",
                "end_type": "3prime", "overhang": "AATT"},
     "compatibility": "sticky", "ligated_overhang": "AATT"}
  ]
}
```

- `pairs[].first` / `second`：参与配对的两个片段末端，含 `digest_id`、
  `fragment_index`（与 `GET /digests/{id}` 结果中的片段索引一致）、
  `end`（`left` / `right`）、`end_type` 与该端的 `overhang`。
- **稳定排序**：全部末端先按（记录请求顺序、片段索引、端方向 left < right）
  排成全序，配对先按 `first` 再按 `second` 的同一键序返回——相同请求（无论
  `digest_ids` 传入顺序，只要内容相同且顺序确定）结果必然一致。
- 环形记录单切口开环产生的全长片段两端均为真实切割末端，同样参与枚举；
  无切口记录（两端为 `none`）不产生任何配对。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 任一记录不存在（按请求顺序报告第一个缺失 id） | **404** |
| `digest_ids` 重复 / 少于 2 条 / 多于 8 条 | **422** |
| `digest_ids` 含非整数、`include_blunt` 为非法布尔值 | **422** |

## 片段拼接链规划

只读规划 **2~8 条已保存记录**之间的最长线性片段拼接链（assembly chain）：
每条记录从其完整酶切结果中**至多贡献一个片段**参与拼接（一条记录在链中
最多出现一次）；片段是无方向元件，可以按原始方向（`forward`，`left→right`）
或翻转方向（`reverse`，`right→left`）摆入链中。服务把"记录 × 片段 ×
摆放方向"展开为**有向图节点**，节点间按 [末端连接兼容规则](#兼容规则)
连边（仅跨记录连边），在图上搜索**节点数最多的简单路径**。**规划过程不写库**：
原记录保持不变，也不会产生新的提交。

### 请求

```bash
curl -X POST http://localhost:8000/ligation-assembly \
  -H 'Content-Type: application/json' \
  -d '{"digest_ids": [1, 2, 3], "include_blunt": true, "max_chain_length": 8}'
```

- `digest_ids`：2~8 个已保存记录 id，不可重复；**记录的请求顺序参与同长度
  方案的稳定择优**（见下）。
- `include_blunt`：**可选**，默认 `true`；为 `false` 时平端不能作为连接末端，
  只允许黏端连接。
- `max_chain_length`：**必填**，拼接链允许包含的最大节点（片段）数，整数闭区间
  **`[2, 8]`**。每条记录最多贡献一个片段，因此链不可能长于记录数。

### 最长链与稳定择优

- 先最大化**节点数**（等价于最大化连接数 = 节点数 - 1）；
- 节点数相同的方案，按节点序列逐键比较
  **（记录请求位置 → 片段索引 → 链左端所用原片段端方向 `left < right`）**，
  取字典序最小者。因此：优先使用请求顺序靠前的记录；同一记录优先使用片段索引
  更小的片段；正向/翻转均可行时优先"链左端取自原片段左端"的摆放；
- 比较结果与内部枚举顺序无关，**相同请求（含 digest_ids 顺序）结果必然一致**；
  调换 `digest_ids` 顺序，选出的链也会相应改变。
- 搜索为带剪枝的记忆化全枚举（状态 = 当前出端 × 未用记录集合 × 剩余链长槽位），
  规模有界：全部记录的片段总数**超过 64 个**时返回 422。

### 响应

以下示例为三条记录：`id=1` 单切口产生 5' 突出端 `AAGT`（f0 右端）；`id=2`
含两个切口，中间片段 f1 左端为 5' `AAGT`、右端为 3' `ACTT`，可作桥接；
`id=3` 单切口产生 3' 突出端 `ACTT`（与 5' `AAGT` 反向互补）：

```json
{
  "digest_ids": [1, 2, 3],
  "include_blunt": true,
  "max_chain_length": 8,
  "record_count": 3,
  "node_count": 3,
  "connection_count": 2,
  "nodes": [
    {"digest_id": 1, "fragment_index": 0, "orientation": "forward",
     "left_end":  {"digest_id": 1, "fragment_index": 0, "end": "left",
                   "end_type": "none", "overhang": ""},
     "right_end": {"digest_id": 1, "fragment_index": 0, "end": "right",
                   "end_type": "5prime", "overhang": "AAGT"}},
    {"digest_id": 2, "fragment_index": 1, "orientation": "reverse",
     "left_end":  {"digest_id": 2, "fragment_index": 1, "end": "right",
                   "end_type": "3prime", "overhang": "ACTT"},
     "right_end": {"digest_id": 2, "fragment_index": 1, "end": "left",
                   "end_type": "5prime", "overhang": "AAGT"}},
    {"digest_id": 3, "fragment_index": 0, "orientation": "reverse",
     "left_end":  {"digest_id": 3, "fragment_index": 0, "end": "right",
                   "end_type": "3prime", "overhang": "ACTT"},
     "right_end": {"digest_id": 3, "fragment_index": 0, "end": "left",
                   "end_type": "none", "overhang": ""}}
  ],
  "connections": [
    {"index": 0,
     "from_end": {"digest_id": 1, "fragment_index": 0, "end": "right",
                  "end_type": "5prime", "overhang": "AAGT"},
     "to_end":   {"digest_id": 2, "fragment_index": 1, "end": "right",
                  "end_type": "3prime", "overhang": "ACTT"},
     "compatibility": "sticky", "ligated_overhang": "AAGT"},
    {"index": 1,
     "from_end": {"digest_id": 2, "fragment_index": 1, "end": "left",
                  "end_type": "5prime", "overhang": "AAGT"},
     "to_end":   {"digest_id": 3, "fragment_index": 0, "end": "right",
                  "end_type": "3prime", "overhang": "ACTT"},
     "compatibility": "sticky", "ligated_overhang": "AAGT"}
  ],
  "used_digest_ids": [1, 2, 3],
  "unused_digest_ids": [],
  "stop_reason": "all_records_used",
  "stop_reason_detail": "every submitted record is used in the chain; ..."
}
```

- `nodes[]`：链上节点，沿链从左到右排列：
  - `fragment_index` 与 `GET /digests/{id}` 结果中的片段索引一致；
  - `orientation` 为 `forward`（按原片段 left→right 摆入）或 `reverse`
    （翻转摆入）；
  - `left_end` / `right_end` 为该节点在**拼接链上**的左/右端点，其中 `end`
    标明该端取自**原片段**的 `left` 还是 `right`（翻转时二者相反）；天然末端
    `none` 只可能出现在整条链的两个最外端。
- `connections[]`：相邻节点间的连接，`index` 从 0 起；`from_end` 恒为左侧
  节点的出端、`to_end` 恒为右侧节点的入端，二者来自不同记录；
  `compatibility` / `ligated_overhang` 的含义与 `/ligation-compatibility`
  完全一致（黏端取接头处退火的 5' 突出端；平端为空字符串）。
- `used_digest_ids` / `unused_digest_ids`：入选 / 未入选链的记录 id，均按
  请求顺序排列。

### 未能延伸的原因（`stop_reason`）

| 值 | 含义 |
| --- | --- |
| `all_records_used` | 提交的全部记录都已在链中，没有剩余记录可继续延伸 |
| `max_chain_length_reached` | 链已达到 `max_chain_length`，且链右端与某条**未用**记录的候选末端仍兼容——继续延伸被上限禁止 |
| `no_compatible_extension` | 仍有未用记录，但链右端与其中任何候选末端都不兼容 |
| `no_compatible_ends` | 全部记录之间不存在任何一对兼容末端，无法形成含至少一处连接的链；此时 `nodes` / `connections` 为空 |

环形记录单切口开环产生的全长片段两端均为真实切割末端，两个方向都可参与
拼接；无切口记录（两端为 `none`）不会产生候选节点，永远不会入选。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 任一记录不存在（按请求顺序报告第一个缺失 id） | **404** |
| `digest_ids` 重复 / 少于 2 条 / 多于 8 条 | **422** |
| 缺少 `max_chain_length`，或其不在 `[2, 8]` / 非整数 | **422** |
| `digest_ids` 含非整数、`include_blunt` 为非法布尔值 | **422** |
| 全部记录的候选片段总数超过 64（枚举上限） | **422** |

## 固定窗口 GC 含量范围

对请求中**直接提交的 IUPAC DNA 序列**（不引用任何已保存记录），按固定长度窗口
计算每个窗口的 **GC 含量范围**——即最低与最高 GC 百分比。**计算过程不写库**：
不会产生新的酶切提交，也不读取数据库。

### 简并碱基如何决定上下界

逐个考察碱基的 IUPAC 可能碱基集合是否落在 `{G, C}` 中：

- **确定碱基**：`G`/`C` 对上下界贡献相同（恒计入），`A`/`T` 恒不计入；
- **简并符号**：
  - 可能集合**全部**为 G/C 时计入**下界**——简并码中只有 `S={G,C}` 必然计入；
  - 可能集合与 `{G,C}` **有交集**时计入**上界**——`R`(AG)、`Y`(CT)、
    `K`(GT)、`M`(AC)、`B`(CGT)、`D`(AGT)、`H`(ACT)、`V`(ACG)、`N`(ACGT)
    均可能计入，仅 `W`(AT) 不可能计入。

因此下界 = "窗口内必然为 GC 的碱基数 / 窗口长度"，
上界 = "窗口内可能为 GC 的碱基数 / 窗口长度"；序列不含简并碱基时
`min == max`，即为普通 GC 百分比。

### 窗口与滑动

- `window_length`：固定窗口长度，正整数，**不得超过序列长度**（否则 422）；
  窗口长度**等于**序列长度时，无论 `step` 为何，只返回唯一一个整长窗口。
- `step`：窗口长度**小于**序列长度时的滑动步长，正整数：
  - **线性序列**：只保留完整落在序列内的窗口，起点为
    `0, step, 2·step, …` 中所有满足 `start + W ≤ L` 者；末端可能剩余不足
    一步的短尾，不另补窗口，窗口总数为 `floor((L - W) / step) + 1`；
  - **环形序列**：起点为 `0, step, 2·step, …` 中所有 `< L` 者，共
    `ceil(L / step)` 个；最后一个窗口**跨越末尾回到开头**，按模取碱基。
    窗口长度恒不超过序列长度，因此**每个窗口最多覆盖一圈**，同一物理位置
    不会在一个窗口内被重复计入。

### 请求

```bash
curl -X POST http://localhost:8000/gc-content \
  -H 'Content-Type: application/json' \
  -d '{"sequence": "ATGCSWNATG", "topology": "circular",
       "window_length": 4, "step": 3}'
```

- `sequence`：仅允许 IUPAC 碱基 `ACGTRYSWKMBDHVN`（大小写均可，按大写处理）。
- `topology`：**可选**，`linear`（默认）或 `circular`。
- `window_length` / `step`：正整数，含义如上。

### 响应

以上例（`L=10`、`W=4`、`step=3`、环形）为例，窗口起点为 0、3、6、9；
起点 9 的窗口覆盖位置 9、0、1、2（`G,A,T,G`）并跨原点（起点 6 的窗口恰好
终止于位置 10=原点，不计为跨原点）：

```json
{
  "topology": "circular",
  "sequence_length": 10,
  "window_length": 4,
  "step": 3,
  "window_count": 4,
  "windows": [
    {"index": 0, "start": 0, "end": 4, "length": 4, "wraps_origin": false,
     "min_gc_count": 2, "max_gc_count": 2,
     "min_gc_percent": 50.0, "max_gc_percent": 50.0},
    {"index": 1, "start": 3, "end": 7, "length": 4, "wraps_origin": false,
     "min_gc_count": 2, "max_gc_count": 3,
     "min_gc_percent": 50.0, "max_gc_percent": 75.0},
    {"index": 2, "start": 6, "end": 0, "length": 4, "wraps_origin": false,
     "min_gc_count": 1, "max_gc_count": 2,
     "min_gc_percent": 25.0, "max_gc_percent": 50.0},
    {"index": 3, "start": 9, "end": 3, "length": 4, "wraps_origin": true,
     "min_gc_count": 2, "max_gc_count": 2,
     "min_gc_percent": 50.0, "max_gc_percent": 50.0}
  ]
}
```

- `windows[]` 按起点升序返回：
  - `start` / `end`：0-based 半开区间端点；环形跨原点窗口的 `end` 为
    `(start + W) mod L`（可能**小于** `start`），`wraps_origin` 为 `true`；
  - `min_gc_count` / `max_gc_count`：窗口内必然 / 可能为 GC 的碱基数，
    恒有 `min ≤ max`；
  - `min_gc_percent` / `max_gc_percent`：对应百分比（`100 × 计数 / 窗口长度`），
    四舍五入保留 **6 位小数**；需要精确值时可用计数除以 `length` 复原。
- 该接口与 `/digests` 体系完全独立：**无需先提交酶切记录**，请求中的序列
  也不会被保存。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| `window_length` 超过序列长度 | **422** |
| `window_length` 或 `step` 小于 1 / 非整数 | **422** |
| 序列为空，或含非 IUPAC 字符 | **422** |
| `topology` 不是 `linear` / `circular` | **422** |

## 未知样本识别（酶切片段长度配对）

只读从**未知样本的观测片段长度**（如凝胶上读得的条带长度）出发，对
**1~8 条已保存记录**逐一做**观测片段 ↔ 预期酶切片段**的一对一配对，按统一
分数排序并给出命中与歧义判定。**识别过程不写库**：原记录保持不变，也不会
产生新的提交。

### 配对与择优口径

- 观测列表与某条记录的预期片段列表各以**列表索引**（0-based）标识片段；
  **重复长度按不同片段计数**（不做多重集折叠）。
- 一个观测与一个预期**至多配对一次**（任一片段都不能被两个配对复用），
  配对当且仅当 **`|observed_length - expected_length| <= tolerance`**
  （整数 bp，闭区间）。
- 对每条候选记录，优化目标依次为：
  1. **最大化配对数**；
  2. 配对数相同则**最小化差值总和**（配对的长度差之和）；
  3. 仍相同时，把配对按（观测索引, 预期索引）排序得到配对序列，取
     **字典序最小**的方案：编号更小的观测优先被配出（未配对因此集中在
     编号最大的观测上），被配出的观测优先占用编号最小的可行预期。
- **未配对总数** = 观测未配对数 + 预期未配对数。
- 候选之间按 **（未配对总数, 差值总和, 请求中的 id 顺序）** 升序排序。
- `max_unpaired` 为命中门槛：未配对总数**不超过**该值的候选才算命中；
  无候选满足时返回**无命中**。命中候选中最优分数
  （未配对总数, 差值总和）被多个候选并列时标记**歧义**，`best_digest_id`
  取请求顺序最靠前的并列候选作为代表。
- 每侧至多 **12** 个片段（观测列表与单条记录的预期片段均受此限），
  超出返回 422。

### 请求

```bash
curl -X POST http://localhost:8000/digest-identification \
  -H 'Content-Type: application/json' \
  -d '{
    "digest_ids": [1, 2],
    "observed_fragment_lengths": [100, 200, 305],
    "tolerance": 10,
    "max_unpaired": 2
  }'
```

- `digest_ids`：1~8 个已保存记录 id，不可重复；顺序参与候选排序与歧义时
  代表 id 的选取。
- `observed_fragment_lengths`：未知样本的观测片段长度（bp），1~12 个
  **正整数**，重复长度按不同片段计数。
- `tolerance`：允许配对的最大长度差（bp，含），**非负整数**。
- `max_unpaired`：命中所允许的两侧未配对片段总数上限（含），**非负整数**。

### 响应

下例中 `id=1` 的记录片段为 `[100, 200, 300]`（3 对、差值 5、0 未配对，
命中且唯一最优），`id=2` 的记录片段为 `[100, 200, 400]`（2 对，观测 305
与预期 400 均未配对，共 2 个未配对但差值为 0）：

```json
{
  "digest_ids": [1, 2],
  "observed_fragment_lengths": [100, 200, 305],
  "tolerance": 10,
  "max_unpaired": 2,
  "candidate_count": 2,
  "hit": true,
  "ambiguous": false,
  "best_digest_id": 1,
  "candidates": [
    {
      "digest_id": 1,
      "fragment_count": 3,
      "pair_count": 3,
      "unpaired_count": 0,
      "diff_sum": 5,
      "eligible": true,
      "best": true,
      "pairs": [
        {"observed_index": 0, "expected_index": 0,
         "observed_length": 100, "expected_length": 100, "length_diff": 0},
        {"observed_index": 1, "expected_index": 1,
         "observed_length": 200, "expected_length": 200, "length_diff": 0},
        {"observed_index": 2, "expected_index": 2,
         "observed_length": 305, "expected_length": 300, "length_diff": 5}
      ],
      "unpaired_observed_indices": [],
      "unpaired_expected_indices": []
    },
    {
      "digest_id": 2,
      "fragment_count": 3,
      "pair_count": 2,
      "unpaired_count": 2,
      "diff_sum": 0,
      "eligible": true,
      "best": false,
      "pairs": [
        {"observed_index": 0, "expected_index": 0,
         "observed_length": 100, "expected_length": 100, "length_diff": 0},
        {"observed_index": 1, "expected_index": 1,
         "observed_length": 200, "expected_length": 200, "length_diff": 0}
      ],
      "unpaired_observed_indices": [2],
      "unpaired_expected_indices": [2]
    }
  ]
}
```

- `hit`：是否存在未配对总数不超过 `max_unpaired` 的候选；无命中时为
  `false`，同时 `ambiguous` 为 `false`、`best_digest_id` 为 `null`
  （候选仍全部回显，各自 `eligible` 为 `false`）。
- `ambiguous`：命中候选中最优分数（未配对总数, 差值总和）被多个候选并列时
  为 `true`；`best_digest_id` 取请求顺序最靠前的并列候选作为代表，每个
  并列候选的 `best` 均为 `true`。
- `candidates[]` 按 **（未配对总数, 差值总和, 请求 ID 顺序）** 稳定排序：
  - `pairs[]` 按观测索引升序，含两侧索引、长度与 `length_diff`
    （恒不超过容差）；
  - `unpaired_observed_indices` / `unpaired_expected_indices`：两侧未配对
    片段的索引，各自按索引升序；
  - `eligible`：未配对总数是否不超过 `max_unpaired`（计入命中）。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 任一记录不存在（按请求顺序报告第一个缺失 id） | **404** |
| `digest_ids` 重复 / 少于 1 条 / 多于 8 条 | **422** |
| `observed_fragment_lengths` 为空、超过 12 个、含非正整数 | **422** |
| `tolerance` 或 `max_unpaired` 为负数 / 非整数 / 缺失 | **422** |
| 某条记录的酶切片段数超过 12（按请求顺序报告第一条） | **422** |
| 无候选满足 `max_unpaired` | **200**，`hit` 为 `false` |

## 行为规则

- **双向搜索**：同时在 top 链搜索 `pattern` 与 `reverse_complement(pattern)`
  （即酶结合于反链的情形）；回文模式只搜索一次。
- **重叠位点**：允许识别位点相互重叠，每个命中都会计算切口。
- **简并碱基**：模式与序列中的 IUPAC 简并碱基按"可能集合有交集"匹配
  （如模式中的 `A` 可匹配序列中的 `R`）。
- **跨原点**：环形序列在扩展序列上搜索，跨越首尾的识别位点与切口均可命中。
- **物理切口去重**：不同酶/不同位点若产生相同的 `(top_cut, bottom_cut)`，
  只保留一条切口记录，酶名合并进 `enzymes` 列表。
- **非法偏移整体拒绝**：任一酶的切口偏移超出 `[0, len(pattern)]`（含负数），
  或序列/模式含非 IUPAC 字符，整个请求以 **422** 拒绝，**不会有任何数据落库**。
- **无切口**：返回一条覆盖整条分子的片段，两端类型为 `none`。

## 项目结构

```
app/
├── main.py      # FastAPI 入口与路由
├── digest.py    # 酶切核心算法（搜索、去重、片段生成）
├── variant.py   # 变异施加与切口/片段差异对比
├── panel.py     # 酶组合面板规划（枚举候选组合、最少酶择优）
├── distribution.py  # 部分酶切概率分布（枚举成功/失败状态、长度多重集合并）
├── methylation.py   # 甲基化阻断分析（敏感位点定位、切口保留/移除、片段重建）
├── pcr.py        # PCR 扩增产物预测（IUPAC 引物匹配、错配/3'端规则、相向位点对枚举）
├── orf.py        # 开放阅读框检索（双链字面量 ATG/终止密码子扫描、线性越界/环形一圈口径）
├── gel.py       # 凝胶条带比较（片段长度降序聚类、诊断条带标记）
├── ligation.py  # 末端连接兼容性分析（跨记录末端对枚举、黏端/平端兼容判定）
├── assembly.py  # 片段拼接链规划（有向图、不重复用记录的最长简单路径与稳定择优）
├── gccontent.py # 固定窗口 GC 含量范围（简并碱基上下界、线性/环形滑窗）
├── identification.py # 未知样本识别（观测/预期片段一对一最优配对、字典序择优、命中/歧义判定）
├── iupac.py     # IUPAC 编码表与反向互补
├── schemas.py   # 请求/响应模型与校验
└── db.py        # SQLite 持久化（启动时自动建表）
tests/
├── test_digest.py       # 酶切算法测试
├── test_api.py          # 提交/查询/校验 API 测试
├── test_variant_api.py  # 变异影响分析 API 测试
├── test_panel_api.py    # 酶组合面板规划 API 测试
├── test_distribution_api.py  # 部分酶切概率分布 API 测试
├── test_methylation_api.py   # 甲基化阻断分析 API 测试
├── test_pcr_api.py       # PCR 扩增产物预测 API 测试
├── test_orf_api.py       # 开放阅读框检索 API 测试
├── test_gel_api.py       # 凝胶条带比较 API 测试
├── test_ligation_api.py  # 末端连接兼容性 API 测试
├── test_assembly_api.py  # 片段拼接链规划 API 测试
├── test_gc_content_api.py  # 固定窗口 GC 含量范围 API 测试
└── test_identification_api.py  # 未知样本识别 API 测试
Dockerfile
docker-compose.yml
requirements.txt
requirements-dev.txt
```

## 运行测试

本地运行：

```bash
pip install -r requirements-dev.txt
pytest
```

容器内运行（无需本地 Python 环境）：

```bash
docker build -t digest-api .
docker run --rm digest-api sh -c "pip install --no-cache-dir -r requirements-dev.txt && pytest"
```
