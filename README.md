# MTJ Gibbs 采样（mtj-gibbs）

用 **MTJ（磁隧道结）p-bit 阵列 + FPGA** 作为物理随机源做 Gibbs 采样：
1L2Y（Trp-cage）侧链 repacking，13 个残基 = 13 次/步的条件分布采样，
随机源可插拔——CPU 参考实现，或合作方的 MTJ 硬件。

**对硬件的全部要求 = 一个函数**：

```python
k = categorical_sample(p)     # p: 离散概率向量 → k: 抽中的下标
```

温度 / 能量 / Boltzmann 折算全部在引擎侧完成，硬件方拿到的 p 已是
归一概率向量，不需要知道 Gibbs 是什么。

**当前口径 = v1 截断池**（packer 剪枝，152 候选，K = 11/4/5/2/34/2/6/
42/3/3/3/34/3，GMEC E = −8.407）。数据自包含（`data/`），不需要
PyRosetta。与 LTSpice 版（旧仓 `Gibbs-Sampling-With-Memristors`）同口径
同数据，可交叉验证。

## 快速开始

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt    # Linux/macOS 用 .venv/bin/pip
.venv/Scripts/python selftest.py                 # 独立复现自检（应全 PASS）
.venv/Scripts/python main.py                     # 打开 GUI
```

GUI 默认"退火"模式 + CPU 后端，点「开始采样」：13 行骰子实时翻转（灰带
= 本 sweep 翻过的行，角上三角 = native 最近候选），能量曲线收敛到
GMEC −8.41（蓝虚线），右轴橙线 = kT 降温日程。

选「硬件（MTJ + FPGA）」再点开始：**未实现/设备未连时弹提示框**，链路
不中断，可改回 CPU 继续。

## 核心接口（合作方只看这两个文件）

| 文件 | 内容 |
|---|---|
| `sampler/base.py` | **契约全文**：`CategoricalSampler` 抽象类、`HardwareUnavailable`、`clip_dyn_range` |
| `sampler/mtj.py` | **硬件桩**：2 个 TODO（类属性 + 建链/采样），接入时整段替换 |

契约速览（细节以 `base.py` 注释为准）：

```python
class CategoricalSampler(ABC):
    blocking  = False   # 一次调用秒级以上 → True（引擎据此上报逐残基进度）
    dyn_range = None    # 概率动态范围上限 p_max/p_min；器件分辨率有限才填，
                        # 引擎会预先截断归一（截断偏置由验收协议量化）

    def categorical_sample(self, p) -> int:
        """p: 一维非负数组 Σp≈1，长度 K 每次可变（本问题 2..42）
           返回 int k ∈ 0..K-1，统计上服从 Categorical(p)"""
```

错误语义：单次失败（超时/校验错）抛 `RuntimeError` → 引擎重试 ≤3 次，
仍败该残基保持原态继续；`HardwareUnavailable`（未实现/断连）→ 立即终止，
GUI 弹提示。

## 文件一览

| 文件 | 作用 |
|---|---|
| `main.py` | GUI 入口 |
| `gui.py` | 主窗口（控制面 / 回放 / 状态栏；后端固定 CPU \| 硬件两项） |
| `engine.py` | Gibbs 引擎（条件能量 → softmax → `categorical_sample(p)`） |
| `worker.py` | 采样线程 + 帧缓冲（采集/回放分离，min-max 对折保峰谷） |
| `sampler/base.py` | ★ 采样器接口契约（合作方看这里） |
| `sampler/cpu.py` | CPU 参考实现（Gumbel-max，数学严格等价） |
| `sampler/mtj.py` | ★ MTJ + FPGA 硬件桩（合作方 2 个 TODO） |
| `canvaspanel.py` | 13 行骰子阵列画布 |
| `plotting.py` | 能量轨迹面板（matplotlib 内嵌） |
| `check_sampler.py` | 采样器验收协议（对拍 + 噪声地板） |
| `selftest.py` | 独立复现自检（数据完整性 / GMEC / 接口冒烟） |
| `tests/` | pytest（引擎核心 + 采样器契约 + GUI 冒烟） |

## 硬件对接步骤（sampler/mtj.py，搜 "TODO" 共 2 处）

1. `dyn_range` / `blocking` 按器件实际填（分辨率有限→截断；往返秒级→
   blocking，GUI 会显示"硬件采样中"高亮行）
2. `__init__` 建立与 FPGA 的连接（以太网 / 串口 / 厂商 SDK…），连不上抛
   `HardwareUnavailable`；`categorical_sample` 里完成
   "写入 p → 触发采样 → 读回获胜下标"
3. K 每次不同（2..42）：设备侧支持变长输入，或一次配置最大 K、尾部补 0
   （补 0 的面永不获胜即可）
4. 验收：`python check_sampler.py --backend mtj --n 100` 全 PASS 再上 GUI

## 验收（本机已全绿）

- `selftest.py` 全 PASS（数据完整性 / GMEC 复现 / CPU uniform TV /
  mtj 桩按契约抛 HardwareUnavailable）
- `python -m pytest tests -q` 全 PASS
- `python check_sampler.py --backend cpu` 全 PASS
