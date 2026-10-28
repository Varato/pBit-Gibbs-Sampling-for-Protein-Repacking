# CLAUDE.md

给在此仓工作的 coding agent 的上下文。人读详情看 README.md，本文件讲
"改代码前必须知道的事"。

## 项目是什么

用 MTJ（磁隧道结）p-bit 阵列 + FPGA 作为**物理随机源**做 Gibbs 采样的
科研演示仓：1L2Y（Trp-cage）侧链 repacking，13 残基，每 sweep 13 次
条件分布采样。科研想法：Gibbs 每步唯一的随机操作 = 从 Categorical(p)
抽一个类别——把这一步交给物理器件，物理熵源就替换了 PRNG。

对硬件合作方的全部要求 = 一个函数 `k = categorical_sample(p)`。
温度/能量/Boltzmann 折算全在引擎侧（`engine.cond_probs`），硬件方拿到
的 p 已是归一概率向量，不需要懂 Gibbs。

前身是 LTSpice 仿真版（`~/Workspace/Gibbs-Sampling-With-Memristors`，
接口是 `roll(v, kT)`）。本仓是其重构：问题实例与数据同口径（可交叉
验证），LTSpice 相关内容本仓一律没有，别往回加。

## 常用命令

```bash
# .venv 已建好（Python 3.14）；重建：
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt    # Linux/macOS 用 .venv/bin/

.venv/Scripts/python selftest.py                 # 独立自检，应 6/6 PASS
.venv/Scripts/python -m pytest tests -q          # 应 10 项全 PASS（含 GUI 冒烟）
.venv/Scripts/python check_sampler.py --backend cpu   # 验收协议，4/4 PASS
.venv/Scripts/python main.py                     # 打开 GUI
```

**任何改动后 selftest + pytest 必须保持全绿。** GMEC 复现（selftest B 节）
是金标准：CPU 后端几何退火 3→0.02、2000 sweep 必须命中 E=−8.407 且组合
与 `data/gibbs_ref.json` 逐残基一致——改引擎/采样逻辑破坏这一点 = 改错了。

Windows 控制台是 GBK，Python 打印中文乱码但不影响判定；要看清输出加
`PYTHONIOENCODING=utf-8` 前缀。

## 架构（数据流）

```
engine.GibbsEngine.run_sweep()          # 纯算法，零 GUI 依赖
  → table.cond(state, i)                # 条件能量 v（REU）
  → cond_probs(v, kT)                   # p = softmax(-v/kT) ← 温度在这里折算
  → sampler.categorical_sample(p)       # ★ 唯一随机原语（可插拔）
  → Frame 快照 → worker.GibbsWorker（线程）→ FrameBuffer → gui._poll 回放
```

| 模块 | 职责 |
|---|---|
| `sampler/base.py` | ★ 接口契约：CategoricalSampler / HardwareUnavailable / clip_dyn_range |
| `sampler/cpu.py` | CPUSampler（Gumbel-max 参考实现，对拍基准） |
| `sampler/mtj.py` | MTJ 硬件桩（**未实现**，构造即抛 HardwareUnavailable，2 个 TODO） |
| `engine.py` | GibbsTable 加载 + Gibbs 引擎 + kT 日程（恒温/几何退火） |
| `worker.py` | 采样线程（采集/回放分离；FrameBuffer min-max 对折保峰谷） |
| `gui.py` | tkinter 主窗口；worker 绝不碰 tk 对象，只写 queue/FrameBuffer |
| `canvaspanel.py` | 13 行骰子画布（itemconfig 差异更新，绝不删除重建） |
| `plotting.py` | matplotlib 内嵌能量面板（draw_idle 节流 ≤12 fps） |
| `data/` | 自包含 v1 截断池（152 候选，K=[11,4,5,2,34,2,6,42,3,3,3,34,3]），无需 PyRosetta |

## 接口契约（改动前必读）

`categorical_sample(p)` 是给外部硬件合作方的承诺，定义在 `sampler/base.py`，
**不要轻改**（改了要同步对方 + 验收协议）：

- 输入 p：一维非负、Σp≈1（±1e-6 容差）、长度 K 每次调用可变（本问题 2..42）
- 返回 Python int k∈0..K-1，统计上服从 Categorical(p)——这是
  `check_sampler.py` 的核心判据（TV 不劣于同 N 的 30 组 CPU 基准最大值）
- 错误语义：`RuntimeError` = 单次失败（引擎重试 ≤3 次后该残基保持原态
  继续）；`HardwareUnavailable` = 未实现/断连（引擎**不重试**立即上抛，
  GUI 弹提示框）——两者不能混
- 可选类属性：`blocking`（一次调用秒级以上→True，引擎上报逐残基进度、
  GUI 高亮"硬件采样中"行）；`dyn_range`（器件概率分辨率 p_max/p_min
  有限才填，引擎预先 clip_dyn_range 截断归一）

## 硬件对接现状

`MTJSampler` 是桩：合作方有成熟的 MTJ + FPGA 任意 Categorical 采样方案，
只需完成 `sampler/mtj.py` 的 TODO①（类属性）与 TODO②（建链 + 采样
实现），`python check_sampler.py --backend mtj --n 100` 全 PASS 即自动
可在 GUI 使用。硬件未接入时 GUI 弹提示、可切回 CPU——这条路径有测试
保护（selftest C / 无头 monkeypatch 验证过），别破坏。

**桩"构造即抛 HardwareUnavailable"是特性不是 bug**，工厂（make_sampler）、
GUI、selftest 都依赖这个语义。

## 约定与坑

- 注释、docstring、GUI 文案、提交信息全部中文（沿旧仓风格）；标识符英文
- Windows + Git Bash，路径用正斜杠；不要提交 `out/`、`.venv/`（已 gitignore）
- Frame 字段是 `n_samples`/`n_fail`（不是旧仓的 `n_rolls`）
- `gui._poll` 是 after(33ms) 自循环，退出必须走 `_on_close`（先 cancel
  定时器再 destroy——旧仓踩过销毁后仍触发的竞态）
- Tk 销毁后 `winfo_exists()` 抛 TclError 而非返回 0（Py3.14），测试里要
  try/except（见 `tests/test_gui_smoke.py`）
- ttk 控件 `['state']` 返回 index 对象不是字符串，别做 `== 'normal'` 比较
