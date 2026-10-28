"""Gibbs 引擎：条件能量 → 概率向量 → categorical_sample（零 GUI 依赖）。

与 LTSpice 版（旧仓 Gibbs-Sampling-With-Memristors）的关键区别：随机
原语的接口从 roll(v, kT) 降维为 categorical_sample(p)——温度折算收进
引擎（cond_probs），硬件方只面对归一概率向量。口径不变 = **v1 截断池**
（packer 剪枝，152 候选，K = [11,4,5,2,34,2,6,42,3,3,3,34,3]）。

  - 表/库/native 一次加载成 GibbsTable（χ 解码与 native 标记在内）
  - kT 日程函数化：恒温 / 几何退火（到 t_lo 后 clamp 恒温续跑，防
    kT→0 时 z=-v/kT 下溢成 NaN）
  - 单次采样失败容错：RuntimeError 重试 ≤3 次，仍败则该残基保持原态
    继续（演示不中断，误差记入 errors）；HardwareUnavailable（未实现/
    断连）不重试、立即上抛终止

Run:  被 worker / tests / selftest 导入，不直接运行
"""
import json
from dataclasses import dataclass, field
from typing import Callable, List

import numpy as np

from sampler import CategoricalSampler, HardwareUnavailable, clip_dyn_range

DATA = "data"

# χ 与 χ+180 同一构象（末端基团二重对称）→ 距离折叠（key = 残基一字母）
DEGEN = {"L": {2}, "F": {2}, "Y": {2}, "W": {2},
         "N": {2}, "D": {2}, "Q": {3}, "E": {3}, "R": {5}}


def circ_dist(a, b, fold):
    d = abs(((a - b + 180) % 360) - 180)
    return min(d, 180 - d) if fold else d


def cond_probs(v, kT):
    """条件能量向量 → 归一概率 p = softmax(-v/kT)。

    Gibbs 每步的全部物理：v 越小越有利，kT 是温度。硬件拿到的 p 已经
    折算完温度，不要再乘 kT。
    """
    z = -np.asarray(v, dtype=float) / kT
    z -= z.max()                     # 数值稳定
    w = np.exp(z)
    return w / w.sum()


@dataclass
class GibbsTable:
    """问题实例本体 + χ 解码 + native 标记 + 参照结果（一次加载）。"""
    res: List[str]
    K: np.ndarray
    ESELF: List[np.ndarray]
    adj: List[List]                       # adj[i] = [(j, M[K_i,K_j]), ...]
    CHIS: List[List[List[float]]]         # CHIS[i][k] = 候选 k 的 χ 向量
    native_k: List[int]                   # 最接近 native 的候选下标
    native_dev: List[float]               # 其最差 χ 偏差（°）
    native_chis: List[List[float]]
    ref: dict = field(default_factory=dict)   # gibbs_ref.json 参照（GMEC 等）

    def __post_init__(self):
        self.K = np.asarray(self.K, dtype=int)
        self.N = len(self.res)

    def cond(self, state, i):
        """残基 i 的条件能量向量 v（= categorical_sample 输入的来源）。"""
        v = self.ESELF[i].copy()
        for j, M in self.adj[i]:
            v += M[:, state[j]]
        return v

    def energy(self, state):
        e = sum(self.ESELF[i][state[i]] for i in range(self.N))
        for i in range(self.N):
            for j, M in self.adj[i]:
                if i < j:
                    e += M[state[i], state[j]]
        return e

    def chis_of(self, state):
        """下标组合 → 各残基 χ 角（度，Dunbrack 峰心）。"""
        return {self.res[i]: self.CHIS[i][state[i]] for i in range(self.N)}

    def cold_start(self):
        """各残基取 E_self 最优（冷启动）。"""
        return [int(np.argmin(v)) for v in self.ESELF]


def load_table(data_dir=DATA) -> GibbsTable:
    """加载 data/ 四件套（v1 截断池口径，152 候选）并预计算 native 标记。"""
    T = json.load(open(f"{data_dir}/energy_table.json"))
    L = {f"{r['aa1']}{r['pdbnum']}": r for r in json.load(
        open(f"{data_dir}/library_chis.json"))["residues"]}
    NAT = json.load(open(f"{data_dir}/native_chis.json"))["residues"]

    res = [n["res"] for n in T["nodes"]]
    idx = {r: i for i, r in enumerate(res)}
    K = [len(n["E_self"]) for n in T["nodes"]]
    ESELF = [np.array(n["E_self"], dtype=float) for n in T["nodes"]]
    adj = [[] for _ in res]
    for e in T["edges"]:
        i, j = idx[e["res_i"]], idx[e["res_j"]]
        M = np.array(e["E_int"], dtype=float)
        adj[i].append((j, M))
        adj[j].append((i, M.T))

    # native 最近候选（对当前池重算；χ 前缀对齐，proton χ 忽略）
    CHIS, native_k, native_dev, native_chis = [], [], [], []
    for r in res:
        chis = L[r]["rotamer_chis"]
        nat = NAT[r]
        degen = DEGEN.get(r[0], set())
        best = (0, 1e9)
        for k, rc in enumerate(chis):
            n = min(len(rc), len(nat))
            worst = max((circ_dist(rc[c], nat[c], (c + 1) in degen)
                         for c in range(n)), default=0.0)
            if worst < best[1]:
                best = (k, worst)
        CHIS.append(chis)
        native_k.append(best[0])
        native_dev.append(round(best[1], 1))
        native_chis.append(nat)

    ref = {}
    try:
        g = json.load(open(f"{data_dir}/gibbs_ref.json"))
        ref = {"gmec_E": g["gmec_E"],                      # v1: -8.407
               "gmec_combo": g["gmec_combo"],
               "marginal": g["marginal"]}                   # kT=0.5 参照链
    except FileNotFoundError:
        pass

    return GibbsTable(res=res, K=K, ESELF=ESELF, adj=adj, CHIS=CHIS,
                      native_k=native_k, native_dev=native_dev,
                      native_chis=native_chis, ref=ref)


# ---- kT 日程 ----------------------------------------------------------------
def kT_const(c):
    """恒温。"""
    return lambda s: c


def kT_geometric(t_hi=3.0, t_lo=0.02, n_anneal=2000):
    """几何降温 t_hi → t_lo，到 t_lo 后 clamp 恒温续跑（防下溢）。"""
    def f(s):
        return t_hi * (t_lo / t_hi) ** (min(s, n_anneal) / n_anneal)
    return f


# ---- 引擎 -------------------------------------------------------------------
class SweepInterrupted(Exception):
    """on_residue 回调抛出以提前结束本轮 sweep（暂停/停止响应粒度 =
    一次采样，硬件 blocking 模式下最多等一个残基的往返）。"""


@dataclass
class Frame:
    """一次 sweep 的快照（GUI 回放的最小单元，~29 B）。"""
    sweep_idx: int
    E: float
    kT: float
    state: np.ndarray        # (13,) rotamer 下标
    changed: int             # 位掩码：本 sweep 哪些行翻转
    best_E: float
    n_samples: int
    n_fail: int


class GibbsEngine:
    """逐 sweep 推进；UI/worker 只读快照，不碰内部。"""

    RETRY = 3                # 单次采样失败的重试上限

    def __init__(self, table: GibbsTable, sampler: CategoricalSampler,
                 kT_schedule=None, seed=0):
        self.table = table
        self.sampler = sampler
        self.kT_schedule = kT_schedule or kT_const(1.0)
        self.state = table.cold_start()
        self.n_sweeps = 0
        self.best_state = list(self.state)
        self.best_E = table.energy(self.state)
        self.errors: List[str] = []
        self._seed = seed

    def set_kT_schedule(self, f: Callable[[int], float]):
        """软更新：下一 sweep 生效（不重置轨迹）。"""
        self.kT_schedule = f

    def run_sweep(self, on_residue=None) -> Frame:
        """一轮系统扫描。on_residue(i, sweep, kT) 在每个残基采样前回调
        （blocking 后端用于上报进度；回调可抛 SweepInterrupted 提前结束
        本轮，用于暂停/停止）。采样失败重试 ≤RETRY 次后该残基保持原态
        继续。"""
        changed = 0
        kT = self.kT_schedule(self.n_sweeps)
        try:
            for i in range(self.table.N):
                if on_residue is not None:
                    on_residue(i, self.n_sweeps, kT)
                p = cond_probs(self.table.cond(self.state, i), kT)
                if self.sampler.dyn_range:
                    p = clip_dyn_range(p, self.sampler.dyn_range)
                k = None
                for attempt in range(self.RETRY):
                    try:
                        k = self.sampler.categorical_sample(p)
                        break
                    except HardwareUnavailable:
                        raise               # 未实现/断连：立即终止，不重试
                    except RuntimeError as err:
                        if attempt == self.RETRY - 1:
                            self.errors.append(f"sweep {self.n_sweeps} 残基 "
                                               f"{self.table.res[i]}: {err}")
                if k is not None and k != self.state[i]:
                    changed |= 1 << i
                    self.state[i] = k
        except SweepInterrupted:
            pass            # 部分完成的 sweep 也出帧（changed 已累计）
        self.n_sweeps += 1
        e = self.table.energy(self.state)
        if e < self.best_E:
            self.best_E, self.best_state = e, list(self.state)
        st = self.sampler.stats()
        return Frame(self.n_sweeps, round(e, 4), kT,
                     np.array(self.state, dtype=np.int16), changed,
                     round(self.best_E, 4), st.get("n_samples", 0),
                     st.get("n_fail", 0))

    def close(self):
        self.sampler.close()
