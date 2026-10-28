"""独立复现自检：新机器上一条命令验证链路完整。

在干净 clone + venv + pip install -r requirements.txt 的机器上运行本
脚本，全 PASS 即证明：本仓自包含（无 PyRosetta / LTSpice 依赖），采样
侧复现旧仓记录结果（v1 截断池口径，GMEC 锚点取自入库 gibbs_ref.json，
其 GMEC 来自 phase1 qubo 的 Potts 退火 5 种子一致解）。

  A 数据完整性   表/库/native 三件套在位且同序，参照 GMEC 在位
  B GMEC 复现    CPU 采样器退火命中 E=-8.407 且组合与记录逐残基一致
  C 接口冒烟     CPU uniform TV < 0.03；mtj 工厂按契约抛 HardwareUnavailable
                 （硬件未实现时的预期行为）

Run:  python selftest.py
Out:  终端报告 + 退出码 0/1
"""
import sys

import numpy as np

from engine import GibbsEngine, kT_geometric, load_table
from sampler import CPUSampler, HardwareUnavailable, make_sampler

checks = []


def check(name, ok, detail=""):
    checks.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


print("A 数据完整性")
t = load_table()
ok = (t.N == 13 and int(t.K.sum()) == 152
      and all(len(c) == k for c, k in zip(t.CHIS, t.K))
      and all(0 <= nk < k for nk, k in zip(t.native_k, t.K)))
check("13 残基 / ΣK=152 / χ 库同序 / native 标记在界内", ok)
check("参照 GMEC 在位（-8.407）", abs(t.ref["gmec_E"] + 8.407) < 1e-3,
      f"E={t.ref['gmec_E']}")

print("B GMEC 复现（CPU 采样器退火 3→0.02，2000 sweep）")
eng = GibbsEngine(t, CPUSampler(seed=0), kT_geometric(n_anneal=2000))
for _ in range(2000):
    eng.run_sweep()
check("退火 best_E == -8.407", abs(eng.best_E - t.ref["gmec_E"]) < 1e-3,
      f"{eng.best_E:.3f} vs 记录 {t.ref['gmec_E']:.3f}")
check("退火组合逐残基一致", eng.best_state == t.ref["gmec_combo"],
      "一致" if eng.best_state == t.ref["gmec_combo"] else
      f"差异 {[(t.res[i], t.ref['gmec_combo'][i], eng.best_state[i]) for i in range(t.N) if t.ref['gmec_combo'][i] != eng.best_state[i]]}")

print("C 接口冒烟")
KU = 11
p0 = np.full(KU, 1.0 / KU)
d = CPUSampler(seed=3)
cnt = np.zeros(KU)
NU = 20000
for _ in range(NU):
    cnt[d.categorical_sample(p0)] += 1
tvu = 0.5 * float(np.abs(cnt / NU - p0).sum())
check("CPU uniform TV < 0.03", tvu < 0.03, f"TV={tvu:.4f}")
try:
    make_sampler("mtj")
    check("mtj 未实现 → HardwareUnavailable", False, "没有抛异常")
except HardwareUnavailable as err:
    check("mtj 未实现 → HardwareUnavailable", True, str(err)[:60] + "…")

n_pass = sum(checks)
print(f"\n{n_pass}/{len(checks)} PASS"
      + (" —— 独立复现成立" if n_pass == len(checks) else " —— 有 FAIL，先排查"))
sys.exit(0 if n_pass == len(checks) else 1)
