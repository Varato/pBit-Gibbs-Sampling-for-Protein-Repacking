"""engine 无 GUI 依赖的核心测试（pytest）。

  数据     表加载：13 残基、v1 截断池 152 候选、K 与库一致、native 标记
           在界内、参照 GMEC 在位
  恒温冒烟  100 sweep 状态合法、能量有限、Frame 字段齐
  退火命中  CPU 采样器几何退火 3→0.02 命中记录 GMEC E=-8.407（同组合）
  容错     HardwareUnavailable 立即上抛（不重试不吞掉）

Run:  python -m pytest tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine import GibbsEngine, kT_const, kT_geometric, load_table  # noqa: E402
from sampler import CategoricalSampler, CPUSampler, HardwareUnavailable  # noqa: E402


@pytest.fixture(scope="module")
def table():
    return load_table()


def test_table_shape(table):
    assert table.N == 13
    assert table.K.sum() == 152 and table.K.max() == 42
    # library_chis（v1）与能量表同序（χ 解码安全性的前提）
    assert all(len(c) == k for c, k in zip(table.CHIS, table.K))
    assert all(0 <= nk < k for nk, k in zip(table.native_k, table.K))
    assert table.ref["gmec_E"] == pytest.approx(-8.407, abs=1e-3)


def test_const_smoke(table):
    eng = GibbsEngine(table, CPUSampler(seed=0), kT_const(1.0))
    for _ in range(100):
        f = eng.run_sweep()
    assert all(0 <= s < k for s, k in zip(eng.state, table.K))
    assert np.isfinite(f.E)
    assert f.changed >> 0 >= 0 and f.n_samples >= 0


def test_anneal_hits_gmec(table):
    eng = GibbsEngine(table, CPUSampler(seed=0), kT_geometric(n_anneal=2000))
    for _ in range(2000):
        eng.run_sweep()
    assert eng.best_E == pytest.approx(table.ref["gmec_E"], abs=1e-3)
    assert eng.best_state == table.ref["gmec_combo"]


class _Unreachable(CategoricalSampler):
    """模拟断连：采样中途抛 HardwareUnavailable。"""

    def __init__(self):
        self.n = 0

    def categorical_sample(self, p):
        self.n += 1
        raise HardwareUnavailable("模拟断连")


def test_hardware_unavailable_propagates(table):
    eng = GibbsEngine(table, _Unreachable(), kT_const(1.0))
    with pytest.raises(HardwareUnavailable):
        eng.run_sweep()
    assert eng.sampler.n == 1          # 没有重试 3 次
