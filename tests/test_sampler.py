"""sampler 包契约测试（pytest）。

  CPU 分布正确性   固定偏斜 p 抽 4 万次，TV < 0.02；p=0 的面永不出现；
                   返回值在界内且为 int
  工厂             make_sampler("cpu") 可用；make_sampler("mtj") 按契约
                   抛 HardwareUnavailable（硬件未实现）
  dyn_range 截断   clip_dyn_range 抬尾 + 重归一

Run:  python -m pytest tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sampler import (CPUSampler, HardwareUnavailable,  # noqa: E402
                     clip_dyn_range, make_sampler)


def test_cpu_matches_categorical():
    p = np.array([0.5, 0.3, 0.15, 0.05, 0.0, 0.0])
    s = CPUSampler(seed=1)
    n, K = 40000, len(p)
    cnt = np.zeros(K)
    for _ in range(n):
        k = s.categorical_sample(p)
        assert isinstance(k, int) and 0 <= k < K
        cnt[k] += 1
    f = cnt / n
    tv = 0.5 * float(np.abs(f - p).sum())
    assert tv < 0.02, f"TV={tv}"
    assert cnt[4] == 0 and cnt[5] == 0     # p=0 的面永不出现


def test_cpu_defensive_normalization():
    s = CPUSampler(seed=2)
    p = np.array([2.0, 1.0, 1.0])          # 未归一输入也能正确处理
    ks = {s.categorical_sample(p) for _ in range(200)}
    assert ks <= {0, 1, 2} and len(ks) >= 2


def test_factory_cpu():
    s = make_sampler("cpu", seed=0)
    assert isinstance(s, CPUSampler)
    k = s.categorical_sample([0.5, 0.5])
    assert k in (0, 1)
    s.close()                              # 幂等 close 不炸


def test_factory_mtj_unimplemented():
    with pytest.raises(HardwareUnavailable):
        make_sampler("mtj")


def test_clip_dyn_range():
    p = np.array([1e-40, 1.0 - 1e-10, 1e-10])
    q = clip_dyn_range(p, 1e3)
    assert q.sum() == pytest.approx(1.0)
    assert q.min() >= q.max() / 1e3 - 1e-15   # 尾部被抬到 p_max/1e3
