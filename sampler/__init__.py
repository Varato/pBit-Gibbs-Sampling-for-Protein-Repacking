"""采样器包：接口契约（base）+ CPU 参考（cpu）+ MTJ 硬件桩（mtj）。

常用入口：

    from sampler import CategoricalSampler, CPUSampler, make_sampler

Run:  被导入，不直接运行
"""
from typing import Literal

from .base import (CategoricalSampler, HardwareUnavailable,  # noqa: F401
                   SamplerStats, clip_dyn_range)
from .cpu import CPUSampler  # noqa: F401

Backend = Literal["cpu", "mtj"]


def make_sampler(backend: Backend, seed: int = 0, **kw) -> CategoricalSampler:
    """名字 → 采样器实例 的统一工厂（GUI / check_sampler / selftest 共用）。

    backend ∈ {cpu, mtj}；mtj 未实现 / 设备未连时抛 HardwareUnavailable。
    """
    if backend == "cpu":
        return CPUSampler(seed=seed)
    if backend == "mtj":
        from .mtj import MTJSampler
        return MTJSampler(seed=seed, **kw)
    raise ValueError(f"未知后端: {backend}（可选 cpu | mtj）")
