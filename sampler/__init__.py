"""采样器包：接口契约（base）+ CPU 参考（cpu）+ MTJ 硬件桩（mtj）。

常用入口：

    from sampler import CategoricalSampler, CPUSampler, make_sampler

    k = categorical_sample(p)      # 便捷函数：默认 CPU 参考实现

Run:  被导入，不直接运行
"""
from .base import (CategoricalSampler, HardwareUnavailable,  # noqa: F401
                   clip_dyn_range)
from .cpu import CPUSampler  # noqa: F401


def make_sampler(backend, seed=0, **kw):
    """名字 → 采样器实例 的统一工厂（GUI / check_sampler / selftest 共用）。

    backend ∈ {cpu, mtj}；mtj 未实现 / 设备未连时抛 HardwareUnavailable。
    """
    if backend == "cpu":
        return CPUSampler(seed=seed)
    if backend == "mtj":
        from .mtj import MTJSampler
        return MTJSampler(seed=seed, **kw)
    raise SystemExit(f"未知后端: {backend}（可选 cpu | mtj）")


_default = CPUSampler()


def categorical_sample(p):
    """模块级便捷函数（默认 CPU 参考实现）。

    等价于 CPUSampler().categorical_sample(p)，供快速试验与合作方对照
    接口签名用；正式链路请用 make_sampler(...) 注入可插拔实例。
    """
    return _default.categorical_sample(p)
