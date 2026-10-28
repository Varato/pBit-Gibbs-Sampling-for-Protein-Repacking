"""采样器接口：整个项目对硬件方暴露的唯一契约。

科研想法：Gibbs 采样每步唯一的随机操作，是从一个 K 维条件分布里抽
一个类别（Categorical）。把这一步交给物理器件（MTJ p-bit 阵列 +
FPGA），物理熵源就替换了伪随机数发生器。因此对硬件的全部要求就是
一个函数：

    k = categorical_sample(p)     # p = 离散概率向量，k = 抽中的下标

温度 / 能量 / Boltzmann 折算全部在引擎侧（engine.py 的 cond_probs）
完成——硬件方拿到的 p 已经是归一概率，不需要知道 Gibbs 是什么。

合作方只需子类化 CategoricalSampler 实现 categorical_sample（见
sampler/mtj.py 的桩与 TODO）；cpu.py 是数学严格的参考实现，用于
对拍与演示。

Run:  被导入（engine / gui / check_sampler），不直接运行
"""
from abc import ABC, abstractmethod
from typing import TypedDict

import numpy as np
import numpy.typing as npt


class HardwareUnavailable(RuntimeError):
    """硬件后端不可用（未实现 / 设备未连接 / 通信失败）。

    GUI 捕获后弹提示框，不进入采样循环；引擎侧不做重试、立即终止。
    """


class SamplerStats(TypedDict, total=False):
    """stats() 的返回类型：运行统计（状态栏与结果 JSON 用）。

    total=False：键均可省略；引擎只读 n_samples / n_fail（缺省按 0）。
    """

    backend: str
    n_samples: int
    n_fail: int


class CategoricalSampler(ABC):
    """任意 Categorical 分布采样器的统一接口（CPU / 硬件可插拔）。

    契约（check_sampler.py 按此验收）：
      输入 p   一维非负 float 数组，Σp ≈ 1（±1e-6 容差；实现可自行
              防御性归一）。长度 K **每次调用可变**——本演示问题
              K ∈ 2..42，实现应支持任意 K ≥ 1；硬件上限不足时抛
              HardwareUnavailable
      返回 k   Python int ∈ 0..K-1，且统计上 P(k) = p[k]——这是
              验收协议（check_sampler.py）的核心判据
      可重入：一次 Gibbs 采样会调用 N_sweep × N_residue 次（百万级），
              每次调用相互独立；p 只读（不要原地修改）

    可选类属性（硬件实现按实际情况覆盖）：
      blocking    一次调用是否秒级以上（FPGA 往返=True）。引擎据此
                  决定是否上报逐残基进度、GUI 显示"硬件采样中"。
      dyn_range   概率动态范围上限 ratio = p_max/p_min。物理器件的
                  概率分辨率有限，而 softmax 概率可横跨几十个数量级
                  （强峰时 1e-40）。设了它，调用方（engine）会先把 p
                  截到 [p_max/ratio, p_max] 并重归一再传入；None =
                  不限（不截断）。截断引入的偏置由 check_sampler.py
                  的强峰探针量化。
    """

    blocking: bool = False
    dyn_range: float | None = None
    name: str = "abstract"

    @abstractmethod
    def categorical_sample(self, p: npt.ArrayLike) -> int:
        """抽一个 k ~ Categorical(p)。"""

    def close(self) -> None:
        """释放资源（设备句柄、socket 等）。必须可重复调用（幂等）。"""

    def stats(self) -> SamplerStats:
        """运行统计（SamplerStats）：引擎只读 n_samples / n_fail。"""
        return {}


def clip_dyn_range(p: npt.ArrayLike, ratio: float) -> np.ndarray:
    """把概率向量截到 [p_max/ratio, p_max] 再重归一（见 dyn_range 注释）。

    截断有偏（强行抬高了尾部概率）——偏多少由 check_sampler.py 的
    强峰探针告诉你。ratio=1 即不截断。
    """
    p = np.asarray(p, dtype=float)
    q = np.maximum(p, p.max() / max(ratio, 1.0))
    return q / q.sum()
