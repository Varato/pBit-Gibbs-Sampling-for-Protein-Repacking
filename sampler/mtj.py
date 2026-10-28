"""MTJ + FPGA 硬件后端（桩）。

接入 = 实现 categorical_sample(p)，接口契约（输入/返回/错误语义、可选
类属性）见 base.py 的注释；怎么实现（通信方式、协议、器件细节）由
硬件侧自定。桩的默认行为：构造即抛 HardwareUnavailable——GUI 弹
"硬件未接入"提示框，可切回 CPU 后端，演示链路不中断。实现完成后用
python check_sampler.py --backend mtj --n 100 验收。
"""
import numpy.typing as npt

from .base import CategoricalSampler, HardwareUnavailable


class MTJSampler(CategoricalSampler):
    """MTJ p-bit 阵列采样器（桩）。"""

    name = "mtj"
    blocking: bool = True             # TODO① 按器件实际情况改
    dyn_range: float | None = None    # TODO① 按器件实际情况填；None = 不限

    def __init__(self, seed: int = 0, **kw):
        """TODO② 建立设备连接；连不上抛 HardwareUnavailable。

        以下 raise 即桩的占位行为，接入时整段替换。
        """
        raise HardwareUnavailable(
            "MTJ 硬件后端尚未实现：请在 sampler/mtj.py 完成 TODO②，"
            "接口契约见 sampler/base.py。")

    def categorical_sample(self, p: npt.ArrayLike) -> int:
        """TODO② 按 base.CategoricalSampler 的契约实现。"""
        raise NotImplementedError(
            "MTJ 硬件后端尚未实现（sampler/mtj.py TODO②）。")
