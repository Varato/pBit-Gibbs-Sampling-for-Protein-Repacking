"""MTJ + FPGA 硬件后端（桩）。★ 合作方看这里：全文搜 "TODO"（共 2 处）

对方已有成熟的 MTJ + FPGA 任意 Categorical 分布采样方案；本仓不关心
方案细节（SOT/HOT 磁化翻转、WTA 读出、FPGA 位流……），只要求实现
base.CategoricalSampler 的 categorical_sample(p)。

桩的默认行为：**构造即抛 HardwareUnavailable**——GUI 会弹"硬件未接入"
提示框，保证无硬件时演示链路不中断。接入步骤：

  ① TODO① dyn_range / blocking 类属性按器件实际填
  ② TODO② __init__ 建立与 FPGA 的连接（以太网 / 串口 / 厂商 SDK…），
     连不上抛 HardwareUnavailable；并在 categorical_sample 里完成
     "写入 p → 触发一次采样 → 读回获胜下标"

实现提示（都是与本演示问题的对接细节）：
  - K 每次调用不同（2..42）：设备侧需支持变长输入；或一次配置最大 K、
    尾部补 0——只要补 0 的面永不获胜即可
  - 器件概率分辨率有限（p_max/p_min 只有几十 dB）就设 dyn_range，
    引擎会预先截断归一，别在硬件里偷偷截
  - 一次往返秒级以上就设 blocking=True：引擎逐残基上报进度，GUI 显示
    "硬件采样中"高亮行
  - 单次失败（超时 / 校验错）抛 RuntimeError：引擎重试 ≤3 次，仍败则
    该残基保持原态继续（演示不中断，误差记入 errors）；
    HardwareUnavailable / NotImplementedError 则立即终止采样
  - 验收：python check_sampler.py --backend mtj --n 100 全 PASS 再上 GUI
"""
import numpy.typing as npt

from .base import CategoricalSampler, HardwareUnavailable


class MTJSampler(CategoricalSampler):
    """MTJ p-bit 阵列采样器（桩）。"""

    name = "mtj"
    blocking: bool = True             # TODO① 一次采样 = 一次 FPGA 往返，按实际改
    dyn_range: float | None = None    # TODO① 器件概率动态范围上限（如 1e3）；None=不限

    def __init__(self, seed: int = 0, **kw):
        """TODO② 连接参数（IP / 端口 / 串口 / 位流路径…）与建链。

        在这里建立与 FPGA 的连接并做一次自检；连不上抛 HardwareUnavailable。
        以下 raise 即为桩的占位行为，接入时整段替换。
        """
        raise HardwareUnavailable(
            "MTJ 硬件后端尚未实现：请在 sampler/mtj.py 完成 TODO②"
            "（连接设备 + 实现 categorical_sample），"
            "详见 README「硬件对接」一节。")

    def categorical_sample(self, p: npt.ArrayLike) -> int:
        """TODO② 写入 p → 触发采样 → 返回获胜下标 int ∈ 0..K-1。

        伪代码（按实际通信协议替换）：
            send(np.asarray(p))          # 或量化后的权重
            trigger()                    # 触发一次 MTJ 阵列弛豫
            return int(read_winner())    # 读回获胜位
        """
        raise NotImplementedError(
            "MTJ 硬件后端尚未实现（sampler/mtj.py TODO②）。")
