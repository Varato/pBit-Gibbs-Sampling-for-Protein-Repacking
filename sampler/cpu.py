"""CPU 参考实现：NumPy Gumbel-max（与 Categorical 采样数学严格等价）。

argmax(log p + Gumbel) 的分布严格等于 Categorical(p)（Gumbel-max
定理），一次向量化、无拒绝；p 中的 0 天然免疫（log 0 = -inf 永不
选中）。这是对拍硬件的 ground-truth 基准，也是 GUI 的默认后端。

Run:  被导入（engine / gui / check_sampler），不直接运行
"""
import numpy as np

from .base import CategoricalSampler


class CPUSampler(CategoricalSampler):
    """软件基准骰：Gumbel-max。"""

    name = "cpu"

    def __init__(self, seed=0):
        self.rng = np.random.default_rng(seed)
        self.n_samples = self.n_fail = 0

    def categorical_sample(self, p):
        q = np.asarray(p, dtype=float)
        q = q / q.sum()                      # 防御性归一（契约允许 ±1e-6）
        with np.errstate(divide="ignore"):
            z = np.log(q)                    # p=0 → -inf，永不选中
        k = int(np.argmax(z + self.rng.gumbel(size=z.size)))
        self.n_samples += 1
        return k

    def stats(self):
        return {"backend": self.name, "n_samples": self.n_samples,
                "n_fail": self.n_fail}
