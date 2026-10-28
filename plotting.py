"""能量轨迹面板：matplotlib Figure 直接嵌入 tk（不用 pyplot）。

样式：
  灰 #cccccc 全轨迹 + 蓝 #2166AC 已播段 + 红 #B2182B 游标点 +
  蓝虚线 GMEC 参考线；退火模式右轴对数 kT（橙 #E08214），恒温模式隐藏。

刷新策略：图元常驻、只 set_data；draw_idle 时间戳节流 ≤12 fps
（全轴重绘一次 30–80 ms，12 fps 占主线程 <50%）。帧缓冲已在上游
min-max 对折到 ≤2000 点，这里不再抽稀。

Run:  被 gui.py 导入
"""
import time

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

C_BLUE, C_ORANGE, C_RED, C_GRAY = "#2166AC", "#E08214", "#B2182B", "#cccccc"
KB = 0.001987                # kcal/(mol·K)：kT(REU) → 开尔文
DRAW_MIN_INTERVAL = 1.0 / 12


class EnergyPlotPanel:
    def __init__(self, parent, gmec_E):
        # constrained_layout：twinx 副轴的刻度不会被推出画布（tight_layout
        # 与 twinx 不兼容，副轴刻度会被裁掉）
        self.fig = Figure(figsize=(7.2, 4.6), dpi=100,
                          constrained_layout=True)
        self.ax = self.fig.add_subplot(111)
        self.axT = self.ax.twinx()
        self.gmec_E = gmec_E

        self.line_gray, = self.ax.plot([], [], lw=0.9, color=C_GRAY, zorder=1)
        self.line_blue, = self.ax.plot([], [], lw=1.4, color=C_BLUE, zorder=2)
        self.ax.axhline(gmec_E, ls="--", lw=1.2, color=C_BLUE, alpha=0.6)
        self.ax.text(0.995, gmec_E, f" best {gmec_E:.2f} ", fontsize=9,
                     color=C_BLUE, ha="right", va="bottom",
                     transform=self.ax.get_yaxis_transform())
        self.cursor, = self.ax.plot([], [], "o", ms=7, color=C_RED, zorder=3)
        self.kt_line, = self.axT.plot([], [], lw=1.6, color=C_ORANGE,
                                      alpha=0.85)
        self.axT.set_yscale("log")
        self.axT.set_ylim(0.015, 5)
        self.axT.set_yticks([0.02, 0.1, 0.616, 1, 3])
        self.axT.set_yticklabels(["0.02", "0.1", "0.62\n(310K)", "1", "3"])
        self.axT.tick_params(labelsize=8, colors=C_ORANGE)
        self.ax.set_xlabel("sweep", fontsize=10)
        self.ax.set_ylabel("E (REU)", fontsize=10)
        self.ax.tick_params(labelsize=9)
        self.set_anneal(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        # 尺寸跟随交给 FigureCanvasTkAgg 内部的 Configure 绑定；
        # 布局交给 constrained_layout——外部再 set_size_inches 会和
        # 内部 resize 打架（画布被压回小尺寸、内容缩在角落）。
        self._last_draw = 0.0

    def set_anneal(self, on):
        """右轴 kT 只在退火模式显示。"""
        self.kt_line.set_visible(on)
        self.axT.set_visible(on)
        self.fig.canvas.draw_idle()

    def update(self, frames, playhead):
        """frames: FrameBuffer 的帧列表；playhead: 当前回放下标。"""
        if not frames:
            return
        ph = min(playhead, len(frames) - 1)
        sw = np.array([f.sweep_idx for f in frames])
        E = np.array([f.E for f in frames])
        kT = np.array([f.kT for f in frames])
        self.line_gray.set_data(sw, E)
        self.line_blue.set_data(sw[:ph + 1], E[:ph + 1])
        self.cursor.set_data([sw[ph]], [E[ph]])
        if self.kt_line.get_visible():
            self.kt_line.set_data(sw[:ph + 1], kT[:ph + 1])
        # ylim：已播段自适应 + 5% pad（GMEC 线保证可见）
        lo = min(E[:ph + 1].min(), self.gmec_E)
        hi = max(E[:ph + 1].max(), self.gmec_E)
        pad = max((hi - lo) * 0.05, 0.1)
        self.ax.set_ylim(lo - pad, hi + pad)
        self.ax.set_xlim(0, max(sw[-1], 10))
        now = time.monotonic()
        if now - self._last_draw >= DRAW_MIN_INTERVAL:
            self._last_draw = now
            self.canvas.draw_idle()
