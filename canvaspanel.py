"""骰子阵列面板：13 行 × K_i 方格，实心 = 当前采中的 rotamer。

视觉约定：
  - 每行 one-hot 恰一个实心格（填残基色），空格白底 + 残基色描边
  - 本 sweep 翻转过的行加灰色背景带（#C9C9C9）
  - 硬件"采样中"的行复用背景带（淡黄）
  - native 最近候选格角上叠小三角（按 χ 偏差着色：黑 ≤30°、橙 ≤45°、
    红 >45° = MISSING）

绘制策略：启动/resize 时一次 create 全量 item，之后只 itemconfig 差异
更新（每 sweep ≤26 格 + ≤13 带切换），绝不删除重建。

Run:  被 gui.py 导入
"""
import tkinter as tk

import numpy as np

C_BAND = "#C9C9C9"          # 翻转行背景带
C_SIM = "#FFF3C2"           # 硬件采样中行的背景带
COV_COLORS = ((30.0, "#333333"), (45.0, "#E08214"), (1e9, "#B2182B"))
LABEL_W = 52                # 行标签列宽（像素）


class DiceCanvasPanel(tk.Canvas):
    def __init__(self, parent, table, palette, **kw):
        super().__init__(parent, background="white", highlightthickness=0, **kw)
        self.table = table
        self.palette = palette
        self.cell_ids = []          # cell_ids[i][k] -> rect item id
        self.band_ids = []          # 13 条背景带（z 序最底）
        self.native_ids = []        # 13 个 native 小三角（z 序最顶）
        self.label_ids = []
        self.last_state = None
        self._sim_rows = set()      # 当前"采样中"高亮的行
        self._pending_resize = None
        self.bind("<Configure>", self._on_configure)

    # ---- 构建 ---------------------------------------------------------------
    def build(self):
        """（重新）创建全部 item。resize debounce 后全量重建也走这里。"""
        self.delete("all")
        self.cell_ids, self.band_ids, self.native_ids, self.label_ids = \
            [], [], [], []
        self.last_state = None
        w, h = self.winfo_width(), self.winfo_height()
        if w < 50 or h < 50:
            return
        n = self.table.N
        kmax = int(self.table.K.max())
        row_h = h / n
        # 正方形格子：边长取"可用宽/最大K"与"行高×0.72"的较小者
        cw = max((w - LABEL_W - 8) / kmax, 2.0)
        s = max(min(cw, row_h * 0.72), 2.0)
        for i in range(n):
            res, Ki = self.table.res[i], int(self.table.K[i])
            c = self.palette[res]
            y0 = i * row_h
            # 背景带（最底层，默认隐藏）
            self.band_ids.append(self.create_rectangle(
                LABEL_W - 6, y0 + 1, LABEL_W + kmax * s + 4, y0 + row_h - 1,
                fill=C_BAND, outline="", state="hidden"))
            # 行标签
            self.label_ids.append(self.create_text(
                LABEL_W - 12, y0 + row_h / 2, text=res, fill=c,
                anchor="e", font=("Segoe UI", 9, "bold")))
            # 方格行（正方形，行内左对齐、行间垂直居中）
            row = []
            y = y0 + (row_h - s) / 2
            for k in range(Ki):
                x = LABEL_W + k * s
                row.append(self.create_rectangle(
                    x, y, x + max(s - 1, 1), y + s,
                    fill="white", outline=c, width=1))
            self.cell_ids.append(row)
            # native 小三角（叠在 native_k 格右下角）
            nk = self.table.native_k[i]
            dev = self.table.native_dev[i]
            mc = next(col for lim, col in COV_COLORS if dev <= lim)
            x = LABEL_W + nk * s + max(s - 1, 1)
            yb = y0 + (row_h + s) / 2
            ts = max(s * 0.45, 3)
            self.native_ids.append(self.create_polygon(
                x, yb, x - ts, yb, x, yb - ts, fill=mc, outline=""))

    # ---- 帧显示 --------------------------------------------------------------
    def show_frame(self, state, changed_rows=None):
        """diff 更新：state 为 (13,) 下标数组；changed_rows 为 None 时
        按与上次显示态的 diff 自动决定灰带行。"""
        state = np.asarray(state)
        if not self.cell_ids:             # 尚未 build（窗口未布局）：
            self.last_state = state.copy()    # 记下态，rebuild 后补画
            return
        if self.last_state is None:
            diff = np.ones(self.table.N, dtype=bool)
        else:
            diff = self.last_state != state
        for i in range(self.table.N):
            if diff[i]:
                c = self.palette[self.table.res[i]]
                if self.last_state is not None:
                    self.itemconfig(self.cell_ids[i][int(self.last_state[i])],
                                    fill="white", outline=c, width=1)
                # 实心格加黑描边：小尺寸下 one-hot 位置更醒目
                self.itemconfig(self.cell_ids[i][int(state[i])], fill=c,
                                outline="#222222", width=1)
        bands = diff if changed_rows is None else np.asarray(changed_rows,
                                                             dtype=bool)
        for i in range(self.table.N):
            sim = i in self._sim_rows
            self.itemconfig(self.band_ids[i],
                            state="normal" if (bands[i] or sim) else "hidden",
                            fill=C_SIM if sim else C_BAND)
        self.last_state = state.copy()

    def set_sim_row(self, i, on):
        """硬件逐残基进度：on=True 该行"采样中"（淡黄带）。"""
        if on:
            self._sim_rows.add(i)
        else:
            self._sim_rows.discard(i)
        if self.last_state is not None:
            sim = i in self._sim_rows
            shown = self.itemcget(self.band_ids[i], "state") == "normal"
            self.itemconfig(self.band_ids[i],
                            state="normal" if (sim or shown) else "hidden",
                            fill=C_SIM if sim else C_BAND)

    def clear_sim_rows(self):
        for i in list(self._sim_rows):
            self.set_sim_row(i, False)

    # ---- resize --------------------------------------------------------------
    def _on_configure(self, _event):
        if self._pending_resize is not None:
            self.after_cancel(self._pending_resize)
        self._pending_resize = self.after(100, self._rebuild)

    def _rebuild(self):
        self._pending_resize = None
        shown = None if self.last_state is None else self.last_state.copy()
        self.build()
        if shown is not None:
            self.show_frame(shown, changed_rows=np.zeros(self.table.N,
                                                         dtype=bool))
            for i in self._sim_rows:
                self.set_sim_row(i, True)
