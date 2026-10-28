"""GUI 冒烟测试：构造主窗口 → 跑一轮 CPU 采样 → 每行恰一个实心格 → 关闭。

无显示环境（无 tk）时自动跳过。

Run:  python -m pytest tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

tk = pytest.importorskip("tkinter")
matplotlib = pytest.importorskip("matplotlib")   # noqa: F841

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gui import CPU_LABEL, GuiApp   # noqa: E402


def test_gui_smoke():
    app = GuiApp()
    app.update_idletasks()
    app.update()
    # 点「开始采样」（CPU 后端）跑一小段
    app.cb_backend.set(CPU_LABEL)
    app._on_start()
    for _ in range(20):
        app.update()
        if app.worker is not None and len(app.worker.buffer) >= 3:
            break
    app._stop_worker()
    # 每行恰一个实心格（one-hot 显示不变式）
    for row in app.canvas.cell_ids:
        fills = [app.canvas.itemcget(cid, "fill") for cid in row]
        solid = [f for f in fills if f not in ("white", "")]
        assert len(solid) == 1
    app._on_close()
    try:                       # 销毁后的 Tk 解释器：winfo 抛 TclError 即已关
        assert not app.winfo_exists()
    except tk.TclError:
        pass
