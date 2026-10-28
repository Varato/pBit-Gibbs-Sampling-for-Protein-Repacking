"""主窗口：控制面 + 骰子阵列 + 能量轨迹 + 状态栏。

线程模型（见 worker.py）：worker 线程全速采样进帧缓冲；UI 线程
after(33ms) 轮询——消费控制消息、按速度档推进回放游标、diff 刷新
Canvas 与 matplotlib（draw_idle 节流 ≤12 fps）。worker 绝不碰 tk 对象。

后端固定两项（不再选网表/配置文件）：
  CPU（NumPy 参考实现）     Gumbel-max，数学精确，演示基准
  硬件（MTJ + FPGA）       sampler/mtj.py 未实现/设备未连时弹提示

控制语义：
  换 kT / 退火开关 = 软更新（下一 sweep 生效，轨迹不重置）
  换后端 / 重置 / seed = 重建（旧 worker join → sampler.close() → 新建）
  暂停 = 残基边界生效（硬件 blocking 下最多等一次采样）
Run:  python main.py
"""
import queue
import tkinter as tk
from tkinter import messagebox, ttk

from canvaspanel import DiceCanvasPanel
from engine import GibbsEngine, kT_const, kT_geometric, load_table
from plotting import EnergyPlotPanel
from sampler import HardwareUnavailable, make_sampler
from worker import GibbsWorker

KB = 0.001987                     # kT(REU) → 开尔文
POLL_MS = 33

# 13 色调色板（沿 phase1，与旧仓 Gibbs-Sampling-With-Memristors 同款）
PALETTE = {"N1": "#D95F02", "L2": "#7570B3", "Y3": "#1B9E77",
           "I4": "#E7298A", "Q5": "#66A61E", "W6": "#E41A1C",
           "L7": "#377EB8", "K8": "#FF7F00", "D9": "#17BECF",
           "S13": "#A6761D", "S14": "#E6AB02", "R16": "#984EA3",
           "S20": "#F781BF"}

CPU_LABEL = "CPU（NumPy 参考实现）"
MTJ_LABEL = "硬件（MTJ + FPGA）"
BACKENDS = [CPU_LABEL, MTJ_LABEL]
SPEEDS = [("跟随", None), ("64×", 64.0), ("32×", 32.0), ("16×", 16.0),
          ("8×", 8.0), ("4×", 4.0), ("2×", 2.0), ("1×", 1.0), ("0.5×", 0.5)]


class GuiApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MTJ Gibbs 采样实验室 —— 1L2Y 侧链 repacking")
        self.geometry("1280x760")
        self.table = load_table()
        self.ui_queue = queue.Queue()
        self.worker = None
        self.engine = None
        self.playhead = -1
        self._speed_acc = 0.0
        self._sim_row_cur = None
        self._build_widgets()
        self.canvas.show_frame(self.table.cold_start(),
                               changed_rows=[False] * self.table.N)
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self._poll_job = self.after(POLL_MS, self._poll)

    # ---- 布局 -----------------------------------------------------------------
    def _build_widgets(self):
        top = ttk.Frame(self, padding=(8, 6))
        top.pack(fill=tk.X)
        ttk.Label(top, text="采样后端").pack(side=tk.LEFT)
        self.cb_backend = ttk.Combobox(top, state="readonly", width=24,
                                       values=BACKENDS)
        self.cb_backend.current(0)
        self.cb_backend.pack(side=tk.LEFT, padx=(4, 10))
        self.cb_backend.bind("<<ComboboxSelected>>", self._on_backend_change)

        ttk.Label(top, text="kT").pack(side=tk.LEFT)
        self.sb_kt = ttk.Spinbox(top, from_=0.02, to=5.0, increment=0.1,
                                 width=5, command=self._live_schedule)
        self.sb_kt.set("1.0")
        self.sb_kt.pack(side=tk.LEFT, padx=(4, 2))
        self.sb_kt.bind("<KeyRelease>", lambda e: self._live_schedule())

        self.var_anneal = tk.BooleanVar(value=True)
        ck = ttk.Checkbutton(top, text="退火", variable=self.var_anneal,
                             command=self._live_schedule)
        ck.pack(side=tk.LEFT, padx=(8, 4))
        ttk.Label(top, text="退火步数").pack(side=tk.LEFT)
        self.sb_ann = ttk.Spinbox(top, from_=200, to=100000, increment=200,
                                  width=6, command=self._live_schedule)
        self.sb_ann.set("2000")
        self.sb_ann.pack(side=tk.LEFT, padx=(4, 10))

        ttk.Label(top, text="回放速度").pack(side=tk.LEFT)
        self.cb_speed = ttk.Combobox(top, state="readonly", width=6,
                                     values=[s for s, _ in SPEEDS])
        self.cb_speed.current(0)
        self.cb_speed.pack(side=tk.LEFT, padx=(4, 10))

        ttk.Label(top, text="seed").pack(side=tk.LEFT)
        self.ent_seed = ttk.Entry(top, width=5)
        self.ent_seed.insert(0, "0")
        self.ent_seed.pack(side=tk.LEFT, padx=(4, 12))

        self.btn_start = ttk.Button(top, text="开始采样", command=self._on_start)
        self.btn_start.pack(side=tk.LEFT, padx=2)
        self.btn_pause = ttk.Button(top, text="暂停", state=tk.DISABLED,
                                    command=self._on_pause)
        self.btn_pause.pack(side=tk.LEFT, padx=2)
        self.btn_reset = ttk.Button(top, text="重置", state=tk.DISABLED,
                                    command=self._on_reset)
        self.btn_reset.pack(side=tk.LEFT, padx=2)

        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.pack(fill=tk.BOTH, expand=True, padx=8, pady=6)
        left = ttk.Labelframe(main, text="骰子状态（13 残基 × 候选 rotamer，"
                                         "实心 = 当前态，角标三角 = native 最近）")
        main.add(left, weight=2)
        self.canvas = DiceCanvasPanel(left, self.table, PALETTE)
        self.canvas.pack(fill=tk.BOTH, expand=True)

        right = ttk.Labelframe(main, text="能量轨迹")
        main.add(right, weight=3)
        self.plot = EnergyPlotPanel(right, self.table.ref["gmec_E"])
        self.plot.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        bar = ttk.Frame(self, padding=(8, 2))
        bar.pack(fill=tk.X)
        self.lbl_stats = ttk.Label(bar, text="未开始", anchor=tk.W)
        self.lbl_stats.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.lbl_msg = ttk.Label(bar, text="", anchor=tk.E, foreground="#B2182B")
        self.lbl_msg.pack(side=tk.LEFT)

    # ---- 控制 -----------------------------------------------------------------
    def _seed(self):
        try:
            return int(self.ent_seed.get())
        except ValueError:
            return 0

    def _current_schedule(self):
        if self.var_anneal.get():
            n = max(int(float(self.sb_ann.get())), 50)
            return kT_geometric(t_hi=3.0, t_lo=0.02, n_anneal=n)
        return kT_const(max(float(self.sb_kt.get()), 1e-3))

    def _live_schedule(self):
        """kT / 退火开关 / 退火步数变更：软更新，不重置。"""
        if self.engine is not None:
            self.engine.set_kT_schedule(self._current_schedule())
            self.plot.set_anneal(self.var_anneal.get())

    def _make_sampler(self):
        """当前下拉选项 → 采样器实例；硬件不可用时弹提示并返回 None。"""
        sel = self.cb_backend.get()
        try:
            if sel == CPU_LABEL:
                return make_sampler("cpu", seed=self._seed())
            return make_sampler("mtj")          # seed：物理熵源用不上
        except HardwareUnavailable as err:
            messagebox.showinfo("硬件未接入", str(err))
            self._show_msg("硬件后端未实现/未连接——可改用 CPU 后端继续演示",
                           error=True)
            return None

    def _on_start(self):
        self._stop_worker()
        try:
            sampler = self._make_sampler()
        except Exception as err:                # 防御：工厂本身出错
            messagebox.showerror("采样后端出错", str(err))
            self._show_msg(f"后端出错：{err}", error=True)
            return
        if sampler is None:
            return
        self.engine = GibbsEngine(self.table, sampler, self._current_schedule(),
                                  seed=self._seed())
        self.worker = GibbsWorker(self.engine, self.ui_queue)
        self.plot.set_anneal(self.var_anneal.get())
        self.playhead, self._speed_acc = -1, 0.0
        self._sim_row_cur = None
        self.canvas.clear_sim_rows()
        self._set_running(True)
        self._show_msg("")
        self.worker.start()

    def _on_pause(self):
        if self.worker is None:
            return
        if self.worker.pause_ev.is_set():
            self.worker.pause_ev.clear()
            self.btn_pause.config(text="暂停")
        else:
            self.worker.pause_ev.set()
            self.btn_pause.config(text="继续" + ("…" if self.worker.blocking
                                                 else ""))
            if self.worker.blocking:
                self._show_msg("暂停将在当前残基采样完成后生效…")

    def _on_reset(self):
        self._stop_worker()
        self.playhead, self._speed_acc = -1, 0.0
        self._sim_row_cur = None
        self.canvas.clear_sim_rows()
        self.canvas.show_frame(self.table.cold_start(),
                               changed_rows=[False] * self.table.N)
        self.plot.update([], 0)
        self._set_running(False)
        self.lbl_stats.config(text="已重置（冷启动态）")
        self._show_msg("")

    def _on_backend_change(self, _event=None):
        if self.worker is not None:
            self._show_msg("后端已切换——下次「开始采样」生效（随机源语义已变，"
                           "不热切换）")

    def _set_running(self, running):
        self.btn_start.config(state=tk.NORMAL)
        self.btn_pause.config(state=tk.NORMAL if running else tk.DISABLED,
                              text="暂停")
        self.btn_reset.config(state=tk.NORMAL if running else tk.DISABLED)
        for w in (self.cb_backend, self.ent_seed):
            w.config(state=tk.DISABLED if running else (
                "readonly" if w is self.cb_backend else tk.NORMAL))

    def _stop_worker(self):
        if self.worker is not None:
            self.worker.shutdown()
            self.worker = None
        if self.engine is not None:
            self.engine.close()
            self.engine = None
        self.canvas.clear_sim_rows()
        self._set_running(False)

    def _on_close(self):
        if self._poll_job is not None:
            self.after_cancel(self._poll_job)
            self._poll_job = None
        self._stop_worker()
        self.destroy()

    # ---- 轮询 -----------------------------------------------------------------
    def _poll(self):
        if not self.winfo_exists():
            return
        try:
            while True:
                msg = self.ui_queue.get_nowait()
                if msg[0] == "residue":
                    _, sweep, i, kT = msg
                    if self._sim_row_cur is not None and self._sim_row_cur != i:
                        self.canvas.set_sim_row(self._sim_row_cur, False)
                    self._sim_row_cur = i
                    self.canvas.set_sim_row(i, True)
                    self._show_msg(f"sweep {sweep}：残基 {i + 1}/13 "
                                   f"({self.table.res[i]}) 硬件采样中… kT={kT:.3g}")
                elif msg[0] == "error":
                    self._show_msg(str(msg[1]).splitlines()[-1][:120],
                                   error=True)
                elif msg[0] == "stopped":
                    self._after_stop()
        except queue.Empty:
            pass

        if self.worker is not None and len(self.worker.buffer):
            buf = self.worker.buffer
            fpp = dict(SPEEDS)[self.cb_speed.get()]
            if fpp is None:                     # 跟随
                ph = len(buf) - 1
            else:
                self._speed_acc += fpp
                step = int(self._speed_acc)
                self._speed_acc -= step
                ph = min(self.playhead + step, len(buf) - 1)
            if ph != self.playhead:
                self.playhead = ph
                f = buf.frames[ph]
                self.canvas.show_frame(f.state)
                self.plot.update(buf.frames, ph)
                self.lbl_stats.config(
                    text=(f"sweep {f.sweep_idx}   E = {f.E:.3f} REU   "
                          f"best = {f.best_E:.3f}   kT = {f.kT:.3g} "
                          f"({f.kT / KB:.0f} K)   采样 {f.n_samples} 次"
                          + (f"（失败 {f.n_fail}）" if f.n_fail else "")))
        self._poll_job = self.after(POLL_MS, self._poll)

    def _after_stop(self):
        self._sim_row_cur = None
        self.canvas.clear_sim_rows()
        self._set_running(False)
        self._show_msg("采样已停止（回放仍可推进缓冲帧）")

    def _show_msg(self, text, error=False):
        self.lbl_msg.config(text=text, foreground="#B2182B" if error else "#333333")
