"""采集/回放分离的采样 worker：全速采样 → 帧缓冲 → GUI 按速度档回放。

数据面：worker 只写 queue（控制消息）和 FrameBuffer（帧），绝不碰
tk / matplotlib 对象。控制面：pause_ev / stop_ev 两个 Event。

FrameBuffer 满后对折 min-max（相邻两帧一桶，保留 E 最小与最大两帧，
顺序不变）——峰谷保真，不会削峰谎报收敛；x 轴用帧自带的 sweep_idx，
对折后间隔不均也能如实显示。硬件后端一帧 = 13 次 FPGA 往返（秒级），
远低于抽稀阈值，同一条代码路径。

Run:  被 gui.py 导入
"""
import queue
import threading
import time
import traceback

from engine import SweepInterrupted

TICK_SLEEP = 0.001          # CPU 后端批间让出 GIL
PAUSE_POLL = 0.05


class FrameBuffer:
    """定长帧缓冲（默认 2000 帧 ≈ 60 KB），满则 min-max 对折。"""

    def __init__(self, cap=2000):
        self.frames = []
        self.cap = cap

    def append(self, frame):
        self.frames.append(frame)
        if len(self.frames) > self.cap:
            self._halve()

    def _halve(self):
        out = []
        fs = self.frames
        for a in range(0, len(fs), 2):
            b = fs[a:a + 2]
            if len(b) == 1:
                out.append(b[0])
                continue
            lo = min(b, key=lambda f: f.E)
            hi = max(b, key=lambda f: f.E)
            out.extend([lo, hi] if lo.sweep_idx <= hi.sweep_idx else [hi, lo])
        self.frames = out

    def clear(self):
        self.frames.clear()

    def __len__(self):
        return len(self.frames)


class GibbsWorker(threading.Thread):
    """引擎驱动线程。生命周期：start() → 循环 → （stop_ev）→ 收尾消息。"""

    def __init__(self, engine, ui_queue):
        super().__init__(daemon=True)
        self.engine = engine
        self.ui_queue = ui_queue            # worker → UI 单向消息
        self.pause_ev = threading.Event()
        self.stop_ev = threading.Event()
        self.buffer = FrameBuffer()
        self.blocking = getattr(engine.sampler, "blocking", False)

    def run(self):
        blocking = self.blocking

        def on_residue(i, sweep, kT):
            # blocking 后端逐残基上报；同时是暂停/停止的响应点
            if self.stop_ev.is_set() or self.pause_ev.is_set():
                raise SweepInterrupted
            self.ui_queue.put(("residue", sweep, i, kT))

        try:
            while not self.stop_ev.is_set():
                if self.pause_ev.is_set():
                    time.sleep(PAUSE_POLL)
                    continue
                f = self.engine.run_sweep(on_residue if blocking else None)
                self.buffer.append(f)
                if self.engine.errors:
                    for e in self.engine.errors:
                        self.ui_queue.put(("error", e))
                    self.engine.errors.clear()
                if not blocking:
                    time.sleep(TICK_SLEEP)   # CPU 后端让出 GIL
        except Exception:
            self.ui_queue.put(("error", traceback.format_exc()))
        finally:
            self.ui_queue.put(("stopped",))

    def shutdown(self, timeout=125.0):
        """停止并回收（绝不强杀；超时靠 daemon 兜底）。"""
        self.stop_ev.set()
        self.pause_ev.clear()
        self.join(timeout=timeout)
        return not self.is_alive()
