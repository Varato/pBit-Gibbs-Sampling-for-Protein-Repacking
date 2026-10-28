"""MTJ Gibbs 采样实验室：入口。

  python main.py

数据自包含（data/*.json 入库），不需要 PyRosetta / LTSpice。
"""
from gui import GuiApp

if __name__ == "__main__":
    GuiApp().mainloop()
