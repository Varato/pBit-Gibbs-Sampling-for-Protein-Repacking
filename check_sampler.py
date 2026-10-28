"""采样器验收协议：任何 CategoricalSampler 都必须服从 Categorical(p)。

方法 = 对拍 + 噪声地板（沿用旧仓 phase1 传统）：

    精确 p（可枚举）                        ← ground truth
    被测采样器抽 N 次 → 频率 f               ← 待验
    CPU 基准（Gumbel-max，数学严格正确）抽 N 次 × R 组 → TV 噪声地板
    判据：TV(f, p) ≤ max(R 组基准的 TV) → PASS
          （读法：没有抽得比"任何一组正确基准"更差。取 max 而非
          95 分位 = 对多枚测试骰做 Bonferroni 式校正）

测试分布自动挑选：本能量表很"陡"（多数残基冷启动条件分布 top1≈1.0，
统计检验力≈0），故按有效面数（1/Σp²）选 spread 最大的两枚真实条件
分布做主力对拍，另保留一枚强峰分布当 dyn_range 截断偏置的探针 + 一枚
合成 uniform。

Run:  python check_sampler.py --backend cpu              # 基准自检
      python check_sampler.py --backend mtj --n 100      # 硬件接入后
Out:  out/check_sampler_<backend>.json
"""
import argparse
import json
from pathlib import Path

import numpy as np

from engine import cond_probs, load_table
from sampler import CPUSampler, make_sampler

table = load_table()
cold = table.cold_start()


def cond(i, kT):
    return cond_probs(table.cond(cold, i), kT)


def eff_faces(p):
    return 1.0 / float(np.sum(p ** 2))


# 各残基在两档 kT 下的条件分布，按有效面数排序，取信息量最大的两枚 + 强峰探针
cands = []
for i in range(table.N):
    for kT in (0.5, 1.0):
        p = cond(i, kT)
        cands.append((eff_faces(p), i, kT))
cands.sort(key=lambda t: -t[0])
SPREAD1, SPREAD2 = cands[0], cands[1]
PEAK = min(cands, key=lambda t: t[0])

TESTS = [
    ("uniform(合成)", np.full(table.K[SPREAD1[1]],
                              1.0 / table.K[SPREAD1[1]])),
    (f"{table.res[SPREAD1[1]]} 冷启动 @kT={SPREAD1[2]}"
     f"（{SPREAD1[0]:.1f} 有效面）", cond(SPREAD1[1], SPREAD1[2])),
    (f"{table.res[SPREAD2[1]]} 冷启动 @kT={SPREAD2[2]}"
     f"（{SPREAD2[0]:.1f} 有效面）", cond(SPREAD2[1], SPREAD2[2])),
    (f"{table.res[PEAK[1]]} 冷启动 @kT=0.5（强峰·截断探针）",
     cond(PEAK[1], 0.5)),
]


def sample_freq(sampler, p, n):
    cnt = np.zeros(len(p))
    for _ in range(n):
        cnt[sampler.categorical_sample(p)] += 1
    return cnt / n


def tv(f, p):
    return float(0.5 * np.abs(np.asarray(f) - np.asarray(p)).sum())


ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--backend", choices=["cpu", "mtj"], default="cpu")
ap.add_argument("--n", type=int, default=None,
                help="每枚分布抽几次（默认 cpu=50000 / 硬件=100）")
ap.add_argument("--floor-repeats", type=int, default=30)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

N = args.n or (50000 if args.backend == "cpu" else 100)
sampler = make_sampler(args.backend, seed=args.seed)

print(f"采样器验收：backend={sampler.name}({args.backend}), 每枚抽 N={N}, "
      f"地板={args.floor_repeats} 组 CPU 基准对拍\n")
results, all_pass = [], True
for label, p in TESTS:
    floor_tvs = [tv(sample_freq(CPUSampler(seed=7777 + r), p, N), p)
                 for r in range(args.floor_repeats)]
    floor95 = float(np.percentile(floor_tvs, 95))
    floor = float(max(floor_tvs))               # 判据地板（见文件头说明）
    f = sample_freq(sampler, p, N)
    tvb = tv(f, p)
    ok = tvb <= floor
    all_pass &= ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:<34} K={len(p):>3}  "
          f"TV={tvb:.4f}  地板max={floor:.4f}(95分位 {floor95:.4f})  "
          f"众数={int(np.argmax(f))}(精确 {int(np.argmax(p))})")
    results.append({"test": label, "K": len(p),
                    "TV": round(tvb, 5), "floor_max": round(floor, 5),
                    "floor95": round(floor95, 5),
                    "freq": f.round(5).tolist(), "p_exact": p.round(5).tolist(),
                    "pass": bool(ok)})
sampler.close()

verdict = ("全部 PASS —— 这个采样器可以接入 Gibbs 引擎（GUI / engine.py）"
           if all_pass else "存在 FAIL —— 先修采样器再上引擎（见 README）")
print(f"\n{verdict}")
outdir = Path("out")
outdir.mkdir(exist_ok=True)
fn = outdir / f"check_sampler_{args.backend}.json"
with open(fn, "w") as fo:
    json.dump({"backend": args.backend, "N": N, "results": results,
               "all_pass": bool(all_pass), "sampler_stats": sampler.stats()},
              fo, indent=1, ensure_ascii=False)
print("wrote", fn)
raise SystemExit(0 if all_pass else 1)
