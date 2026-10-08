import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, json, os

OUT = os.path.dirname(os.path.abspath(__file__)) + "/fig"
os.makedirs(OUT, exist_ok=True)

# ---------- parametros del diseno ----------
G = 1 + 39e3/10e3            # 4.9
VSAT = 10.5
R7, R9, C1 = 100.0, 100e3, 1e-6
TAU = R9*C1                  # 0.1 s
TAUC = R7*C1                 # 0.1 ms
R10, R11 = 10e3, 430e3
VREF = 12*1e3/(15e3+1e3)     # 0.75
k = R10/R11
VTH = VREF*(1+k) + VSAT*k
VTL = VREF*(1+k) - VSAT*k
FS = 20000.0
DT = 1/FS

BLUE, ORANGE, AQUA, YEL, MAG, VIO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e6e5e0"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.8, "figure.facecolor": "white", "axes.facecolor": "white",
                     "legend.frameon": False, "lines.linewidth": 1.8})

def chain(vin):
    vg = np.clip(G*vin, -VSAT, VSAT)
    vhw = np.where(vg > 0, -vg, 0.0)
    vr = np.clip(-(vg + 2*vhw), -VSAT, VSAT)
    env = np.zeros_like(vr); vo = np.zeros_like(vr)
    e = 0.0; s = -VSAT
    dec = np.exp(-DT/TAU); chg = 1-np.exp(-DT/TAUC)
    for i, x in enumerate(vr):
        e *= dec
        if x > e: e += (x-e)*chg
        env[i] = e
        if s < 0 and e > VTH: s = VSAT
        elif s > 0 and e < VTL: s = -VSAT
        vo[i] = s
    return vg, vhw, vr, env, vo

def pulse(t, f, A, w=0.075, pol=1):
    return pol*A*(((t % (1/f)) < w).astype(float))

def stats(t, vo, nbeats_expected, f):
    hi = vo > 0
    rises = np.where(np.diff(hi.astype(int)) == 1)[0]
    falls = np.where(np.diff(hi.astype(int)) == -1)[0]
    widths = []
    for r in rises:
        fl = falls[falls > r]
        if len(fl): widths.append((fl[0]-r)*DT)
    stuck = hi[-1] and len(falls) < len(rises)-1
    return len(rises), (np.median(widths)*1e3 if widths else float("nan")), stuck

# ---------- tabla de casos con pulso ----------
rows = []
T_SIM = 10.0
t = np.arange(0, T_SIM, DT)
for bpm in [40, 72, 120, 150, 180]:
    f = bpm/60
    for A in [0.5, 1.0, 2.0]:
        for pol in [1, -1]:
            vin = pulse(t+0.2, f, A, pol=pol)   # +0.2 s: primer pulso no en t=0
            vg, vhw, vr, env, vo = chain(vin)
            n_exp = int(np.sum(np.diff((((t+0.2) % (1/f)) < 0.075).astype(int)) == 1))
            n, wms, stuck = stats(t, vo, n_exp, f)
            # minimo de la envolvente justo antes de cada latido (ultimo latido completo)
            starts = np.where(np.diff((((t+0.2) % (1/f)) < 0.075).astype(int)) == 1)[0]
            pre = [env[s0-1] for s0 in starts[1:]]
            rows.append(dict(bpm=bpm, Vpp=A, pol=pol, esperados=n_exp, detectados=n,
                             ancho_ms=round(wms, 0), pico=round(env.max(), 2),
                             env_min_antes=round(min(pre), 2) if pre else None))

# ---------- tiempo de descarga teorico vs simulado ----------
theo = []
for A in [0.5, 1.0, 2.0]:
    Vp = min(G*A, VSAT)
    td = TAU*np.log(Vp/VTL)
    theo.append(dict(Vpp=A, Vp=round(Vp, 2), t_desc_ms=round(td*1e3, 0),
                     ancho_teo_ms=round((0.075+td)*1e3, 0),
                     fmax_lpm=round(60/(0.075+td), 0)))

# ---------- sinusoide ----------
sine_rows = []
for fs_ in [0.5, 1.0, 1.5]:
    for A in [0.5, 1.0, 2.0]:
        ts = np.arange(0, 8, DT)
        vin = (A/2)*np.sin(2*np.pi*fs_*ts)
        *_, env, vo = chain(vin)
        n, wms, stuck = stats(ts, vo, 0, fs_)
        sine_rows.append(dict(f=fs_, Vpp=A, ciclos=int(8*fs_), pulsos=n,
                              pulsos_por_ciclo=round(n/(8*fs_), 2),
                              env_min=round(env[int(2/DT):].min(), 2)))

# ---------- tolerancias del rectificador (Monte Carlo + peor caso) ----------
def rect_gains(R2, R3, R4, R5, R6):
    gpos = (R6/R5)*(R3/R2) - R6/R4
    gneg = R6/R4
    return gpos, gneg
worst = {}
for tol in [0.01, 0.05]:
    best = 0
    for s2 in [-1, 1]:
        for s3 in [-1, 1]:
            for s4 in [-1, 1]:
                for s5 in [-1, 1]:
                    for s6 in [-1, 1]:
                        gp, gn = rect_gains(10e3*(1+s2*tol), 10e3*(1+s3*tol), 20e3*(1+s4*tol),
                                            10e3*(1+s5*tol), 20e3*(1+s6*tol))
                        best = max(best, abs(gp-gn)/gn)
    rng = np.random.default_rng(1)
    sd = tol/3
    R = [10e3, 10e3, 20e3, 10e3, 20e3]
    m = np.array([R]*20000)*(1+rng.normal(0, sd, (20000, 5)))
    gp, gn = rect_gains(*m.T)
    worst[tol] = dict(peor=round(best*100, 1), p95=round(np.percentile(np.abs(gp-gn)/gn, 95)*100, 2))

res = dict(G=G, TAU=TAU, VTH=VTH, VTL=VTL, VREF=VREF, k=k, theo=theo, rows=rows,
           sine=sine_rows, tol=worst)
json.dump(res, open(OUT + "/results.json", "w"), indent=1)

def thr(ax):
    ax.axhline(VTH, color=INK2, ls="--", lw=1)
    ax.axhline(VTL, color=INK2, ls=":", lw=1.2)
    ax.text(ax.get_xlim()[1], VTH, f" VTH {VTH:.2f} V".replace(".",","), va="center", fontsize=8.5, color=INK2)
    ax.text(ax.get_xlim()[1], VTL, f" VTL {VTL:.2f} V".replace(".",","), va="center", fontsize=8.5, color=INK2)

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

# F1 transferencia del rectificador
vgx = np.linspace(-10, 10, 401)
_, vhw, vr, *_ = chain(vgx/G) if False else (None, np.where(vgx > 0, -vgx, 0), np.abs(vgx))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(vgx, np.where(vgx > 0, -vgx, 0), color=ORANGE, label="Vhw (TP2)")
ax[0].plot(vgx, np.abs(vgx), color=BLUE, label="Vr (TP3)")
ax[0].set_xlabel("Vg (V)"); ax[0].set_ylabel("Salida (V)")
ax[0].set_title("Curva de transferencia", fontsize=10.5, color=INK, loc="left")
ax[0].legend(loc="lower right")
ts = np.arange(0, 2, DT); vin = 0.5*np.sin(2*np.pi*1*ts)
vg, vhw2, vr2, *_ = chain(vin)
ax[1].plot(ts, vg, color=MUTED, lw=1.2, label="Vg (TP1)")
ax[1].plot(ts, vhw2, color=ORANGE, label="Vhw (TP2)")
ax[1].plot(ts, vr2, color=BLUE, label="Vr (TP3)")
ax[1].set_xlabel("Tiempo (s)"); ax[1].set_ylabel("Tensión (V)")
ax[1].set_title("Sinusoide 1 Vpp, 1 Hz", fontsize=10.5, color=INK, loc="left")
ax[1].set_ylim(-3.7, 2.9); ax[1].legend(loc="lower center", ncol=3, fontsize=8.5)
fig.tight_layout(); save(fig, "f1_rectificador")

# F2 formas de onda pulso 1 Vpp 72 lpm
ts = np.arange(0, 3.4, DT); f = 1.2
vin = pulse(ts+0.6, f, 1.0)
vg, vhw, vr, env, vo = chain(vin)
fig, ax = plt.subplots(5, 1, figsize=(9, 8.6), sharex=True)
for a, y, c, lab in zip(ax, [vin, vg, vr, env, vo], [INK2, BLUE, BLUE, ORANGE, AQUA],
                        ["TP0 entrada (V)", "TP1 Vg (V)", "TP3 Vr (V)", "TP4 Venv (V)", "TP5 Vo (V)"]):
    a.plot(ts, y, color=c); a.set_ylabel(lab)
ax[2].plot(ts, vhw, color=ORANGE, lw=1.2, label="TP2 Vhw"); ax[2].legend(loc="lower right", fontsize=8.5)
ax[3].set_xlim(0, 3.4); thr(ax[3])
ax[-1].set_xlabel("Tiempo (s)")
ax[0].set_title("Pulso de 1 Vpp a 72 lpm (1,2 Hz), ancho 75 ms", fontsize=11, color=INK, loc="left")
fig.tight_layout(); save(fig, "f2_pulso_72lpm")

# F3 descarga de la envolvente
fig, ax = plt.subplots(figsize=(9, 3.6))
tt = np.linspace(0, 0.6, 600)
for A, c in zip([0.5, 1.0, 2.0], [AQUA, BLUE, VIO]):
    Vp = G*A; v = Vp*np.exp(-tt/TAU); td = TAU*np.log(Vp/VTL)
    ax.plot(tt*1e3, v, color=c, label=f"{A} Vpp (pico {Vp:.2f} V)")
    ax.plot([td*1e3], [VTL], "o", color=c, ms=6, mec="white", mew=1.5)
    ax.annotate(f"{td*1e3:.0f} ms", (td*1e3, VTL), textcoords="offset points", xytext=(4, 8),
                fontsize=8.5, color=INK2)
ax.set_xlim(0, 600); ax.set_ylim(0, 10.5); thr(ax)
ax.set_xlabel("Tiempo desde el final del pulso (ms)"); ax.set_ylabel("Venv (V)")
ax.set_title("Descarga de C1 con τ = 100 ms hasta VTL", fontsize=11, color=INK, loc="left")
ax.legend(loc="upper right"); fig.tight_layout(); save(fig, "f3_descarga")

# F4 frecuencia maxima vs amplitud y vs tau
fig, ax = plt.subplots(figsize=(9, 3.6))
A = np.linspace(0.25, 2.1, 200)
for tau_, c, lab in [(0.068, AQUA, "τ = 68 ms"), (0.1, BLUE, "τ = 100 ms (diseño)"), (0.15, ORANGE, "τ = 150 ms")]:
    Vp = np.minimum(G*A, VSAT)
    fmax = 60/(0.075 + tau_*np.log(np.maximum(Vp/VTL, 1.0001)))
    ax.plot(A, fmax, color=c, label=lab)
ax.axhline(150, color=INK2, ls="--", lw=1); ax.text(2.1, 150, " 150 lpm", va="center", fontsize=8.5, color=INK2)
ax.axvline(VTH/G, color=MUTED, ls=":", lw=1.2)
ax.text(VTH/G, 30, " umbral de\n detección", fontsize=8.5, color=INK2)
ax.set_ylim(0, 420); ax.set_xlim(0.2, 2.1)
ax.set_xlabel("Amplitud de entrada (Vpp)"); ax.set_ylabel("Frecuencia máxima (lpm)")
ax.set_title("Frecuencia cardiaca máxima que se libera a tiempo: 60 / (w + τ·ln(Vp/VTL))", fontsize=10.5, color=INK, loc="left")
ax.legend(loc="upper right"); fig.tight_layout(); save(fig, "f4_fmax")

# F5 lazo de histeresis
fig, ax = plt.subplots(figsize=(6, 3.6))
x = np.concatenate([np.linspace(0, 2, 400), np.linspace(2, 0, 400)])
s = -VSAT; y = []
for xi in x:
    if s < 0 and xi > VTH: s = VSAT
    elif s > 0 and xi < VTL: s = -VSAT
    y.append(s)
ax.plot(x, y, color=BLUE)
for xv, lab, ha, dx in [(VTL, f"VTL = {VTL:.2f} V".replace(".",","), "right", -0.08), (VTH, f"VTH = {VTH:.2f} V".replace(".",","), "left", 0.08)]:
    ax.axvline(xv, color=MUTED, ls=":", lw=1.2); ax.text(xv+dx, 0.6, lab, fontsize=8.5, color=INK2, ha=ha)
ax.annotate("", xy=(VTH+0.06, 6), xytext=(VTH+0.06, -6), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4))
ax.annotate("", xy=(VTL-0.06, -6), xytext=(VTL-0.06, 6), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4)); ax.text(1.25, VSAT-2.2, "Venv sube: dispara en VTH", fontsize=8.5, color=INK2); ax.text(0.03, -VSAT+1.4, "Venv baja: libera en VTL", fontsize=8.5, color=INK2)
ax.set_xlabel("Venv (TP4) (V)"); ax.set_ylabel("Vo (TP5) (V)")
ax.set_title("Lazo de histéresis del Schmitt (U6)", fontsize=11, color=INK, loc="left")
fig.tight_layout(); save(fig, "f5_histeresis")

# F6 sinusoide 0.5 Hz y 1.5 Hz, 2 Vpp
fig, ax = plt.subplots(3, 2, figsize=(9, 6), sharex="col")
for j, fs_ in enumerate([0.5, 1.5]):
    ts = np.arange(0, 4, DT); vin = 1.0*np.sin(2*np.pi*fs_*ts)
    vg, vhw, vr, env, vo = chain(vin)
    ax[0, j].plot(ts, vg, color=BLUE); ax[0, j].set_title(f"Sinusoide 2 Vpp, {fs_} Hz".replace(".", ","), fontsize=10.5, color=INK, loc="left")
    ax[1, j].plot(ts, vr, color=MUTED, lw=1.2, label="Vr"); ax[1, j].plot(ts, env, color=ORANGE, label="Venv")
    ax[1, j].set_xlim(0, 4); thr(ax[1, j]) if j == 1 else (ax[1, j].axhline(VTH, color=INK2, ls="--", lw=1), ax[1, j].axhline(VTL, color=INK2, ls=":", lw=1.2))
    ax[2, j].plot(ts, vo, color=AQUA); ax[2, j].set_xlabel("Tiempo (s)")
ax[0, 0].set_ylabel("TP1 Vg (V)"); ax[1, 0].set_ylabel("TP3/TP4 (V)"); ax[2, 0].set_ylabel("TP5 Vo (V)")
ax[1, 0].legend(loc="upper right", fontsize=8.5)
fig.tight_layout(); save(fig, "f6_sinusoide")

# F7 latido sintetico PWLIN, ambas polaridades
pts = np.array([[0, 0], [0.300, 0], [0.320, -0.10], [0.345, 1.00], [0.370, -0.25], [0.390, 0],
                [0.500, 0], [0.580, 0.20], [0.660, 0], [0.833, 0]])
def pwl(t):
    return np.interp(t % 0.8333, pts[:, 0], pts[:, 1])
fig, ax = plt.subplots(3, 2, figsize=(9, 6), sharex=True)
for j, pol in enumerate([1, -1]):
    ts = np.arange(0, 2.5, DT); vin = pol*pwl(ts)
    vg, vhw, vr, env, vo = chain(vin)
    ax[0, j].plot(ts, vin, color=INK2)
    ax[0, j].set_title("QRS positivo" if pol > 0 else "QRS invertido", fontsize=10.5, color=INK, loc="left")
    ax[1, j].plot(ts, vr, color=MUTED, lw=1.2, label="Vr"); ax[1, j].plot(ts, env, color=ORANGE, label="Venv")
    ax[1, j].axhline(VTH, color=INK2, ls="--", lw=1); ax[1, j].axhline(VTL, color=INK2, ls=":", lw=1.2)
    ax[2, j].plot(ts, vo, color=AQUA); ax[2, j].set_xlabel("Tiempo (s)")
ax[0, 0].set_ylabel("TP0 (V)"); ax[1, 0].set_ylabel("TP3/TP4 (V)"); ax[2, 0].set_ylabel("TP5 Vo (V)")
ax[1, 0].legend(loc="upper right", fontsize=8.5)
fig.tight_layout(); save(fig, "f7_qrs_pwlin")

print(json.dumps(dict(G=G, VTH=VTH, VTL=VTL, theo=theo, tol=worst), indent=1))
for r in rows: print(r)
for r in sine_rows: print(r)
