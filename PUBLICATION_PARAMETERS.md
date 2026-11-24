# Parameters of the numerical experiments

Reference: Vincent Bansaye, Vianney de la Salle and Gaël Raoul, *Metastability of a stochastic competitive model with large phenotypic diversity*, manuscript dated 24 November 2026, Section 4 (pp. 27–31), Figures 1–2. The final bibliographic reference can be added when available.

## 1. Model specified in the manuscript

The trait space is $\mathcal X=[-1,1]$. For $x,y\in\mathcal X$,

$$
a(x)=A-x^2,\qquad b(x)=1+A-x^2,\qquad d(x)=1,\qquad c(x,y)=\frac{1}{1+(x-y)^2}.
$$

The limiting initial measure is $\xi_0(dx)=f_0(x)\,dx$, with $f_0(x)=1/4$ on $[-1,1]$. Its total mass is $1/2$. In each stochastic experiment, $N_0^K=\lfloor K/2\rfloor$, and initial traits are independent and uniform on $[-1,1]$. Offspring inherit parental traits without mutation.

The empirical population is $\nu_t^K=K^{-1}\sum_{i=1}^{N_t^K}\delta_{X_t^{K,i}}$. For a trait $x_i$ with $n_i$ individuals, the total birth rate is $n_i b(x_i)$ and the total death rate is

$$
n_i\left(1+\frac{1}{K}\sum_j c(x_i,x_j)n_j\right).
$$

The deterministic comparison solves

$$
\partial_t f(t,x)=f(t,x)\left(A-x^2-\int_{-1}^1\frac{f(t,y)}{1+(x-y)^2}\,dy\right).
$$

The ESD is $\bar\xi=A\delta_0$ for $0<A\leq1$. For $1<A<2$, it is $m(\delta_{-z}+\delta_z)$, where

$$
z=\frac{\sqrt{-7+\sqrt{17+32A}}}{4},\qquad m=\frac{13+16A-3\sqrt{17+32A}}{16}.
$$

The energy (equation (18) in the manuscript) is

$$
\mathcal E(\mu)=\frac12\iint c(x,y)\,\mu(dx)\mu(dy)-\int a(x)\,\mu(dx)+\frac12\iint c(x,y)\,\bar\xi(dx)\bar\xi(dy).
$$

## 2. Publication experiments

| Figure | A | K | Initial individuals | Final time | Snapshot times | Energy display range |
|---|---:|---:|---:|---:|---|---|
| 1(a–e) | 0.5 | 10000 | 5000 | 5000 | 0, 50, 887, 3000 | [0, 0.005] |
| 2(a,b) | 1.5 | 600 | 300 | 1000 | 1000 | [0, 0.08] |
| 2(c,d) | 1.5 | 10000 | 5000 | 3000 | 3000 | [0, 0.006] |

Figure 1(b) uses a Gaussian convolution of the empirical measure with standard deviation 0.05 (variance 0.0025), for visualization only. Figure 1 uses the neighborhood half-width h=0.2 for the support-concentration illustration. Figure 2 uses the energy tick intervals 0.02 and 0.002 respectively. The manuscript's Figure 2(d) time axis extends to 5000, whereas its caption and Section 4.2 specify a final time of 3000; the reproduction uses 3000. Energy values above the zoomed limits are not visible in those panels; the initial energy is annotated and the CSV includes the full series.

## 3. Numerical choices of the supplied reproduction code

The following are choices made for this repository, not claimed to be the original publication settings.

| Setting | Python reproduction |
|---|---|
| Stochastic solver | Exact direct Gillespie simulation, no time discretization |
| Random seeds (Fig. 1; Fig. 2 K=600; Fig. 2 K=10000) | 20261124; 20261125; 20261126 |
| Deterministic spatial grid | 801 midpoint nodes on [-1,1] |
| Deterministic time solver | SciPy DOP853, adaptive |
| Relative / absolute tolerances | 1e-7 / 1e-10 |
| Energy and mass recording interval | 1 time unit, including snapshots and final time |
| Smoothing standard deviation | 0.05 |

The original trajectories, original random seeds and original deterministic solver settings were not supplied. The examples here are new stochastic realizations using the manuscript's model parameters. The browser uses a different random generator and deterministic solver (see [WEB_README.md](WEB_README.md)); the browser presets select the publication's model parameters, not all of the Python numerical settings.

## 4. Reproduction commands

After installing the dependencies, run from the repository root:

```sh
python simulate.py --case mono --seed 20261124
python simulate.py --case dimorphic --seed 20261124
# Or both in one run:
python simulate.py --case all --seed 20261124
```

The dimorphic runs automatically add 1 and 2 to the base seed. `--quick` runs smaller demonstrations and does not reproduce the parameter table above. See [README.md](README.md) for setup and outputs. The [examples](examples/README.md) directory contains the resulting CSV files, parameter records and figures from the new reference runs.
