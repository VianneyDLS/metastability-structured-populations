"""Numerical experiments for Metastability of a stochastic competitive model
with large phenotypic diversity. Run python simulate.py --help."""
import argparse
import json
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def growth(x, A):
    return A - np.asarray(x)**2


def competition(x, y):
    return 1 / (1 + (np.asarray(x)-np.asarray(y))**2)


def esd(A):
    """Locations and masses of the evolutionarily stable distribution."""
    if not 0 < A < 2:
        raise ValueError('Require 0 < A < 2')
    if A <= 1:
        return np.array([0.]), np.array([A])
    x = np.sqrt(-7+np.sqrt(17+32*A))/4
    m = (13+16*A-3*np.sqrt(17+32*A))/16
    return np.array([-x,x]), np.array([m,m])


def equilibrium_constant(A):
    x,m = esd(A)
    return float(m @ competition(x[:,None],x[None,:]) @ m / 2)


def energy(A, x, m):
    """Energy of a discrete measure, with masses m (not densities). Eq. (18)."""
    x,m = np.asarray(x),np.asarray(m)
    return equilibrium_constant(A) - m @ growth(x,A) + m @ competition(x[:,None],x[None,:]) @ m / 2


@njit(cache=True)
def _gillespie(x0, K, A, times, snapshot_indices, const, seed):
    """Exact direct Gillespie algorithm; one initial individual per trait.
    Only living traits are kept in the active part of the arrays. Competition
    is updated after every jump, in O(number of surviving traits) operations.
    """
    np.random.seed(seed)
    x = x0.copy()
    active = len(x)
    ids = np.arange(active)
    n = np.ones(active, dtype=np.int64)
    b = 1+A-x*x
    comp = np.zeros(active)
    for i in range(active):
        for j in range(active):
            comp[i] += 1/(K*(1+(x[i]-x[j])**2))
    cum = np.zeros(2*active)
    E = np.empty(len(times))
    total = np.empty(len(times))
    diversity = np.empty(len(times), dtype=np.int64)
    counts = np.zeros((len(snapshot_indices),active),dtype=np.int64)
    t, r, s, events = 0., 0, 0, 0
    while r < len(times):
        rate = 0.
        for j in range(active):
            rate += n[j]*b[j]
            cum[2*j] = rate
            rate += n[j]*(1+max(0.,comp[j]))
            cum[2*j+1] = rate
        next_t = t + np.random.exponential(1/rate) if rate > 0 else np.inf
        # Record the state BEFORE the next jump at all intervening times.
        while r < len(times) and (times[r] < next_t or rate == 0):
            e, mass = const, 0.
            for j in range(active):
                m = n[j]/K
                e += m*(comp[j]/2-(b[j]-1))
                mass += m
            E[r],total[r],diversity[r] = e,mass,active
            if s < len(snapshot_indices) and r == snapshot_indices[s]:
                for j in range(active):
                    counts[s,ids[j]] = n[j]
                s += 1
            r += 1
        if r == len(times):
            break
        k = np.searchsorted(cum[:2*active], np.random.random()*rate, side='right')
        k = min(k,2*active-1)
        j = k//2
        delta = 1 if k%2 == 0 else -1
        n[j] += delta
        for i in range(active):
            comp[i] += delta/(K*(1+(x[i]-x[j])**2))
        if n[j] == 0:
            active -= 1
            x[j],n[j],b[j],ids[j],comp[j] = x[active],n[active],b[active],ids[active],comp[active]
        t = next_t
        events += 1
        # Limit floating point drift from incremental updates.
        if events%100000 == 0:
            for i in range(active):
                comp[i] = 0.
                for j in range(active):
                    comp[i] += n[j]/(K*(1+(x[i]-x[j])**2))
    return counts,E,total,diversity,events


def stochastic(A,K,times,snapshots,seed):
    """Initial population floor(K/2), i.i.d. uniform on [-1,1]."""
    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(-1,1,K//2))
    idx = np.searchsorted(times,snapshots)
    if not np.array_equal(times[idx],snapshots):
        raise ValueError('Include snapshots in recording times.')
    event_seed = int(rng.integers(0,2**32,dtype=np.uint32))
    return (x,)+_gillespie(x,K,A,times,idx,equilibrium_constant(A),event_seed)


def deterministic(A,times,grid_size=801):
    """Midpoint quadrature in x and adaptive DOP853 integration in t."""
    dx = 2/grid_size
    x = -1 + dx*(np.arange(grid_size)+.5)
    C = competition(x[:,None],x[None,:])
    a = growth(x,A)
    def rhs(t,f):
        return f*(a-dx*(C@f))
    with threadpool_limits(limits=1,user_api='blas'):
        sol = solve_ivp(rhs,(0.,float(times[-1])),np.full(grid_size,.25),
                        t_eval=times,method='DOP853',rtol=1e-7,atol=1e-10)
        if not sol.success:
            raise RuntimeError(sol.message)
        f = sol.y.T
        m = dx*f
        E = equilibrium_constant(A)-m@a+.5*np.einsum('ij,ij->i',m@C,m)
    return x,f,E


def convolution(x,m,grid,eta=.05):
    """Gaussian convolution for visualization only, not mutation."""
    keep = m>0
    return np.exp(-.5*((grid[:,None]-x[None,keep])/eta)**2)@m[keep]/(np.sqrt(2*np.pi)*eta)


def snapshot_plot(ax,x,n,K,grid,f,A,t,smooth=False,initial=False,show_esd=True):
    ax.set(xlim=(-1,1),xlabel='Phenotypic trait',ylabel='Number of individuals / $K$',title=f'$t={t:g}$')
    handles = []
    if not initial:
        keep = n>0
        ax.vlines(x[keep],0,n[keep]/K,color='blue',lw=.9)
        ax.plot(x[keep],n[keep]/K,'.',color='blue',ms=3)
        handles.append(Line2D([0],[0],color='blue',marker='.',label=r'$\nu_t^K$'))
    if show_esd:
        z,m = esd(A)
        ax.vlines(z,0,m,color='black',linestyle='--',lw=1.5)
        ax.plot(z,m,'ko',ms=4)
        handles.append(Line2D([0],[0],color='black',ls='--',marker='o',label=r'$\bar\xi$'))
    ax.set_ylim(bottom=0)
    if handles:
        ax.legend(handles=handles,loc='upper left',fontsize=8)
    if f is not None:
        twin = ax.twinx()
        twin.plot(grid,f,color='red',lw=1.5,label=r'$f_t$')
        if smooth:
            twin.plot(grid,convolution(x,n/K,grid),color='blue',alpha=.7,label=r'$\Gamma*\nu_t^K$')
        twin.set_ylabel('Population density')
        twin.set_ylim(bottom=0)
        twin.legend(loc='upper right',fontsize=8)
    if initial:
        ax.set_yticks([])
        ax.set_ylabel('')


def energy_plot(ax,t,E,Ed=None,markers=(),zoom=False,energy_limits=None):
    ax.plot(t,E,color='blue',lw=1,label=r'$\mathcal{E}(\nu_t^K)$')
    if Ed is not None:
        ax.plot(t,Ed,color='red',lw=1.2,label=r'$\mathcal{E}(\xi_t)$')
    for s in markers:
        y = E[np.searchsorted(t,s)]
        ax.vlines(s,0,y,color='black',ls='--',lw=.8)
        ax.plot(s,y,'ko',ms=3)
    ax.set(xlabel='Time $t$',ylabel='Energy',xlim=(0,t[-1]))
    ax.set_ylim(bottom=0)
    if zoom:
        ax.set_ylim(0,.005)
        ax.annotate(f'Initial energy = {E[0]:.3g} (outside scale)',xy=(.02,.92),xycoords='axes fraction',fontsize=8)
    if energy_limits is not None:
        ymax, tick = energy_limits
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(round(ymax/tick)+1)*tick)
        if E[0] > ymax:
            ax.annotate(r'$\mathcal{E}(\nu_0^K)\approx$' + f' {E[0]:.2f}',
                        xy=(.035,.97), xytext=(.12,.80),
                        xycoords='axes fraction', textcoords='axes fraction',
                        color='blue', fontsize=9,
                        arrowprops=dict(arrowstyle='->', color='blue', linestyle='--'))
    ax.legend(fontsize=8, loc='upper right' if energy_limits is not None else 'best')


def dimorphic_energy_limits(K):
    """Energy-axis range and tick interval of manuscript Figure 2."""
    return {600: (.08, .02), 10000: (.006, .002)}.get(K)


def save(fig,path):
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'),dpi=160,bbox_inches='tight')
    plt.close(fig)


def csv(path,columns,arrays):
    np.savetxt(path,np.column_stack(arrays),delimiter=',',header=','.join(columns),comments='',fmt='%.16g')


def run_case(name,A,K,T,snapshots,args,seed):
    s = np.asarray(snapshots,dtype=float)
    t = np.unique(np.r_[np.arange(0.,T,args.record_step),T,s])
    print(f'{name}: A={A}, K={K}, T={T}, seed={seed}',flush=True)
    x,n,E,mass,div,events = stochastic(A,K,t,s,seed)
    print(f'  {events:,} exact jumps',flush=True)
    g,f,Ed = deterministic(A,t,args.grid_size)
    idx = np.searchsorted(t,s)
    out = args.output/name
    out.mkdir(parents=True,exist_ok=True)
    labels = [f't{v:g}' for v in s]
    csv(out/'empirical.csv',['traits']+labels,[x]+list(n/K))
    csv(out/'deterministic.csv',['x']+labels,[g]+list(f[idx]))
    csv(out/'energy.csv',['t','energy_sto','energy_det','total_mass','n_traits'],[t,E,Ed,mass,div])
    np.savez_compressed(out/'results.npz',times=t,snapshot_times=s,traits=x,counts=n,grid=g,density=f,energy_sto=E,energy_det=Ed,total_mass=mass,diversity=div,A=A,K=K)
    (out/'parameters.json').write_text(json.dumps(dict(A=A,K=K,final_time=T,snapshots=s.tolist(),seed=seed,grid_size=args.grid_size,record_step=args.record_step,events=events,eta=.05,rtol=1e-7,atol=1e-10),indent=2)+'\n')
    for i,v in enumerate(s):
        fig,ax = plt.subplots(figsize=(5,4),layout='constrained')
        snapshot_plot(ax,x,n[i],K,g,f[idx[i]] if name=='mono' else None,A,v,smooth=name=='mono' and i==1,initial=v==0,show_esd=i==len(s)-1)
        save(fig,out/f'snapshot_t{v:g}')
    fig,ax = plt.subplots(figsize=(9,3),layout='constrained')
    energy_plot(ax,t,E,Ed,markers=s[1:] if name=='mono' else (),zoom=name=='mono' and not args.quick)
    save(fig,out/'energy')
    return x,n,E,t,g,f,Ed,idx,s,K,A


def figure_one(r,args):
    x,n,E,t,g,f,Ed,idx,s,K,A = r
    fig = plt.figure(figsize=(10,10),layout='constrained')
    gs = fig.add_gridspec(3,2,height_ratios=[1,1,.65])
    for i,v in enumerate(s):
        ax = fig.add_subplot(gs[i//2,i%2])
        snapshot_plot(ax,x,n[i],K,g,f[idx[i]],A,v,smooth=i==1,initial=i==0,show_esd=i==3)
        ax.set_title(f'({chr(97+i)})  $t={v:g}$')
    ax = fig.add_subplot(gs[2,:]); energy_plot(ax,t,E,Ed,s[1:],not args.quick); ax.set_title('(e)')
    save(fig,args.output/'figure_1')


def figure_two(results,args):
    fig,axs = plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for i,r in enumerate(results):
        x,n,E,t,g,f,Ed,idx,s,K,A = r
        snapshot_plot(axs[i,0],x,n[-1],K,g,None,A,s[-1])
        axs[i,0].set_title(f'({chr(97+2*i)})  $K={K}$, $t={s[-1]:g}$')
        energy_plot(axs[i,1],t,E,energy_limits=None if args.quick else dimorphic_energy_limits(K))
        axs[i,1].set_title(f'({chr(98+2*i)})  $K={K}$')
    save(fig,args.output/'figure_2')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=['mono','dimorphic','all'],default='mono')
    p.add_argument('--quick',action='store_true',help='Small demonstration, NOT article parameters')
    p.add_argument('--seed',type=int,default=20261124)
    p.add_argument('--grid-size',type=int,default=801)
    p.add_argument('--record-step',type=float,default=1.)
    p.add_argument('--output',type=Path,default=Path('output'))
    args = p.parse_args()
    if args.grid_size<3 or args.record_step<=0 or not np.isfinite(args.record_step) or not 0<=args.seed<2**32-2:
        p.error('Invalid grid size, recording step or seed.')
    args.output.mkdir(parents=True,exist_ok=True)
    if args.case in ('mono','all'):
        r = run_case('mono',.5,200 if args.quick else 10000,30 if args.quick else 5000,[0,5,15,25] if args.quick else [0,50,887,3000],args,args.seed)
        figure_one(r,args)
    if args.case in ('dimorphic','all'):
        results = []
        for i,(K,T) in enumerate([(100,20),(200,30)] if args.quick else [(600,1000),(10000,3000)]):
            results.append(run_case(f'dimorphic_K{K}',1.5,K,T,[T],args,args.seed+i+1))
        figure_two(results,args)

if __name__ == '__main__':
    main()
