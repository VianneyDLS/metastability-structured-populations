# Interactive browser simulator

Companion to *Metastability of a stochastic competitive model with large phenotypic diversity* (Vincent Bansaye, Vianney de la Salle and Gaël Raoul).

## Open locally

Open `docs/index.html` in a modern browser. The page is self-contained: no CDN, Python, package installation, backend or external requests. The calculations run in a Web Worker. If a browser blocks workers on `file://`, serve the folder locally instead:

```sh
python3 -m http.server 8000
```

Then open `http://localhost:8000/`.

## Publish on GitHub Pages

1. Add the `docs/index.html` file to the `main` branch of the existing `metastability-structured-populations` repository (preserve the `docs/` directory).
2. Open the repository **Settings → Pages**.
3. Under Build and deployment, choose **Deploy from a branch**, branch **main**, folder **/docs**; click **Save**.
4. The Pages settings will display the public URL once the deployment is complete. For the account `gaelraoul` and that repository name it would be `https://gaelraoul.github.io/metastability-structured-populations/`.

## Using the simulator

- Select A in (0,2), integer K, final time, and initial population mass M.
- Select the shape of f0: uniform, truncated Gaussian, two Gaussian peaks, linear tilt, or custom piecewise linear density (enter nonnegative heights at equally spaced points of [-1,1]). Shapes are normalized to mass M. The article's default is M=0.5, uniform f0=0.25.
- The stochastic initial population has floor(K M) individuals, each with a trait sampled independently from the normalized shape. The initial number is limited to 10,000, K to 20,000, and final time to 10,000. Custom input is limited to 2–1001 heights.
- Presets load the paper's parameters for Figures 1 and 2; they are not the original realizations. Click Run simulation after selecting a preset. Stop cancels an ongoing computation.
- The timeline displays the empirical atoms, reference ESD, deterministic density, and Gaussian-smoothed empirical distribution. Changing the smoothing width only changes the display.
- The energy plot supports a full linear range, logarithmic scale, article-style zoom (0.08 for dimorphic K=600, 0.006 for dimorphic K=10000, 0.005 for the monomorphic article preset), and custom upper limit. The time series has roughly unit-time spacing for long runs up to 5000, independently of the snapshot count; longer runs use a coarser time grid. A line connects discrete observations, not every individual jump.
- Exports: CSV time series, CSV selected snapshot (empirical and deterministic rows with different trait grids), all results and parameters in JSON, and energy plot in SVG. In JSON, `times` labels snapshots, `seriesTimes` labels energy/mass/diversity, and `stride` maps snapshot index to series index. `snap` stores integer counts, `x` the initial traits, `det` the deterministic density on `g`; `energy`, `Ed`, `mass`, `Md`, and `div` are the two energies, two masses, and surviving trait count respectively.

## Mathematical and numerical implementation

For x,y in [-1,1], a(x)=A-x², b(x)=1+a(x), d(x)=1, and c(x,y)=1/(1+(x-y)²). The direct Gillespie algorithm samples the exact jump process with density-dependent competition divided by K. It uses incremental competition updates and removes extinct traits. A seeded 32-bit PRNG makes results repeatable within this implementation; it does not match NumPy's random sequence. The deterministic problem is a midpoint-quadrature discretization in space (cell-averaged custom initial data) and a second-order exponential midpoint update with a configurable maximum time step and rate-based step cap. The deterministic initial density is normalized discretely. Numerical errors can be checked by increasing spatial resolution and reducing the time step; the Python solver remains the reference for publication-quality reproduction. The browser plots and snapshots are not intended to demonstrate the asymptotic theorem by themselves. Initial profiles lacking mass near the reference ESD need not converge to it. Long exact stochastic runs may be computationally intensive, especially on a phone.

## Verification

The page was exercised in headless Chromium with uniform, bimodal and custom initial profiles; the default K=600, A=1.5, T=100 experiment completed. Initial masses, decreasing deterministic energy on the tested trajectory, nonincreasing surviving traits, seed reproducibility, snapshot indexing, and CSV export were checked. Full article presets have not been benchmarked in the browser. No license is imposed by this addition; choose one with the authors of the repository.
