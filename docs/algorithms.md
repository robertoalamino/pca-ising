# Simulation and Analysis Algorithms

- This page describes the algorithms implemented in the companion code.
- Workflows that perform simulations generate fresh samples using fixed seeds.
- Figure numbers identify the analyses in the paper.
- The Figure 6 workflow reuses archived simulation profiles for its sampled points.
- All implementation details given below refer to the companion workflows unless explicitly attributed to the archived research notebooks.

## Square Lattice: Metropolis Sampling

The lattice has periodic boundaries and zero external field. For a spin
$s_i$ with four neighbours, the energy change for flipping that spin is

$$
\Delta E_i = 2 J s_i \sum_{j\sim i} s_j .
$$

The Metropolis acceptance probability is $\min(1,\exp[-\beta\Delta E_i])$;
see [Landau and Binder (2005)](#landau-and-binder-2005) for the
standard Metropolis method (reference [32] in the paper).
Sites are grouped so that no two neighbours belong to the same group.
This allows the flip probabilities to be computed for every spin in a
group at once. Even side lengths use two groups in a checkerboard
pattern. Odd side lengths use three groups because the checkerboard
pattern would put some neighbours across the periodic boundary in
the same group. Each sweep updates the groups in turn, using the
latest spin values to compute the next group's flip probabilities.

The temperature sweep uses **Metropolis updates with per-spin laziness**, meaning that
each proposed spin flip is independently skipped with probability $1/2$.
Otherwise, the usual Metropolis acceptance rule is applied. This halves
the probability of a flip, preserves the target distribution, and
prevents the deterministic alternation at $\beta=0$, when every spin
would otherwise flip on every sweep. After thermalisation, configurations
are recorded at a fixed interval measured in sweeps. The recorded
configurations can still be autocorrelated.

In the companion, this sampler is used for the workflows corresponding to Figures 1 and 2 and the square-lattice panel of Figure 8.

## Cayley Tree: Independent Equilibrium Sampling

For a finite tree with free boundaries and zero external field, draw the
root spin uniformly from $\{-1,+1\}$. Visit the remaining sites in
parent-before-child order. Each child copies its parent with probability

$$
p_{\mathrm{same}} = \frac{1+\tanh(\beta J)}{2}
$$

and takes the opposite value otherwise. The edge choices are independent.
This works because the zero-field Gibbs distribution factorises in the
root spin and the edge variables $s_{\mathrm{parent}}s_{\mathrm{child}}$.
A fresh traversal gives an independent equilibrium configuration, so
burn-in and thinning are unnecessary.

The parent-to-child sampling rule is the broadcast representation of
the free-boundary Ising model described by
[Evans et al. (2000)](#evans-et-al-2000). For related background on
the Ising model on trees and the role of boundary conditions, see
[Pemantle and Peres (2010)](#pemantle-and-peres-2010), reference [35]
in the paper.

In the companion, this sampler is used for the workflows corresponding
to Figures 3–5 and the Cayley-tree panel of Figure 8.

## PCA and the Overlap with Magnetisation

Each row of the sample matrix $X$ is a spin configuration; each column
is a site. The companion workflows use the following two covariance
estimates. The figure numbers refer to the analyses each workflow
corresponds to.

- **Ordinary sample centring** (Figures 3–5): subtract the observed mean
  of each column. Figure 3 extracts the leading right singular vector
  of the centred matrix; Figures 4–5 diagonalise its sample covariance.
- **Spin-flip symmetrisation** (Figures 1, 2 and 8): treat every
  configuration and its global flip as a pair. The resulting empirical
  mean is zero and the covariance estimate is $X^{\mathsf T}X/n$, where
  $n$ is the number of recorded configurations before pairing. This
  imposes the zero-field symmetry without physically duplicating
  the data.

For the standard covariance-eigenvector formulation of PCA, see
[Bishop (2006, section 12.1)](#bishop-2006), reference [23] in the paper.

For a unit leading eigenvector $v_1$ on $N$ sites, the sign-independent
overlap with the uniform magnetisation direction is

$$
\alpha = |v_1^{\mathsf T}\eta|,\qquad
\eta = \frac{(1,\ldots,1)}{\sqrt N}.
$$

The square-lattice size sweep uses power iteration with products
$X^{\mathsf T}(Xv)$, avoiding construction of the full covariance
matrix. After choosing the sign of $v_1$ so that $v_1^{\mathsf T}\eta
\geq 0$, it computes the residual as
$1-\alpha=\tfrac12\lVert v_1-\eta\rVert^2$. This avoids subtracting two nearly equal numbers. See the
[Figure 1 note](#figure-1-square-lattice-size-sweep) for how the
companion plot and the archived paper routine aggregate this residual.

## Exact Covariance on a Binary Tree

At zero field the population mean vanishes, and the exact covariance is

$$
C_{ij} = [\tanh(\beta J)]^{d(i,j)},
$$

where $d(i,j)$ is the number of edges on the unique path between sites
$i$ and $j$. The 255-site example diagonalises this matrix directly
to obtain a population reference for the sampled PC1.

The appendix to the companion paper derives this distance formula
for the two-point function. The reduction to generation indicators
uses the symmetry of the rooted binary tree, as described below.

For the size sweep and boundary calculations, the code works in the
orthonormal basis of generation indicators. At the positive coupling
and inverse temperature used here, the leading eigenvector is unique
and invariant under permutations of sites induced by rooted-tree
symmetries. Its entries therefore agree within each generation.

For a tree of depth $D$, the reduced covariance has only $D+1$ rows
and columns. If $u_\ell$ is the basis coefficient at generation
$\ell$, the PC1 value per site in that generation is
$v_\ell=u_\ell/\sqrt{2^\ell}$. The boundary analysis plots
$R_\ell=v_\ell/v_0$ for the first four non-root generations
as the buffer grows.

For Figure 6, each archived repetition is normalised by its **own**
root value before averaging the ratios across repetitions. The plotted
error bars are sample standard deviations across repetitions. The
exact curves are recalculated independently. Figure 7 uses only the
exact calculation.

## Temperature Sweep Observables

At each inverse temperature, the code computes the PC1 overlap, the
relative spectral gap $(\lambda_1-\lambda_2)/\lambda_1$ for the two
largest covariance eigenvalues, and the sample mean of the absolute
magnetisation per site. The companion plot for Figure 8 shows averages across repetitions.
The shaded bands show one sample standard deviation across runs.
The reference inverse temperatures used to rescale the horizontal
axes are specified in the plotting script; they are not fitted to
the finite-system results.

## Reproducing the Figures

### Figure 1: Square-lattice size sweep

The companion plot takes the natural logarithm of the **mean residual**
at each size, $\ln(\overline{1-\alpha})$. The archived plotting routine
for the paper instead takes $\ln(1-\alpha)$ for each repetition and
then computes the **mean and standard deviation of those log values**
at each size. The two operations generally give different results.
The archived sweep covers $L=3,\ldots,39$. The finite-sample overlap
stays close to one over the sizes tested; the plot does not establish
that the residual vanishes as the lattice grows.

### Figure 3: Cayley-tree size sweep

The companion reproduces the decline in PC1 overlap with tree size,
but the decline is smaller. At $N=2047$, the companion gives
$0.870 \pm 0.006$ (mean ± one standard deviation across 20 runs),
compared with about 0.6 read from the paper figure. The exact
population overlap at this size is 0.906.

The companion uses 500 independent equilibrium configurations per
run. In a separate check at $N=255$, the sampled overlap approached
the exact value as the number of configurations increased. A smaller
effective sample size in the original sweep could contribute to
lower overlaps and greater variation across runs. The original
sampling and PCA settings have not been verified, so the cause of
the difference at $N=2047$ remains unresolved.

### Figure 8: Temperature sweeps

The plotted magnetisation statistic is the mean absolute
magnetisation per site, $\langle |m| \rangle$. The gap is measured
relative to the largest covariance eigenvalue. The vertical lines
mark the reference critical inverse temperatures used to scale the
horizontal axes; they do not indicate singular transitions in the
finite graphs.

The archived research notebooks, which are separate from this
companion, contain candidate runs for both panels. A match to every
setting used for the final paper graphic has not been established.
To run the companion version, use the Figure 8 commands in the
[README's reproduction table](../README.md#reproduce-the-analyses).

## References

- <a id="landau-and-binder-2005" name="landau-and-binder-2005"></a>Landau, D. P., and Binder, K. (2005).
  *A Guide to Monte Carlo Simulations in Statistical Physics*.
  Cambridge University Press.
- <a id="evans-et-al-2000" name="evans-et-al-2000"></a>Evans, W., Kenyon, C., Peres, Y., and Schulman, L. J. (2000).
  [Broadcasting on trees and the Ising model](https://doi.org/10.1214/aoap/1019487349).
  *The Annals of Applied Probability*, 10(2), 410–433.
- <a id="pemantle-and-peres-2010" name="pemantle-and-peres-2010"></a>Pemantle, R., and Peres, Y. (2010).
  [The critical Ising model on trees, concave recursions and nonlinear capacity](https://doi.org/10.1214/09-AOP482).
  *The Annals of Probability*, 38(1), 184–206.
- <a id="bishop-2006" name="bishop-2006"></a>Bishop, C. M. (2006).
  *Pattern Recognition and Machine Learning*. Springer. Section 12.1.