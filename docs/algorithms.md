# Simulation and analysis algorithms

This page describes the algorithms implemented in the companion code.
The stochastic workflows use fixed seeds and generate fresh samples.
Figure 6 reuses archived simulation profiles for its sampled points.
The methods below describe the released workflows; they do not establish
the original settings where those remain unverified.

## Square lattice: Metropolis sampling

The lattice has periodic boundaries and zero external field. For a spin
$s_i$ with four neighbours, the energy change for flipping that spin is

$$
\Delta E_i = 2 J s_i \sum_{j\sim i} s_j .
$$

The Metropolis acceptance probability is
$\min(1,\exp[-\beta\Delta E_i])$. Sites are partitioned into sets of
nonadjacent sites, allowing updates within each set to be vectorised.
Even side lengths use two checkerboard colours. Odd side lengths use
three colours to account for the periodic wraparound. Each sweep
visits all colour sets in sequence, using the updated spins when
computing the next set's acceptance probabilities.

The temperature sweep uses **lazy Metropolis**: it halves each flip
probability, which is equivalent to independently skipping half the
proposals on average. This preserves the target distribution and avoids
the deterministic alternation at $\beta=0$, when every proposed flip
would otherwise be accepted. After thermalisation, configurations are
recorded at a fixed sweep interval. They can still be autocorrelated.

Used for Figures 1, 2 and the square-lattice panel of Figure 8.

## Cayley tree: independent equilibrium sampling

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

Used for Figures 3–5 and the Cayley-tree panel of Figure 8.

## PCA and the overlap with magnetisation

Each row of the sample matrix $X$ is a spin configuration; each column
is a site. The scripts use two covariance estimates:

- **Ordinary sample centring** (Figures 3–5): subtract the observed mean
  of each column. Figure 3 extracts the leading right singular vector
  of the centred matrix; Figures 4–5 diagonalise its sample covariance.
- **Spin-flip symmetrisation** (Figures 1, 2 and 8): treat every
  configuration and its global flip as a pair. The resulting empirical
  mean is zero and the covariance estimate is $X^{\mathsf T}X/n$, where
  $n$ is the number of recorded configurations before pairing. This
  imposes the zero-field symmetry without physically duplicating
  the data.

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
$1-\alpha=\tfrac12\lVert v_1-\eta\rVert^2$. This avoids subtracting
two nearly equal numbers. The plot takes the natural log of the
**mean residual** at each size.

## Exact covariance on a binary tree

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

## Temperature sweep observables

At each inverse temperature, the code computes the PC1 overlap, the
relative spectral gap $(\lambda_1-\lambda_2)/\lambda_1$ for the two
largest covariance eigenvalues, and the sample mean of the absolute
magnetisation per site. Figure 8 plots averages across repetitions.
The shaded bands show one sample standard deviation across runs.
The reference inverse temperatures used to rescale the horizontal
axes are specified in the plotting script; they are not fitted to
the finite-system results.

## References

- Landau, D. P., and Binder, K. (2005). *A Guide to Monte Carlo
  Simulations in Statistical Physics*. Cambridge University Press.
- Evans, W., Kenyon, C., Peres, Y., and Schulman, L. J. (2000).
  [Broadcasting on trees and the Ising model](https://doi.org/10.1214/aoap/1019487349).
  *The Annals of Applied Probability*, 10(2), 410–433.
- Pemantle, R., and Peres, Y. (2010).
  [The critical Ising model on trees, concave recursions and nonlinear capacity](https://doi.org/10.1214/09-AOP482).
  *The Annals of Probability*, 38(1), 184–206.
- Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*.
  Springer. Section 12.1.