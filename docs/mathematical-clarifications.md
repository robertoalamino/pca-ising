# Mathematical Clarifications

These notes accompany *When does PCA recover the magnetisation? The principal component of the ferromagnetic Ising model*.

The details presented below are intended to make the assumptions and arguments more explicit and better justified, but were not ready at the time of publication and are included here for the reader's benefit. Equation numbers are those found in the published version.

## Finite Trees: Scope of the Non-alignment Result

For the homogeneous zero-field ferromagnetic Ising model on a finite
tree, the non-alignment result applies to trees with at least three
vertices and at fixed $0<\beta J<\infty$. A one-vertex tree is trivial,
and the two-vertex tree has equal covariance row sums and a uniform
leading eigenvector.

Let $C$ be the exact spin covariance matrix and define the row sum as in the article by
$C_i=\sum_j C_{ij}$. If $C\eta=\lambda\eta$ with
$\eta=N^{-1/2}(1,\ldots,1)^{\mathsf T}$, each row sum equals $\lambda$.
For a leaf $\ell$ and its unique neighbour $p$, with $t=\tanh(\beta J)$, the row sums satisfy
$C_\ell-C_p=(1-t)(1+t-C_p)$. Because $p$ has at least two neighbours
in a tree with at least three vertices, $C_p\ge 1+2t$. The row sums
therefore differ, so $\eta$ cannot be an eigenvector.

## Non-alignment along a Sequence of Rooted Trees

Let $T_n$ be a sequence of finite rooted trees with
$N=|V(T_n)|\to\infty$. Write $\partial T_n$ for the set of non-root
leaves and $P_n$ for the set of their parents. On each tree, consider
the exact covariance of the homogeneous zero-field ferromagnetic
Ising model with free boundaries. Let $v$ be its normalised positive
leading eigenvector and let $\eta$ be the uniform unit vector on
that tree.

At fixed $t=\tanh(\beta J)\in(0,1)$, the leaf–parent ratio bound in
equation (48) gives $v_p\ge c\,v_\ell$, with $c>1$ independent of
the tree size. Choose one leaf child $\ell(p)$ per parent in $P_n$.
The selected pairs are disjoint. With $a=N^{-1/2}$, the squared
distance from $(a,a)$ to the half-plane $x\ge cy$ gives

$$
(v_p-a)^2+(v_{\ell(p)}-a)^2
\ge \frac{a^2(c-1)^2}{1+c^2}.
$$

Summing over the selected pairs gives

$$
\|v-\eta\|^2
\ge \frac{|P_n|}{N}\frac{(c-1)^2}{1+c^2}.
$$

Thus $\liminf_{n\to\infty}|P_n|/N>0$ rules out alignment in the limit.
A positive lower limiting leaf density and a uniform bound $k$ on the
number of leaf children per parent suffice, since
$|\partial T_n|\le k|P_n|$. This is a bound on the full vector norm that
does not require pointwise convergence of the rescaled components.

## Translation Invariance and Combining Data at Different Temperatures

On a periodic lattice with translation-invariant couplings, for a
Gibbs ensemble that shares the translation symmetry, the covariance
kernel depends on the displacement vector $\mathbf r_j-\mathbf r_i$.
Translation invariance alone does not imply that it depends only on
the scalar distance. With Fourier modes $e^{-i\mathbf k\cdot\mathbf r_i}$
and displacement $\mathbf r=\mathbf r_j-\mathbf r_i$, the Fourier
sum uses $e^{-i\mathbf k\cdot\mathbf r}$. The kernel is even, so changing
the sign convention leaves the eigenvalues unchanged.

The between-temperature outer product in equation (19) is positive
semidefinite, but can have negative off-diagonal entries. In the
finite zero-field Gibbs population, spin-flip symmetry makes the mean
configuration zero at every temperature. The between-temperature
term then vanishes, and the weighted population covariances preserve
the uniform leading direction on the periodic ferromagnetic lattice.
Empirical covariance matrices need not preserve that symmetry. The
quantity $\eta^{\mathsf T}C\eta$ is always the variance along the
magnetisation direction. It is the top eigenvalue when alignment holds.

## Temperature Limits and the Infinite-tree Comparison

For a finite connected ferromagnetic graph in the symmetric
zero-field Gibbs limit as $\beta\to\infty$, the covariance matrix is
the all-ones matrix. Its leading direction is uniform, with eigenvalue
$N$ and with the other eigenvalues being zero. If the ensemble is restricted
to a single fully ordered configuration, the centred covariance becomes zero.

Vertex transitivity makes sites equivalent when the Gibbs state
respects that symmetry. On the exactly infinite Bethe lattice,
however, the constant vector is not square-summable and does not
define a normalised PC1 in $\ell^2$. The non-alignment result above
concerns a sequence of finite rooted trees and their finite-dimensional
leading eigenvectors.

## Appendix - Edge Expansion

In equation (A.5), the spin sum factors over all $N$ sites of the
original graph. A site incident to no selected edge has degree zero
and contributes a factor of two. Every nonempty selected edge set of a tree with at least one edge contains a
vertex of odd degree. For the two-point function, the insertion
$s_i s_j$ changes the parity at $i,j$, leaving only the unique path
between them. The resulting partition function and two-point function are unchanged.