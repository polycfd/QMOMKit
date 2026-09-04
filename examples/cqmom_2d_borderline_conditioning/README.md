# 2D CQMOM Borderline Conditioning

## Purpose

This deterministic, static example studies a two-dimensional CQMOM
reconstruction as two x-nodes approach one another. It explains three related
effects:

- conditioning order changes where a weak inversion appears;
- comparing four solvers shows whether a failure depends on the linear solver;
- Wheeler reduction can construct a lower-node quadrature when the requested
  fixed two-node representation is no longer supported.

The central question is:

> Does the observed failure depend on the linear solver? If Dense, BP, pivoted
> QR, and SVD all fail the same two-node reconstruction, the two-node system
> has lost numerical independence; changing the solver cannot recover the
> missing information.

There is no physical model, time integration, regularization, or hidden tolerance relaxation. 
The sweep and both output files are deterministic.

## Prescribed Population

For separation parameter $\delta$, the four exact atoms are

$$
\left(-\frac{\delta}{2},-1\right),\quad w=0.10,
\qquad
\left(-\frac{\delta}{2},+1\right),\quad w=0.30,
$$

$$
\left(+\frac{\delta}{2},-1\right),\quad w=0.15,
\qquad
\left(+\frac{\delta}{2},+1\right),\quad w=0.45.
$$

The marginals are

$$
x=\left(-\frac{\delta}{2},+\frac{\delta}{2}\right),
\qquad w_x=(0.4,0.6),
$$

$$
y=(-1,+1),
\qquad w_y=(0.25,0.75).
$$

All raw moments with $0\le p,q\le3$ are computed directly from the atoms:

$$
M_{p,q}=\sum_j w_j x_j^p y_j^q.
$$

The tensor uses phase-space order $(p,q)$ and extents $(4,4)$. Checks include

$$
M_{0,0}=1,\qquad M_{1,0}=0.1\delta,\qquad M_{0,1}=0.5,
\qquad M_{2,0}=0.25\delta^2.
$$

## Where The Numerical Difficulty Comes From

In the x-first order, the weighted parent system is

$$
A_x=
\begin{bmatrix}
0.4 & 0.6\\
-0.2\delta & 0.3\delta
\end{bmatrix}.
$$

Its infinity-norm condition number behaves asymptotically as

$$
\kappa_\infty(A_x)\sim\frac{2.5}{\delta}.
$$

The example also computes the local 2-norm condition number from singular
values. The two condition numbers answer different questions: the infinity
norm is the explanatory asymptotic estimate, while $\kappa_2$ is the
scale-independent metric shown in the report.

For the reversed order, the parent matrix is

$$
A_y=
\begin{bmatrix}
0.25 & 0.75\\
-0.25 & 0.75
\end{bmatrix}.
$$

It remains well conditioned. Reversing the conditioning order relocates the
weak x-inversion to later Wheeler or conditional stages; it does not remove
the underlying degeneracy.

## Why Compare Several Solvers?

The fixed-count sweep requests four explicit policies:

- Dense;
- Björck–Pereyra (BP);
- pivoted QR; and
- SVD.

The four direct policies—Dense, BP, QR, and SVD—are compared to determine
whether a failure depends on the linear solver. Agreement shows that the same
formulation is being resolved consistently. Disagreement can reveal
solver-specific sensitivity, such as pivot detection or backend rank handling.
A common failure across all four means that the two-node system has lost
numerical independence, so selecting a different solver cannot restore the
missing information.

## What The Diagnostics Mean

### Relative singular value

The local matrix diagnostic reports

$$
r_\sigma=\frac{\sigma_{\min}}{\sigma_{\max}}.
$$

This is a scale-independent measure of numerical independence. A smaller
$r_\sigma$ means the columns are closer to dependence. The local rank indicator
uses the relative threshold $1\times10^{-14}$.

### Condition number

The local 2-norm condition number is

$$
\kappa_2=\frac{\sigma_{\max}}{\sigma_{\min}}=\frac{1}{r_\sigma}.
$$

It measures sensitivity to perturbations. A large condition number does not
automatically mean that the current solve must fail, and a finite condition
number is not a robustness guarantee.

### Moment error

For a successful quadrature, the report shows the worst scaled raw-moment
error

$$
E_M=\max_{p,q}
\frac{\left|\widehat M_{p,q}-M_{p,q}\right|}
{\max\left(1,\left|M_{p,q}\right|\right)}.
$$

This measures accuracy of the returned quadrature for the prescribed moments.
It is distinct from conditioning: a small error for one exact input does not
prove insensitivity to perturbations.

### Residual

The residual is a backward-error-style check of linear-system consistency. A
small residual means the computed weights satisfy the represented system
closely; it does not establish independent columns or a stable solution.
Residuals remain available in the complete flat diagnostics file.

The local matrix quantities are calculated by the example from the prescribed
population. They are not production QR or SVD factorization metadata.

## Conditioning Orders

### x_then_y — (0,1)

- reconstructs x first;
- places the close x-nodes in the parent quadrature; and
- exposes the weak system in conditional-moment diagnostics.

### y_then_x — (1,0)

- reconstructs y first;
- begins with the well-separated parent matrix; and
- encounters the weak x-reconstruction later in the recursive process.

The second order is not unconditionally better. It can report a later
Wheeler/QMOM failure instead of an earlier conditional-moment failure.

## Fixed And Adaptive Reconstructions

The fixed-count sweep uses node counts $(2,2)$, fixed branch counts, and

    1e-1, 1e-6, 1e-10, 1e-12, 1e-13, 2.5e-14, 2e-14, 0

It produces

$$
8\ \delta\text{ values}\times2\ \text{orders}\times5\ \text{policies}=80
$$

fixed-count rows.

Two separate adaptive references use $\delta=0$ with adaptive branch counts
enabled. Adaptive reduction occurs during Wheeler inversion, before the
conditional Vandermonde solve. It is not triggered by QR or SVD. The adaptive
rows represent lower-order quadratures and are not solver alternatives for the
fixed two-node formulation.

## Build And Run

Install QMOMKit, then configure this directory as a standalone installed
package consumer:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
    cmake --build build
    cmake --install build

    cmake -S examples/cqmom_2d_borderline_conditioning -B examples/cqmom_2d_borderline_conditioning/build -DCMAKE_PREFIX_PATH="$PWD/install"
    cmake --build examples/cqmom_2d_borderline_conditioning/build
    
    cd examples/cqmom_2d_borderline_conditioning/build
    ./cqmom_2d_borderline_conditioning_example

## Output Files

### output/report.md

This is the recommended starting point for interpretation. It explains
the formulation, defines the metrics, compares Dense/BP/QR/SVD by conditioning
order, lists the observed failed stages, and shows the lower-node cases. 

### output/diagnostics.txt

This is the complete 66-row dataset. It retains all four
requested fixed-count policies, the reported diagnostics, local matrix values,
residuals, moment errors, failure messages, and adaptive rows. It is
whitespace-delimited with a commented header and uses a dash for unavailable
values.

## How To Read The Report

1. Read the central solver-comparison message.
2. Compare the relative singular value and condition number at one $\delta$.
3. Compare Dense, BP, QR, and SVD at that same $\delta$ and conditioning order.
4. Use moment error for successful reconstruction accuracy and residual for
   linear-system consistency.
5. Read the failed-stage table to distinguish a conditional Vandermonde
   failure from a Wheeler/QMOM failure.
6. Compare fixed-count failure with the adaptive $\delta=0$ result as a
   lower-node reconstruction, not as a solver rescue.

The example does not provide a universal solver ranking. See
[the internal linear-algebra guide](../../docs/developer/linear_algebra.md) for
the broader solver and tolerance contract.

## Expected Behavior

Large separations should be successful for all direct policies. As $\delta$
decreases, $r_\sigma$ shrinks and $\kappa_2$ grows. Near the numerical
threshold, fixed-count runs can fail. At exact degeneracy, the requested
two-node representation is unsupported, while adaptive Wheeler reduction can
recover a lower-order representation.

The exact transition is platform- and backend-sensitive. The report describes
the observed run; its threshold values are not portable guarantees.

## Limitations

This is an artificial static four-atom population. The local diagnostic covers
only the known weighted $2\times2$ parent matrix and cannot expose production
QR/SVD factorization data. The example does not implement pseudoinverse,
truncation, regularization, minimum-norm recovery, automatic QR/SVD selection,
or adaptive node-count changes based on QR/SVD rank.
