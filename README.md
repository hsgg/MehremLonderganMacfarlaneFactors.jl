# MehremLonderganMacfarlaneFactors.jl
<!---
[![Build Status](https://github.com/hsgg/MehremLonderganMacfarlaneFactors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hsgg/MehremLonderganMacfarlaneFactors.jl/actions/workflows/CI.yml?query=branch%3Amain)
-->

This [Julia](https://julialang.org/) package implements the *Mehrem-Londergan-Macfarlane* factor as derived in
[Mehrem *etal.* (1991)](https://ui.adsabs.harvard.edu/abs/1991JPhA...24.1435M/abstract), published [here](https://iopscience.iop.org/article/10.1088/0305-4470/24/7/018).
It is defined as the integral over three spherical Bessel functions,
```math
M_{\ell_1 \ell_2 \ell_3}(k_1, k_2, k_3) =
\int_0^\infty dx\,x^2
\,j_{\ell_1}(k_1x)
\,j_{\ell_2}(k_2x)
\,j_{\ell_3}(k_3x)
\,.
```
We implement this factor multiplied by a Wigner $3j$-symbol, as written in
[Mehrem & Hohenegger (2010)](https://arxiv.org/abs/1006.2108),
```math
\begin{align}
\begin{pmatrix} \ell_1 & \ell_2 & \ell_3 \\ 0 & 0 & 0 \end{pmatrix}
&M_{\ell_1 \ell_2 \ell_3}(k_1, k_2, k_3)
=
\frac{\pi\beta(\Delta)}{4 k_1 k_2 k_3} \,i^{\ell_1 + \ell_2 - \ell_3}
\\
&\times
\sqrt{2\ell_3 + 1}
\,\left(\frac{k_1}{k_3}\right)^{\ell_3}
\sum_{L=0}^{\ell_3}
\sqrt{\begin{pmatrix} 2\ell_3 \\ 2L \end{pmatrix}}
\,\left(\frac{k_2}{k_1}\right)^{L}
\sum_l (2l+1)
\begin{pmatrix} \ell_1 & \ell_3 - L & l \\ 0 & 0 & 0 \end{pmatrix}
\\
&\times
\begin{pmatrix} \ell_2 & L & l \\ 0 & 0 & 0 \end{pmatrix}
\begin{Bmatrix} \ell_1 & \ell_2 & \ell_3 \\ L & \ell_3 - L & l \end{Bmatrix}
P_l(\Delta)\,,
\end{align}
```
where
```math
\Delta = \frac{k_1^2 + k_2^2 - k_3^2}{2 k_1 k_2}\,,
```
and
```math
\beta(\Delta) = \theta(1 - \Delta)\,\theta(1 + \Delta)\,,
```
with $\theta$ the Heaviside function in half-maximum convention.


## Installation

To install this Julia package, press `]` and run
```julia
pkg> add https://github.com/hsgg/MehremLonderganMacfarlaneFactors.jl
```

If you have access to [WATCosmologyJuliaRegistry](https://github.com/Wide-Angle-Team/WATCosmologyJuliaRegistry),
then you can install it as usual: press `]` and run
```julia
pkg> add MehremLonderganMacfarlaneFactors
```


## Usage

Two functions are exported. First, there is the integral directly, with call signature
```julia
mehremlonderganmacfarlanefactor(l1, l2, l3, k1, k2, k3)
```
which will return $M_{\ell_1 \ell_2 \ell_3}(k_1, k_2, k_3)$.

Second, there is
```julia
mehremlonderganmacfarlanefactor3j(l1, l2, l3, k1, k2, k3)
```
which will return

$$\begin{pmatrix} \ell_1 & \ell_2 & \ell_3 \\\ 0 & 0 & 0 \end{pmatrix}
M_{\ell_1 \ell_2 \ell_3}(k_1, k_2, k_3)$$

and so it won't blow up in cases where the triangle condition on the $\ell_i$ is violated.


## TODO

- Choose $\ell_3$ to be the smallest angular momentum.
- Choose $k_3 > 0.1 min(k_1, k_2)$.
- Provide convenience function for choosing $\ell_3$ when it can be chosen freely.
- Add docstrings.
