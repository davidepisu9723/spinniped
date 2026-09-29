# Chapter 3 -- Single-degree-of-freedom vibration

## Table of contents

- [Chapter content](#chapter-content)

## Chapter content

Let a point mass $m$ move in coordinate $q(t)$ against a linear spring $k$.
With no applied force and no damping, Newton's law gives

$$
m\ddot q+kq=0.
$$

Trying harmonic motion,

$$
q(t)=Qe^{i\omega t},
$$

gives $\ddot q=-\omega^2q$ and hence

$$
(k-m\omega^2)Q=0.
$$

A nonzero amplitude $Q$ requires the first term to be zero

$$
\omega_n=\sqrt{\frac{k}{m}},
\qquad
f_n=\frac{1}{2\pi}\sqrt{\frac{k}{m}}.
$$

This is the smallest eigenvalue problem. A continuous shaft follows the same
logic, but its displacement depends on both $z$ and $t$, so it has infinitely
many natural frequencies and mode shapes.
