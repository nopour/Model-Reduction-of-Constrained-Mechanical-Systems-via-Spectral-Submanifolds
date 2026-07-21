# Model Reduction of Constrained Mechanical Systems via Spectral Submanifolds

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![SSMTool](https://img.shields.io/badge/Powered_by-SSMTool_v2.2%2B-green.svg)](https://github.com/e3s-ethz/SSMTool)

A draft study on exact model reduction for mechanical systems subject to holonomic algebraic configuration constraints, building upon Spectral Submanifold (SSM) theory developed by Prof. George Haller's group (ETH Zürich).

---

## 📌 Executive Summary

Dynamical systems in engineering are frequently subject to algebraic configuration constraints alongside their governing differential equations, forming **Differential-Algebraic Equations (DAEs)**. Simulating high-dimensional constrained systems is computationally expensive and can obscure underlying nonlinear modal interactions.

This repository explores the mathematical reduction of constrained mechanical systems. By casting the DAEs into a first-order descriptor form, **Spectral Submanifold (SSM)** parameterization can process these constraints natively—without explicit index-reduction algorithms or state-dependent coordinate partitioning.

### Key Theoretical Aspects
* **Descriptor-Form Processing:** Operates directly on index-3 DAEs.
* **Internal Resonance Resolution:** Demonstrates why classical 2D invariant manifolds fail during $1:2$ internal resonance ($\omega_2 \approx 2\omega_1$) due to small divisors, and resolves the issue using a **4D Spectral Hyper-Manifold**.
* **Constraint Force Extraction:** Directly recovers dynamic Lagrange multipliers ($\lambda$) from the geometry of the manifold.

---

## 📐 Mathematical Formulation

### 1. Descriptor DAE System

The governing equations for a mechanical system with $n_c$ holonomic constraints $g(x) = 0$ are defined as:

$$M \ddot{x} + C \dot{x} + K x + f_{\text{nl}}(x,\dot{x}) + G(x)^T \mu = \epsilon f_{\text{ext}}(x,\dot{x},\Omega t)$$

$$g(x) = 0$$

By defining the augmented state vector $z = (x, \dot{x}, \mu)^T$, the system maps into a first-order descriptor form:

$$B \dot{z} = A z + F(z) + \epsilon F_{\text{ext}}(z, \Omega t)$$

where the global linear matrices $A$ and $B$ are defined as:

$$A = \begin{bmatrix} -K & 0 & -G_0^T \\ 0 & M & 0 \\ G_0 & 0 & 0 \end{bmatrix}, \quad B = \begin{bmatrix} C & M & 0 \\ M & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}$$

### 2. The Invariance Equation

The autonomous SSM is defined by an embedding function $W_0(p)$ and its reduced internal dynamics $R_0(p)$, which satisfy the invariance equation:

$$B D W_0(p) R_0(p) = A W_0(p) + F(W_0(p))$$

Solving this equation via multivariate Taylor series expansion yields a sequence of linear cohomological equations, allowing for analytical extraction of the manifold up to arbitrary high polynomial orders.

---

## 🗺️ Reduction Architecture
The logic flow from the full-order DAE to the low-dimensional ODE Reduced-Order Model (ROM) is visualized below:
```mermaid
graph TD
    A([Full DAE System<br/>Index-3]) --> B[Assemble Descriptor Matrices<br/>A, B, F]
    B --> C{Generalized Eigenvalue Problem<br/>Av = λBv}
    
    C -->|Infinite λ| D[Constraint Modes<br/>Excluded]
    C -->|Finite λ| E[Master Spectral Subspace<br/>Selected Modes]
    
    E --> F{Check Commensurate Frequencies}
    F -->|No Resonance| G[2D SSM Master Space<br/>p = p₁, p₁_bar]
    F -->|1:2 Resonance| H[4D SSM Hyper-Manifold<br/>p = p₁, p₁_bar, p₂, p₂_bar]
    
    G --> I[Solve Cohomological Equations]
    H --> I
    
    I --> J([Polynomial Reduced Dynamics<br/>ROM])
    I --> K([Phase Space Embedding W₀<br/>Recovers Displacements & Forces])
```
## 🧪 Benchmark: Hyperbolic Paraboloid Oscillator
To demonstrate the necessity of correct subspace selection, the framework is applied to a 3-DOF oscillator strictly confined to a hyperbolic paraboloid manifold:
The system is specifically tuned such that the principal curvature modes exhibit a 1:2 internal resonance (\omega_2 = 4.0, \omega_1 = 2.0).

> [!WARNING]
> **Topological Breakdown of 2D Manifolds**
> During 1:2 internal resonance, standard 2D manifolds suffer from "small divisor" singularities because the matrix [(2\lambda_1)B - A] becomes strictly singular. The trajectory physically escapes the 2D plane, necessitating the 4D SSM approach used in this repository.
> 
### Verification Metrics
Comparison of the 4D SSM prediction (\mathcal{O}(9) expansion) against a direct Index-1 stabilized numerical integration of the full-scale system reveals extreme precision:

| State Variable | Max Absolute Error | RMS Error | Relative L_2 Norm |
| :--- | :--- | :--- | :--- |
| **x_1 (Mode 1)** | 1.34 \times 10^{-6} | 7.19 \times 10^{-7} | 1.42 \times 10^{-3} |
| **x_2 (Resonant Mode 2)** | 4.68 \times 10^{-9} | 2.43 \times 10^{-9} | 1.75 \times 10^{-6} |
| **x_3 (Constrained Axis)** | 1.79 \times 10^{-7} | 5.14 \times 10^{-8} | 1.75 \times 10^{-5} |
| **\lambda (Constraint Force)** | 1.09 \times 10^{-5} | 1.12 \times 10^{-6} | 1.84 \times 10^{-5} |

## 📚 Acknowledgments & Literature
This study relies on the fundamental advancements in spectral submanifold theory and constrained mechanical reduction. Key literature includes:
 1. **Haller, G., & Ponsioen, S. (2016).** *Nonlinear normal modes and spectral submanifolds: existence, uniqueness and use in model reduction.* Nonlinear Dynamics, 86(3), 1493-1534.
 2. **Li, M., Jain, S., & Haller, G. (2023).** *Model reduction for constrained mechanical systems via spectral submanifolds.* Nonlinear Dynamics, 111, 8881-8911.
 3. **Shaw, S. W., & Pierre, C. (1993).** *Normal modes for non-linear vibratory systems.* Journal of Sound and Vibration, 164(1), 85-124.


Distributed under the MIT License.

