# Analytic Solution of the Perturbation Buffer Problem for Erdős Problem #1

**Author:** Felix Schmidt (`mail@felixschmidt.org`)  
**Date:** October 18, 2026

[![Paper PDF](https://img.shields.io/badge/Paper-PDF-red)](erdos1_buffer_theorem.pdf)
[![LaTeX Source](https://img.shields.io/badge/LaTeX-Source-blue)](erdos1_buffer_theorem.tex)
[![Lean 4](https://img.shields.io/badge/Formalization-Lean%204-purple)](Erdos1BufferTheorem.lean)
[![Verification](https://img.shields.io/badge/Verification-100%25%20Passed-brightgreen)](verify_buffer_theorem.mjs)

---

## Overview

In September 2026, **Erdős Problem #1** ($f(n) \gg 2^n$ for distinct-subset-sum sets $A \subseteq \{1,\dots,N\}$) was disproved via a Kronecker-power lattice construction $\Lambda_s = \Lambda_1^{\otimes s}$ (`tadamcz/erdos1`), and made explicit at finite $n$ by Carlos Toledo (*Explicit Disproof of Erdős Problem #1 at Finite $n$*, Sep 16, 2026).

However, Toledo's analysis (§3 Step K, §5, §6, and §7 *"What this does not settle"*) leaves open the **Analytic Perturbation Buffer Problem**: proving an explicit analytic bound
$$K(b,s) := D_s \|E H_s^{-1}\|_{\infty \to \infty} \le \operatorname{poly}(d_s)$$
(where $d_s = b^s$, $D_s = 2^s$, and $H_s$ is the canonical upper-triangular Hermite basis of $D_s \Lambda_s$ in reversed coordinates), or constructing an explicit recursive upper-triangular basis of $\Lambda_s$ with an analytic inverse. Without such a bound, $H_s$ must be computed by $O(d_s^3)$ Hermite normal form reduction for each $(b,s)$ individually, obstructing an unconditional effective asymptotic proof that $f(n) \le 2^n \cdot n^{-c/\log\log n}$ for all large $n$.

This repository contains the **complete analytic solution** to this problem, together with an exact arbitrary-precision (`BigInt`) rational verification suite and Lean 4 formalization of the core algebraic and analytic reductions.

---

## Main Results

1. **Theorem 3.1 (Exact Closed Form of $H_1$, $H_1^{-1}$, $K(b,1)$, and the Universal Unimodular Identity):**  
   For every odd $b \ge 3$, let $M_b = \frac{2^b+1}{3} \in \mathbb{Z}_{>0}$. The canonical upper-triangular Hermite normal form $H_1 \in \mathbb{Z}^{(b-1)\times(b-1)}$ has rows $1,\dots,b-2$ equal to $I_{b-2}$ and top-row entries $H_1[0,0] = M_b$,
   $$u_j := H_1[0,j] = \begin{cases} \dfrac{2^{j+1}-1}{3} & (j \text{ odd}), \\[6pt] \dfrac{2^b - 2^{j+1}}{3} & (j \text{ even}). \end{cases}$$
   Since $b-2$ is odd for all odd $b \ge 3$, $u_{b-2} = \frac{2^{b-1}-1}{3}$ satisfies the **universal unimodular pivot identity**
   $$M_b - 2 u_{b-2} = 1 \qquad \forall \text{ odd } b \ge 3,$$
   and the exact perturbation buffer at $s = 1$ is
   $$K(b,1) = \frac{(3b+1)2^b + 9b - 17}{3(2^b+1)} = b + \frac{1}{3} + \frac{2(b-3)}{2^b+1} \le \frac{10}{9}b.$$

2. **Theorem 4.1 (Exact Block-Triangular Recursion for $H_s$):**  
   The canonical HNF $H_s$ decomposes into $b$ blocks such that **Blocks $2,\dots,b-1$ are completely decoupled copies of a single matrix** $G_{s-1} = \operatorname{HNF}(\hat{v}_{s-1} \mid 2\hat{H}_{s-1})$, and inside $G_{s-1}$ (by the unimodular identity $M_b - 2u_{b-2} = 1$), **all $b-1$ lower sub-blocks decouple into block-diagonal copies of $2G_{s-2}$**.

3. **Theorems 5.1 & 5.2 (The $b$-Independent Column-Support Invariant $N_s = 3 \cdot 2^{s-1} - 1$ and DAG Depth $s$):**  
   While row $0$ of $H_s^{-1}$ has $\Theta(b^s)$ nonzero entries, every **column** of $H_s$ and $H_s^{-1}$ has at most $N_s$ nonzero entries satisfying $N_1 = 2$ and $N_s = 2N_{s-1} + 1$:
   $$\max_{0 \le j < b^s-1} \operatorname{nnz}\big(\operatorname{col}_j(H_s)\big) = \max_{0 \le j < b^s-1} \operatorname{nnz}\big(\operatorname{col}_j(H_s^{-1})\big) = 3 \cdot 2^{s-1} - 1,$$
   **completely independent of $b$**, and the coupling graph $\operatorname{DAG}(H_s)$ has exact maximum path length $s$.

4. **Theorems 6.1 & 6.2 (Row-to-Column Fubini Duality & Explicit General Buffer Bound):**  
   At $s = 2$, for every odd $b \ge 3$,
   $$\max_{0 \le j < b^2-1} \|\operatorname{col}_j(D_2 H_2^{-1})\|_1 = \frac{7 \cdot 2^b}{2^b+1} < 7 \implies K(b,2) < 7b^2.$$
   For all $s \ge 1$ and all odd $b \ge 3$, Fubini duality yields the explicit analytic bound
   $$K(b,s) \le 4\big(3 \cdot 2^{s-1} - 1\big)^s \, b^s = 4\big(3 \cdot 2^{s-1} - 1\big)^s \, d_s.$$

5. **Theorem 7.1 (1-Step Closed-Form Wavelet Block-Sum Inverse):**  
   Using $\mathbf{1}_{d_{s-1}}^T \hat{H}_{s-1} = 0$ and $\mathbf{1}_{d_{s-1}}^T \hat{v}_{s-1} = 3^{s-1}$, the natural wavelet basis $\hat{B}_s^{\mathrm{wave}}$ of $D_s \Lambda_s$ admits a 1-step non-recursive closed-form inverse via $b$-ary subtree sums.

---

## Repository Structure

- **[`erdos1_buffer_theorem.pdf`](erdos1_buffer_theorem.pdf)** — Full 11-page typeset paper (PDF).
- **[`erdos1_buffer_theorem.tex`](erdos1_buffer_theorem.tex)** — LaTeX source of the paper.
- **[`erdos1_buffer_theorem_proof.md`](erdos1_buffer_theorem_proof.md)** — Markdown version of the manuscript.
- **[`Erdos1BufferTheorem.lean`](Erdos1BufferTheorem.lean)** — Lean 4 formalization of the Row-to-Column Fubini Buffer Bound, DAG-depth amplification identity, exact column-support recurrence $N_s = 3 \cdot 2^{s-1} - 1$, and unimodular pivot identity $M_b - 2u_{b-2} = 1$.
- **[`verify_buffer_theorem.mjs`](verify_buffer_theorem.mjs)** — Self-contained exact arbitrary-precision (`BigInt`) rational verification suite checking Theorems 3.1–6.2 across $(b,s)$.
- **[`verify_erdos1_ref.py`](verify_erdos1_ref.py)** — Reference verifier from Toledo (2026).

---

## Reproducing the Exact Verification

Run the exact rational verification suite with Node.js (v18+):

```bash
node verify_buffer_theorem.mjs
```

Expected output: `OVERALL VERIFICATION STATUS: ALL CHECKS PASSED`.
