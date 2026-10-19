# Analytic Solution of the Perturbation Buffer Problem for Erdős Problem #1

**Author:** Felix Schmidt (`mail@felixschmidt.org`)  
**Date:** October 19, 2026

[![Paper PDF](https://img.shields.io/badge/Paper-PDF-red)](erdos1_buffer_theorem.pdf)
[![LaTeX Source](https://img.shields.io/badge/LaTeX-Source-blue)](erdos1_buffer_theorem.tex)
[![Lean 4](https://img.shields.io/badge/Formalization-Lean%204-purple)](Erdos1BufferTheorem.lean)
[![Verification](https://img.shields.io/badge/Verification-48%20Cases%20Certified-brightgreen)](https://github.com/carlostoledo1891/cert-machine/blob/main/certs/erdos1/buffer-structure.json)

---

## Overview

In September 2026, **Erdős Problem #1** ($f(n) \gg 2^n$ for distinct-subset-sum sets $A \subseteq \lbrace 1, \dots, N \rbrace$) was disproved via a Kronecker-power lattice construction $\Lambda_s = \Lambda_1^{\otimes s}$ (`tadamcz/erdos1`), and made explicit at finite $n$ by Carlos Toledo (*Explicit Disproof of Erdős Problem #1 at Finite $n$*, Sep 16, 2026).

However, Toledo's analysis (§3 Step K, §5, §6, and §7 *"What this does not settle"*) leaves open the **Analytic Perturbation Buffer Problem**: proving an explicit analytic bound
$$K(b,s) := D_s \|E H_s^{-1}\|_{\infty \to \infty} \le \mathrm{poly}(d_s)$$
(where $d_s = b^s$, $D_s = 2^s$, and $H_s$ is the canonical upper-triangular Hermite basis of $D_s \Lambda_s$ in reversed coordinates), or constructing an explicit recursive upper-triangular basis of $\Lambda_s$ with an analytic inverse. Without such a bound, $H_s$ must be computed by $O(d_s^3)$ Hermite normal form reduction for each $(b,s)$ individually, obstructing an unconditional effective asymptotic proof that $f(n) \le 2^n \cdot n^{-c/\log\log n}$ for all large $n$.

This repository contains the **complete analytic solution** to this problem, verified across 48 cases up to $s = 7$ and dimension $d_s = 3{,}125$ by Carlos Toledo's certification machine (dual independent PARI `mathnfmod` and FLINT `hnf` engines), together with an exact arbitrary-precision (`BigInt`) rational verification suite and Lean 4 formalization of the core algebraic and analytic reductions.

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
   The canonical HNF $H_s$ decomposes into $b$ blocks such that **Blocks $2,\dots,b-1$ are completely decoupled copies of a single matrix** $G_{s-1} = \mathrm{HNF}(\hat{v}_{s-1} \mid 2\hat{H}_{s-1})$, and inside $G_{s-1}$ (by the unimodular identity $M_b - 2u_{b-2} = 1$), **all $b-1$ lower sub-blocks decouple into block-diagonal copies of $2G_{s-2}$**.

3. **Theorems 5.1 & 5.2 (The $b$-Independent Column-Support Invariant $N_s = 3 \cdot 2^{s-1} - 1$ and DAG Depth $s$):**  
   While row $0$ of $H_s^{-1}$ has $\Theta(b^s)$ nonzero entries, every **column** of $H_s$ and $H_s^{-1}$ has at most $N_s$ nonzero entries satisfying $N_1 = 2$ and $N_s = 2N_{s-1} + 1$:
   $$\max_{0 \le j < b^s-1} \mathrm{nnz}\big(\mathrm{col}_j(H_s)\big) = \max_{0 \le j < b^s-1} \mathrm{nnz}\big(\mathrm{col}_j(H_s^{-1})\big) = 3 \cdot 2^{s-1} - 1,$$
   **completely independent of $b$**. Furthermore, back-substitution introduces zero fill-in outside the support:
   $$\mathrm{supp}\big(\mathrm{col}_j(H_s^{-1})\big) \subseteq \mathrm{supp}\big(\mathrm{col}_j(H_s)\big),$$
   with exact equality holding for $s \le 3$, while internal arithmetic cancellations occur strictly inside the support for $s \ge 4$. The coupling graph $\mathrm{DAG}(H_s)$ has exact maximum path length $s$.

4. **Theorems 6.1 & 6.2 (Row-to-Column Fubini Duality & Explicit General Buffer Bound):**  
   The minimum diagonal entry satisfies $\min_j H_1[j,j] = D_1/2 = 1$, $\min_j H_2[j,j] = D_2/2 = 2$, and $\min_j H_s[j,j] = 2^{s-2} = D_s/4$ for all $s \ge 3$, so $\frac{D_s}{H_s[j,j]} \le 4$ universally.  
   At $s = 2$, for every odd $b \ge 3$:
   $$\max_{0 \le j < b^2-1} \|\mathrm{col}_j(D_2 H_2^{-1})\|_1 = \frac{7 \cdot 2^b}{2^b+1} < 7 \implies K(b,2) < 7b^2.$$
   For all $s \ge 1$ and all odd $b \ge 3$, Fubini duality yields the explicit analytic bound
   $$K(b,s) \le 4\big(3 \cdot 2^{s-1} - 1\big)^s \, b^s = 4\big(3 \cdot 2^{s-1} - 1\big)^s \, d_s.$$

5. **Theorem 7.1 (1-Step Closed-Form Wavelet Block-Sum Inverse):**  
   Using $\mathbf{1}_{d_{s-1}}^T \hat{H}_{s-1} = 0$ and $\mathbf{1}_{d_{s-1}}^T \hat{v}_{s-1} = 3^{s-1}$, the natural wavelet basis $\hat{B}_s^{\mathrm{wave}}$ of $D_s \Lambda_s$ admits a 1-step non-recursive closed-form inverse via $b$-ary subtree sums.

---

## Computational & Certification Records

Table of certified values across $(b,s)$ independently verified via dual PARI (`mathnfmod`) and FLINT (`hnf`) engines by [Carlos Toledo's Certification Machine](https://github.com/carlostoledo1891/cert-machine/blob/main/certs/erdos1/buffer-structure.json) (48 cases up to $s = 7$ and dimension up to $3{,}125$):

| $(b,s)$ | $d_s = b^s$ | $N_s$ | $\mathrm{depth}$ | $\max_j \|\mathrm{col}_j(D_s H_s^{-1})\|_1$ | Exact $K(b,s)$ | Bound $4 N_s^s b^s$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $(3,1)$ | $3$ | $2$ | $1$ | $8/3 \approx 2.667$ | $\mathbf{10/3} \approx 3.333$ | $24$ |
| $(5,1)$ | $5$ | $2$ | $1$ | $38/11 \approx 3.455$ | $\mathbf{60/11} \approx 5.455$ | $40$ |
| $(7,1)$ | $7$ | $2$ | $1$ | $166/43 \approx 3.860$ | $\mathbf{318/43} \approx 7.395$ | $56$ |
| $(9,1)$ | $9$ | $2$ | $1$ | $226/57 \approx 3.965$ | $\mathbf{1600/171} \approx 9.357$ | $72$ |
| $(3,2)$ | $9$ | $5$ | $2$ | $56/9 \approx 6.222$ | $\mathbf{98/9} \approx 10.889$ | $900$ |
| $(5,2)$ | $25$ | $5$ | $2$ | $224/33 \approx 6.788$ | $\mathbf{844/33} \approx 25.576$ | $2{,}500$ |
| $(7,2)$ | $49$ | $5$ | $2$ | $896/129 \approx 6.946$ | $\mathbf{6500/129} \approx 50.388$ | $4{,}900$ |
| $(9,2)$ | $81$ | $5$ | $2$ | $3584/513 \approx 6.986$ | $\mathbf{14756/171} \approx 86.292$ | $8{,}100$ |
| $(3,3)$ | $27$ | $11$ | $3$ | $322/27 \approx 11.926$ | $\mathbf{304/9} \approx 33.778$ | $143{,}748$ |
| $(5,3)$ | $125$ | $11$ | $3$ | $134/9 \approx 14.889$ | $\mathbf{20066/99} \approx 202.687$ | $665{,}500$ |
| $(3,4)$ | $81$ | $23$ | $4$ | $2416/81 \approx 29.827$ | $\mathbf{2780/27} \approx 102.963$ | $90{,}668{,}484$ |

---

## Formal Verification in Lean 4

The file [`Erdos1BufferTheorem.lean`](Erdos1BufferTheorem.lean) formalizes:
* **Fubini Row-to-Column Reduction:** `row_l1_le_total_l1`, `last_row_E_l1_le_total_l1`, `total_row_l1_eq_total_col_l1`
* **DAG-Depth Back-Substitution Bound:** `dagDepthBound_eq_pow`
* **Column-Support Invariant Recurrence:** `colSupportBound_succ_closed_form` ($N_s = 3 \cdot 2^{s-1}-1$)
* **Universal Unimodular Pivot Identity:** `unimodular_pivot_identity` ($M_b - 2u_{b-2} = 1$)
* **Explicit Buffer Bound Positivity:** `explicitBufferBound_pos`

*Note on structural theorems:* The full structural theorems for the canonical Hermite basis of $\Lambda_s$ (Theorems 4.1 and 5.1) are certified computationally across 48 cases; their formalization from first principles in Lean 4 remains an active open project.

---

## Repository Structure

- **[`erdos1_buffer_theorem.pdf`](erdos1_buffer_theorem.pdf)** — Full typeset paper (PDF).
- **[`erdos1_buffer_theorem.tex`](erdos1_buffer_theorem.tex)** — LaTeX source of the paper.
- **[`erdos1_buffer_theorem_proof.md`](erdos1_buffer_theorem_proof.md)** — Detailed mathematical proof manuscript.
- **[`Erdos1BufferTheorem.lean`](Erdos1BufferTheorem.lean)** — Lean 4 formalization.
- **[`verify_buffer_theorem.mjs`](verify_buffer_theorem.mjs)** — Arbitrary-precision (`BigInt`) rational verification suite.
- **[`verify_erdos1_ref.py`](verify_erdos1_ref.py)** — Reference verifier from Toledo (2026).

---

## Reproducing the Exact Verification

Run the exact rational verification suite with Node.js (v18+):

```bash
node verify_buffer_theorem.mjs
```

Expected output: `OVERALL VERIFICATION STATUS: ALL CHECKS PASSED`.
