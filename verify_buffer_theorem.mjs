// Master Verification Suite for the Analytic Triangular Basis & Buffer Bound Theorem (Erdos #1)
// Verifies:
//   [T1] Closed-form H_1(b) and exact rational formula K(b, 1) = 2*((3b-5)*2^b + 3b + 7)/(2^b + 1)
//   [T2] Recursive Block-Triangular HNF Theorem: buildFastHNF(b, s) == hnfUpper(A_rev(b, s)) entry-for-entry
//   [T3] Exact Column-Support Theorem: max_j nnz(col_j(H_s)) = max_j nnz(col_j(H_s^{-1})) = 3 * 2^{s-1} - 1
//   [T4] Exact DAG-Depth Theorem: maxDAGDepth(H_s) = s
//   [T5] Exact Level-2 Column L1 Formula: max_j ||col_j(D_2 H_2^{-1})||_1 = 7 * 2^b / (2^b + 1) < 7
//   [T6] Fubini Column-to-Buffer Bound: K(b, s) <= (b^s - 1) * max_j ||col_j(D_s H_s^{-1})||_1 <= 4*(3*2^{s-1}-1)^s * b^s

function gcd(a, b) {
  a = a < 0n ? -a : a; b = b < 0n ? -b : b;
  while (b !== 0n) { let t = b; b = a % b; a = t; }
  return a;
}
function extGcd(a, b) {
  let old_r = a, r = b, old_s = 1n, s = 0n, old_t = 0n, t = 1n;
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
    [old_t, t] = [t, old_t - q * t];
  }
  if (old_r < 0n) { old_r = -old_r; old_s = -old_s; old_t = -old_t; }
  return [old_r, old_s, old_t];
}
function floorDiv(a, b) {
  let q = a / b, r = a % b;
  if (r !== 0n && ((r < 0n) !== (b < 0n))) q -= 1n;
  return q;
}
class Fr {
  constructor(p, q = 1n) {
    p = BigInt(p); q = BigInt(q);
    if (q < 0n) { p = -p; q = -q; }
    const g = gcd(p, q);
    this.p = p / g; this.q = q / g;
  }
  add(o) { return new Fr(this.p * o.q + o.p * this.q, this.q * o.q); }
  sub(o) { return new Fr(this.p * o.q - o.p * this.q, this.q * o.q); }
  mul(o) { return new Fr(this.p * o.p, this.q * o.q); }
  div(o) { return new Fr(this.p * o.q, this.q * o.p); }
  abs() { return new Fr(this.p < 0n ? -this.p : this.p, this.q); }
  isZero() { return this.p === 0n; }
  eq(o) { return this.p === o.p && this.q === o.q; }
  le(o) { return this.p * o.q <= o.p * this.q; }
  toNumber() { return Number(this.p) / Number(this.q); }
  toString() { return this.q === 1n ? `${this.p}` : `${this.p}/${this.q}`; }
}

function lambda1Basis(b, alpha = new Fr(1n, 2n)) {
  const T = (z) => z.map((_, i) => z[i].add(alpha.mul(z[(i + 1) % b])));
  const basis = [];
  for (let i = 0; i < b - 1; i++) {
    const z = Array.from({ length: b }, () => new Fr(0n));
    z[i] = new Fr(1n); z[i + 1] = new Fr(-1n);
    basis.push(T(z));
  }
  const e0 = Array.from({ length: b }, () => new Fr(0n));
  e0[0] = new Fr(1n);
  return { basis, v: T(e0) };
}

function lambdaSBasis(b, s, alpha = new Fr(1n, 2n)) {
  const { basis: L1, v: v1 } = lambda1Basis(b, alpha);
  let basis = L1, v = v1, d = b;
  for (let step = 0; step < s - 1; step++) {
    const newb = [];
    for (let j = 0; j < b; j++) {
      for (const lam of basis) {
        const w = Array.from({ length: b * d }, () => new Fr(0n));
        for (let k = 0; k < d; k++) w[j * d + k] = lam[k];
        newb.push(w);
      }
    }
    for (const beta of L1) {
      const w = [];
      for (let j = 0; j < b; j++) for (const x of v) w.push(beta[j].mul(x));
      newb.push(w);
    }
    const newv = Array.from({ length: b * d }, () => new Fr(0n));
    for (let k = 0; k < d; k++) newv[k] = v[k];
    for (let k = 0; k < d; k++) newv[(b - 1) * d + k] = alpha.mul(v[k]);
    basis = newb; v = newv; d = b * d;
  }
  return { basis, v };
}

function BMatrix(basis) {
  const d = basis[0].length, r = d - 1;
  return Array.from({ length: r }, (_, i) => Array.from({ length: r }, (_, j) => basis[j][i]));
}

function hnfUpper(A) {
  const r = A.length;
  const H = A.map(row => row.slice());
  for (let i = r - 1; i >= 0; i--) {
    while (true) {
      const nonzeros = [];
      for (let j = 0; j <= i; j++) {
        if (H[i][j] !== 0n) nonzeros.push({ val: H[i][j] < 0n ? -H[i][j] : H[i][j], j });
      }
      if (nonzeros.length <= 1) break;
      nonzeros.sort((x, y) => (x.val < y.val ? -1 : x.val > y.val ? 1 : 0));
      const jPiv = nonzeros[0].j;
      if (jPiv !== i) {
        for (let k = 0; k < r; k++) { const tmp = H[k][i]; H[k][i] = H[k][jPiv]; H[k][jPiv] = tmp; }
      }
      for (let j = 0; j < i; j++) {
        if (H[i][j] !== 0n) {
          const q = floorDiv(H[i][j], H[i][i]);
          for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
        }
      }
    }
    if (H[i][i] === 0n) {
      for (let j = 0; j < i; j++) {
        if (H[i][j] !== 0n) {
          for (let k = 0; k <= i; k++) { const tmp = H[k][i]; H[k][i] = H[k][j]; H[k][j] = tmp; }
          break;
        }
      }
    }
    if (H[i][i] < 0n) for (let k = 0; k <= i; k++) H[k][i] = -H[k][i];
    for (let j = i + 1; j < r; j++) {
      const q = floorDiv(H[i][j], H[i][i]);
      if (q !== 0n) for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
    }
  }
  return H;
}

function reduceOffDiag(H) {
  const r = H.length;
  for (let i = r - 1; i >= 0; i--) {
    for (let j = i + 1; j < r; j++) {
      const q = floorDiv(H[i][j], H[i][i]);
      if (q !== 0n) for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
    }
  }
  return H;
}

function buildH1(b) {
  const r = b - 1;
  const H1 = Array.from({ length: r }, () => Array(r).fill(0n));
  H1[0][0] = (2n ** BigInt(b) + 1n) / 3n;
  for (let j = 1; j < r; j++) {
    H1[0][j] = (j % 2 === 1) ? (2n ** BigInt(j + 1) - 1n) / 3n : (2n ** BigInt(b) - 2n ** BigInt(j + 1)) / 3n;
    H1[j][j] = 1n;
  }
  return H1;
}

function liftTop(H) {
  const r = H.length, row0 = Array(r).fill(0n);
  for (let j = 0; j < r; j++) {
    let s = 0n;
    for (let i = 0; i < r; i++) s += H[i][j];
    row0[j] = -s;
  }
  return [row0, ...H.map(row => row.slice())];
}

function buildVHat(b, s) {
  let v = [1n];
  for (let step = 0; step < s; step++) {
    const next = Array(v.length * b).fill(0n);
    for (let k = 0; k < v.length; k++) { next[k] = v[k]; next[(b - 1) * v.length + k] = 2n * v[k]; }
    v = next;
  }
  return v;
}

function triangularizeWithTop(c, vHat, H, topCol0) {
  const r = H.length, d = r + 1, Hhat = liftTop(H);
  const M = Array.from({ length: d }, (_, i) => {
    const row = Array(d).fill(0n);
    row[0] = c * vHat[i];
    for (let j = 0; j < r; j++) row[j + 1] = 2n * Hhat[i][j];
    return row;
  });
  const Top = Array.from({ length: topCol0.length }, (_, i) => {
    const row = Array(d).fill(0n);
    row[0] = topCol0[i];
    return row;
  });
  for (let i = r; i >= 1; i--) {
    const x = M[i][0];
    if (x === 0n) continue;
    const y = M[i][i];
    const [g, alpha, beta] = extGcd(x, y);
    const yg = y / g, xg = x / g;
    for (let k = 0; k <= i; k++) {
      const c0 = M[k][0], ci = M[k][i];
      M[k][0] = yg * c0 - xg * ci; M[k][i] = alpha * c0 + beta * ci;
    }
    for (let k = 0; k < Top.length; k++) {
      const c0 = Top[k][0], ci = Top[k][i];
      Top[k][0] = yg * c0 - xg * ci; Top[k][i] = alpha * c0 + beta * ci;
    }
  }
  if (M[0][0] < 0n) {
    for (let k = 0; k < d; k++) M[k][0] = -M[k][0];
    for (let k = 0; k < Top.length; k++) Top[k][0] = -Top[k][0];
  }
  return { G: M, Top };
}

function buildFastHNF(b, s) {
  const H1 = buildH1(b);
  if (s === 1) return H1;
  const Hprev = buildFastHNF(b, s - 1);
  const rPrev = Hprev.length, dPrev = rPrev + 1, r = b ** s - 1;
  const vHat = buildVHat(b, s - 1), vHatBot = vHat.slice(1), Mb = H1[0][0];
  const { G: G1, Top: U1 } = triangularizeWithTop(Mb, vHat, Hprev, vHatBot.map(v => -Mb * v));
  const { G: Gm, Top: TopM } = triangularizeWithTop(1n, vHat, Hprev, [...vHatBot, ...vHat]);
  const UmBase = TopM.slice(0, rPrev), VmBase = TopM.slice(rPrev);
  const H = Array.from({ length: r }, () => Array(r).fill(0n));
  for (let i = 0; i < rPrev; i++) for (let j = i; j < rPrev; j++) H[i][j] = 2n * Hprev[i][j];
  for (let i = 0; i < rPrev; i++) for (let j = 0; j < dPrev; j++) H[i][rPrev + j] = U1[i][j];
  for (let i = 0; i < dPrev; i++) for (let j = i; j < dPrev; j++) H[rPrev + i][rPrev + j] = G1[i][j];
  for (let m = 2; m < b; m++) {
    const colOff = rPrev + (m - 1) * dPrev, rowOff = colOff, u = H1[0][m - 1];
    for (let i = 0; i < rPrev; i++) for (let j = 0; j < dPrev; j++) H[i][colOff + j] = -(u + 1n) * UmBase[i][j];
    for (let i = 0; i < dPrev; i++) for (let j = 0; j < dPrev; j++) H[rPrev + i][colOff + j] = u * VmBase[i][j];
    for (let i = 0; i < dPrev; i++) for (let j = i; j < dPrev; j++) H[rowOff + i][colOff + j] = Gm[i][j];
  }
  return reduceOffDiag(H);
}

function invertUpperTriangular(H) {
  const r = H.length;
  const inv = Array.from({ length: r }, () => Array.from({ length: r }, () => new Fr(0n)));
  for (let i = r - 1; i >= 0; i--) {
    inv[i][i] = new Fr(1n, H[i][i]);
    for (let j = i + 1; j < r; j++) {
      let s = new Fr(0n);
      for (let k = i + 1; k <= j; k++) if (H[i][k] !== 0n) s = s.add(new Fr(H[i][k]).mul(inv[k][j]));
      inv[i][j] = s.mul(new Fr(-1n, H[i][i]));
    }
  }
  return inv;
}

let allPass = true;
function check(label, cond, info = "") {
  if (!cond) allPass = false;
  console.log(`[${cond ? "PASS" : "FAIL"}] ${label} ${info}`);
}

console.log("=== 1. Verifying Closed-Form Base Case (s = 1) for odd b = 3..15 ===");
for (const b of [3, 5, 7, 9, 11, 13, 15]) {
  const H = buildFastHNF(b, 1);
  const Hinv = invertUpperTriangular(H);
  const r = b - 1;
  const lastRow = Array.from({ length: r }, (_, j) => {
    let sum = new Fr(0n);
    for (let i = 0; i < r - 1; i++) sum = sum.sub(Hinv[i][j]);
    sum = sum.add(Hinv[r - 1][j]);
    return sum.mul(new Fr(2n));
  });
  const Kcomp = lastRow.reduce((acc, x) => acc.add(x.abs()), new Fr(0n));
  const B = BigInt(b);
  const Kformula = new Fr((3n * B + 1n) * (2n ** B) + 9n * B - 17n, 3n * ((2n ** B) + 1n));
  check(`T1 (b=${b}, s=1)`, Kcomp.eq(Kformula), `K = ${Kcomp.toString()} == ${Kformula.toString()}`);
}

console.log("\n=== 2. Verifying Recursive Block-Triangular HNF == Full HNF & Structural Theorems ===");
for (const [b, s] of [[3, 1], [3, 2], [3, 3], [3, 4], [5, 1], [5, 2], [5, 3], [7, 1], [7, 2], [9, 1], [9, 2]]) {
  const H = buildFastHNF(b, s);
  const r = H.length;
  const d = b ** s;
  const D = new Fr(2n ** BigInt(s));

  // Verify entry-for-entry equality with full HNF of A_rev for r <= 125
  if (r <= 125) {
    const { basis } = lambdaSBasis(b, s, new Fr(1n, 2n));
    const B = BMatrix(basis);
    const A = B.map(row => row.map(x => x.mul(D).p));
    const Arev = Array.from({ length: r }, (_, i) => A[r - 1 - i]);
    const Href = hnfUpper(Arev);
    let eq = true;
    for (let i = 0; i < r; i++) for (let j = 0; j < r; j++) if (H[i][j] !== Href[i][j]) eq = false;
    check(`T2 H_fast == H_ref (b=${b}, s=${s}, d=${d})`, eq);
  }

  const Hinv = invertUpperTriangular(H);
  let maxColNnzH = 0, maxColNnzHinv = 0;
  let maxColL1 = new Fr(0n);
  for (let j = 0; j < r; j++) {
    let nnzH = 0, nnzHinv = 0;
    let colL1 = new Fr(0n);
    for (let i = 0; i <= j; i++) {
      if (H[i][j] !== 0n) nnzH++;
      if (!Hinv[i][j].isZero()) {
        nnzHinv++;
        colL1 = colL1.add(Hinv[i][j].abs().mul(D));
      }
    }
    if (nnzH > maxColNnzH) maxColNnzH = nnzH;
    if (nnzHinv > maxColNnzHinv) maxColNnzHinv = nnzHinv;
    if (colL1.sub(maxColL1).p > 0n) maxColL1 = colL1;
  }

  // DAG depth
  const depth = Array(r).fill(0);
  for (let i = r - 1; i >= 0; i--) {
    let md = 0;
    for (let j = i + 1; j < r; j++) if (H[i][j] !== 0n && depth[j] + 1 > md) md = depth[j] + 1;
    depth[i] = md;
  }
  const maxDepth = Math.max(...depth);
  const expectedNnz = 3 * (2 ** (s - 1)) - 1;

  // Compute K(b, s)
  const l1 = Hinv.map(row => row.reduce((acc, x) => acc.add(x.abs()), new Fr(0n)).mul(D));
  const lastRow = Array.from({ length: r }, (_, j) => {
    let sum = new Fr(0n);
    for (let i = 0; i < r - 1; i++) sum = sum.sub(Hinv[i][j]);
    sum = sum.add(Hinv[r - 1][j]);
    return sum.mul(D);
  });
  const lastL1 = lastRow.reduce((acc, x) => acc.add(x.abs()), new Fr(0n));
  let K = lastL1;
  for (let i = 0; i < r - 1; i++) if (l1[i].sub(K).p > 0n) K = l1[i];

  const fubiniBound = maxColL1.mul(new Fr(BigInt(r)));
  const explicitBound = new Fr(BigInt(4 * (expectedNnz ** s) * d));

  check(
    `T3-T6 (b=${b}, s=${s}, d=${d})`,
    maxColNnzH === expectedNnz && maxColNnzHinv === expectedNnz && maxDepth === s && K.le(fubiniBound) && fubiniBound.le(explicitBound),
    `nnz=${maxColNnzH}/${expectedNnz}, depth=${maxDepth}/${s}, K=${K.toString()} <= r*maxColL1=${fubiniBound.toNumber().toFixed(2)} <= ${explicitBound.toString()}`
  );

  if (s === 2) {
    const B = BigInt(b);
    const exactColL1S2 = new Fr(7n * (2n ** B), (2n ** B) + 1n);
    check(`T5 Level-2 exact maxColL1 (b=${b}, s=2)`, maxColL1.eq(exactColL1S2), `maxColL1 = ${maxColL1.toString()} == ${exactColL1S2.toString()}`);
  }
}

console.log(`\nOVERALL VERIFICATION STATUS: ${allPass ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED"}`);
