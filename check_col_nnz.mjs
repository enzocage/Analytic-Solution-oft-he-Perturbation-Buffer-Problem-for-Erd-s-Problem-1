// Print the exact induced DAG on S_j (the support of col_j(H_s^{-1})) and its edge weights H[i][k]/H[i][i]!
import fs from 'fs';

// Reuse FastHNF and invertUpperTriangular from check_col_nnz.mjs
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
  toNumber() { return Number(this.p) / Number(this.q); }
  toString() { return this.q === 1n ? `${this.p}` : `${this.p}/${this.q}`; }
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

// Check the maximum out-degree of ANY vertex i inside S_j!
for (const [b, s] of [[3, 2], [3, 3], [3, 4], [5, 2], [5, 3], [7, 2], [7, 3]]) {
  const H = buildFastHNF(b, s);
  const Hinv = invertUpperTriangular(H);
  const r = H.length;
  let maxOutDegInSj = 0;
  let maxOutDegToNonJ = 0;
  for (let j = 0; j < r; j++) {
    const Sj = [];
    for (let i = 0; i <= j; i++) if (!Hinv[i][j].isZero()) Sj.push(i);
    const SjSet = new Set(Sj);
    for (const i of Sj) {
      if (i === j) continue;
      let deg = 0, degNonJ = 0;
      for (const k of Sj) {
        if (k > i && H[i][k] !== 0n) {
          deg++;
          if (k < j) degNonJ++;
        }
      }
      if (deg > maxOutDegInSj) maxOutDegInSj = deg;
      if (degNonJ > maxOutDegToNonJ) maxOutDegToNonJ = degNonJ;
    }
  }
  console.log(`b=${b}, s=${s}: maxOutDegInSj = ${maxOutDegInSj}, maxOutDegToNonJ (k < j in S_j) = ${maxOutDegToNonJ}`);
}
