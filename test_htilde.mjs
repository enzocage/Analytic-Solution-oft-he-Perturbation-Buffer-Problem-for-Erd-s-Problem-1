// Test the pure recursive upper-triangular basis Htilde_s (without cross-block Euclidean reduction)
// and compare K(Htilde_s) and D_s * Htilde_s^{-1} with H_s!

function gcd(a, b) {
  a = a < 0n ? -a : a; b = b < 0n ? -b : b;
  while (b !== 0n) { let t = b; b = a % b; a = t; }
  return a;
}

// Extended GCD: returns [g, alpha, beta] such that alpha*a + beta*b = g >= 0
function extGcd(a, b) {
  let old_r = a, r = b;
  let old_s = 1n, s = 0n;
  let old_t = 0n, t = 1n;
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
    [old_t, t] = [t, old_t - q * t];
  }
  if (old_r < 0n) {
    old_r = -old_r; old_s = -old_s; old_t = -old_t;
  }
  return [old_r, old_s, old_t];
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

// Closed-form H_1(b) for any odd b >= 3
function buildH1(b) {
  const r = b - 1;
  const H1 = Array.from({ length: r }, () => Array(r).fill(0n));
  const Mb = (2n ** BigInt(b) + 1n) / 3n;
  H1[0][0] = Mb;
  for (let j = 1; j < r; j++) {
    if (j % 2 === 1) {
      H1[0][j] = (2n ** BigInt(j + 1) - 1n) / 3n;
    } else {
      H1[0][j] = (2n ** BigInt(b) - 2n ** BigInt(j + 1)) / 3n;
    }
    H1[j][j] = 1n;
  }
  return H1;
}

// Lift an r x r matrix H to (r+1) x r matrix Hhat by prepending row 0 = -1^T H
function liftTop(H) {
  const r = H.length;
  const row0 = Array(r).fill(0n);
  for (let j = 0; j < r; j++) {
    let s = 0n;
    for (let i = 0; i < r; i++) s += H[i][j];
    row0[j] = -s;
  }
  return [row0, ...H.map(row => row.slice())];
}

// Build vHat_s = (1, 0, ..., 0, 2)^{\otimes s} of length b^s
function buildVHat(b, s) {
  let v = [1n];
  for (let step = 0; step < s; step++) {
    const next = Array(v.length * b).fill(0n);
    for (let k = 0; k < v.length; k++) {
      next[k] = v[k];
      next[(b - 1) * v.length + k] = 2n * v[k];
    }
    v = next;
  }
  return v;
}

// Triangularize M = [ c * vHat | 2 * Hhat ] where H is r x r upper-triangular
// Returns { G, steps } where G is (r+1) x (r+1) upper-triangular
function triangularizeBlock(c, vHat, H) {
  const r = H.length;
  const d = r + 1;
  const Hhat = liftTop(H);
  // M has columns: col 0 = c * vHat, cols 1..r = 2 * Hhat
  const M = Array.from({ length: d }, (_, i) => {
    const row = Array(d).fill(0n);
    row[0] = c * vHat[i];
    for (let j = 0; j < r; j++) row[j + 1] = 2n * Hhat[i][j];
    return row;
  });
  let steps = 0;
  for (let i = r; i >= 1; i--) {
    const x = M[i][0];
    if (x === 0n) continue;
    steps++;
    const y = M[i][i];
    const [g, alpha, beta] = extGcd(x, y);
    // new_col_0 = (y/g)*col_0 - (x/g)*col_i
    // new_col_i = alpha*col_0 + beta*col_i
    const yg = y / g, xg = x / g;
    for (let k = 0; k <= i; k++) {
      const c0 = M[k][0], ci = M[k][i];
      M[k][0] = yg * c0 - xg * ci;
      M[k][i] = alpha * c0 + beta * ci;
    }
  }
  if (M[0][0] < 0n) {
    for (let k = 0; k < d; k++) M[k][0] = -M[k][0];
  }
  return { G: M, steps };
}

// Wait: What if we also track the effect of the column operations of triangularizeBlock on an external top block?
function triangularizeWithTop(c, vHat, H, topCol0) {
  // topCol0 is a vector of length rowsTop (the entries in column 0 of the top block; cols 1..r of top block are 0!)
  const r = H.length;
  const d = r + 1;
  const Hhat = liftTop(H);
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
  let steps = 0;
  for (let i = r; i >= 1; i--) {
    const x = M[i][0];
    if (x === 0n) continue;
    steps++;
    const y = M[i][i];
    const [g, alpha, beta] = extGcd(x, y);
    const yg = y / g, xg = x / g;
    for (let k = 0; k <= i; k++) {
      const c0 = M[k][0], ci = M[k][i];
      M[k][0] = yg * c0 - xg * ci;
      M[k][i] = alpha * c0 + beta * ci;
    }
    for (let k = 0; k < Top.length; k++) {
      const c0 = Top[k][0], ci = Top[k][i];
      Top[k][0] = yg * c0 - xg * ci;
      Top[k][i] = alpha * c0 + beta * ci;
    }
  }
  if (M[0][0] < 0n) {
    for (let k = 0; k < d; k++) M[k][0] = -M[k][0];
    for (let k = 0; k < Top.length; k++) Top[k][0] = -Top[k][0];
  }
  return { G: M, Top, steps };
}

// Build Htilde_s recursively WITHOUT ANY cross-block Euclidean reduction!
function buildHtilde(b, s) {
  const H1 = buildH1(b);
  if (s === 1) return { H: H1, totalSteps: 0 };
  const { H: Hprev, totalSteps: prevSteps } = buildHtilde(b, s - 1);
  const rPrev = Hprev.length;
  const dPrev = rPrev + 1;
  const r = b ** s - 1;
  const vHat = buildVHat(b, s - 1);
  const vHatBot = vHat.slice(1);
  const Mb = H1[0][0];

  // Block 1: c = Mb, topCol0 in Block 0 is -Mb * vHatBot
  const top1 = vHatBot.map(v => -Mb * v);
  const { G: G1, Top: U1, steps: steps1 } = triangularizeWithTop(Mb, vHat, Hprev, top1);

  // Block m (m = 2..b-1): c = 1, top in Block 0 is -(u_{m-1}+1)*vHatBot, in Block 1 is +u_{m-1}*vHat
  // Notice we only need to run triangularizeWithTop ONCE with topCol0 = [vHatBot, vHat] and then scale!
  const { G: Gm, Top: TopM, steps: stepsM } = triangularizeWithTop(1n, vHat, Hprev, [...vHatBot, ...vHat]);
  const UmBase = TopM.slice(0, rPrev);
  const VmBase = TopM.slice(rPrev);

  const H = Array.from({ length: r }, () => Array(r).fill(0n));
  // Fill Block 0 diagonal: 2 * Hprev
  for (let i = 0; i < rPrev; i++) {
    for (let j = i; j < rPrev; j++) H[i][j] = 2n * Hprev[i][j];
  }
  // Fill Block 1 (cols rPrev .. rPrev + dPrev - 1)
  for (let i = 0; i < rPrev; i++) {
    for (let j = 0; j < dPrev; j++) H[i][rPrev + j] = U1[i][j];
  }
  for (let i = 0; i < dPrev; i++) {
    for (let j = i; j < dPrev; j++) H[rPrev + i][rPrev + j] = G1[i][j];
  }
  // Fill Blocks m = 2 .. b-1
  for (let m = 2; m < b; m++) {
    const colOff = rPrev + (m - 1) * dPrev;
    const rowOff = rPrev + (m - 1) * dPrev;
    const u = H1[0][m - 1];
    for (let i = 0; i < rPrev; i++) {
      for (let j = 0; j < dPrev; j++) H[i][colOff + j] = -(u + 1n) * UmBase[i][j];
    }
    for (let i = 0; i < dPrev; i++) {
      for (let j = 0; j < dPrev; j++) H[rPrev + i][colOff + j] = u * VmBase[i][j];
    }
    for (let i = 0; i < dPrev; i++) {
      for (let j = i; j < dPrev; j++) H[rowOff + i][colOff + j] = Gm[i][j];
    }
  }
  return { H, totalSteps: prevSteps + steps1 + stepsM, steps1, stepsM };
}

function invertUpperTriangular(H) {
  const r = H.length;
  const inv = Array.from({ length: r }, () => Array.from({ length: r }, () => new Fr(0n)));
  for (let i = r - 1; i >= 0; i--) {
    inv[i][i] = new Fr(1n, H[i][i]);
    for (let j = i + 1; j < r; j++) {
      let s = new Fr(0n);
      for (let k = i + 1; k <= j; k++) {
        if (H[i][k] !== 0n) s = s.add(new Fr(H[i][k]).mul(inv[k][j]));
      }
      inv[i][j] = s.mul(new Fr(-1n, H[i][i]));
    }
  }
  return inv;
}

for (const [b, s] of [[3, 1], [3, 2], [3, 3], [3, 4], [5, 1], [5, 2], [5, 3], [7, 1], [7, 2], [9, 2]]) {
  const { H, steps1, stepsM } = buildHtilde(b, s);
  const r = H.length;
  const D = 2n ** BigInt(s);
  const Hinv = invertUpperTriangular(H);
  const l1 = Hinv.map(row => row.reduce((acc, x) => acc.add(x.abs()), new Fr(0n)).mul(new Fr(D)));
  let maxRowL1 = new Fr(0n);
  for (const val of l1) if (val.sub(maxRowL1).p > 0n) maxRowL1 = val;
  const lastRow = Array.from({ length: r }, (_, j) => {
    let sum = new Fr(0n);
    for (let i = 0; i < r - 1; i++) sum = sum.sub(Hinv[i][j]);
    sum = sum.add(Hinv[r - 1][j]);
    return sum.mul(new Fr(D));
  });
  const lastL1 = lastRow.reduce((acc, x) => acc.add(x.abs()), new Fr(0n));
  const K = maxRowL1.sub(lastL1).p > 0n ? maxRowL1 : lastL1;
  console.log(`b=${b}, s=${s}, d=${b**s}: steps1=${steps1||0}, stepsM=${stepsM||0}, K(Htilde) = ${K.toString()} (${K.toNumber().toFixed(4)}), K/d = ${(K.toNumber()/(b**s)).toFixed(4)}, maxRowL1(D*Htilde^-1) = ${maxRowL1.toString()} (${maxRowL1.toNumber().toFixed(4)})`);
}
