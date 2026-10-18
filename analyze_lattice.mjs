// Exact BigInt Rational arithmetic & HNF analysis for Lambda_s (Erdos #1 Buffer Theorem)

function gcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    let t = b;
    b = a % b;
    a = t;
  }
  return a;
}

class Fr {
  constructor(p, q = 1n) {
    p = BigInt(p);
    q = BigInt(q);
    if (q < 0n) { p = -p; q = -q; }
    const g = gcd(p, q);
    this.p = p / g;
    this.q = q / g;
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

function lambda1Basis(b, alpha = new Fr(1n, 2n)) {
  const T = (z) => z.map((_, i) => z[i].add(alpha.mul(z[(i + 1) % b])));
  const basis = [];
  for (let i = 0; i < b - 1; i++) {
    const z = Array.from({ length: b }, () => new Fr(0n));
    z[i] = new Fr(1n);
    z[i + 1] = new Fr(-1n);
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
      for (let j = 0; j < b; j++) {
        for (const x of v) w.push(beta[j].mul(x));
      }
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
  const B = Array.from({ length: r }, () => Array(r));
  for (let i = 0; i < r; i++) {
    for (let j = 0; j < r; j++) {
      B[i][j] = basis[j][i];
    }
  }
  return B;
}

// Euclidean floor division for BigInt (a = q*b + r, 0 <= r < |b|)
function floorDiv(a, b) {
  let q = a / b;
  let r = a % b;
  if (r !== 0n && ((r < 0n) !== (b < 0n))) {
    q -= 1n;
  }
  return q;
}

function hnfUpper(A) {
  const r = A.length;
  const H = A.map(row => row.slice());
  // Eliminate bottom-up: for row i = r-1 down to 0, make H[i][0..i-1] = 0 and H[i][i] > 0
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
        for (let k = 0; k < r; k++) {
          const tmp = H[k][i]; H[k][i] = H[k][jPiv]; H[k][jPiv] = tmp;
        }
      }
      for (let j = 0; j < i; j++) {
        if (H[i][j] !== 0n) {
          const q = floorDiv(H[i][j], H[i][i]);
          for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
        }
      }
    }
    // Ensure pivot is at column i
    if (H[i][i] === 0n) {
      for (let j = 0; j < i; j++) {
        if (H[i][j] !== 0n) {
          for (let k = 0; k <= i; k++) {
            const tmp = H[k][i]; H[k][i] = H[k][j]; H[k][j] = tmp;
          }
          break;
        }
      }
    }
    if (H[i][i] < 0n) {
      for (let k = 0; k <= i; k++) H[k][i] = -H[k][i];
    }
    // Reduce columns j > i in row i so that 0 <= H[i][j] < H[i][i]
    for (let j = i + 1; j < r; j++) {
      const q = floorDiv(H[i][j], H[i][i]);
      if (q !== 0n) {
        for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
      }
    }
  }
  return H;
}

function invertUpperTriangular(H) {
  const r = H.length;
  const inv = Array.from({ length: r }, () => Array.from({ length: r }, () => new Fr(0n)));
  for (let i = r - 1; i >= 0; i--) {
    inv[i][i] = new Fr(1n, H[i][i]);
    for (let j = i + 1; j < r; j++) {
      let s = new Fr(0n);
      for (let k = i + 1; k <= j; k++) {
        if (H[i][k] !== 0n) {
          s = s.add(new Fr(H[i][k]).mul(inv[k][j]));
        }
      }
      inv[i][j] = s.mul(new Fr(-1n, H[i][i]));
    }
  }
  return inv;
}

for (const [b, s] of [[3, 1], [3, 2], [3, 3], [5, 1], [5, 2], [7, 1]]) {
  const { basis } = lambdaSBasis(b, s, new Fr(1n, 2n));
  const B = BMatrix(basis);
  const r = B.length;
  const D = 2n ** BigInt(s);
  const A = B.map(row => row.map(x => x.mul(new Fr(D)).p));
  const Arev = Array.from({ length: r }, (_, i) => A[r - 1 - i]);
  const H = hnfUpper(Arev);
  const Hinv = invertUpperTriangular(H);

  const l1 = Hinv.map(row => row.reduce((acc, x) => acc.add(x.abs()), new Fr(0n)));
  const lastRow = Array.from({ length: r }, (_, j) => {
    let sum = new Fr(0n);
    for (let i = 0; i < r - 1; i++) sum = sum.sub(Hinv[i][j]);
    sum = sum.add(Hinv[r - 1][j]);
    return sum;
  });
  const lastL1 = lastRow.reduce((acc, x) => acc.add(x.abs()), new Fr(0n));
  let maxL1 = lastL1;
  let maxRowIdx = r;
  for (let i = 0; i < r - 1; i++) {
    if (l1[i].sub(maxL1).p > 0n) {
      maxL1 = l1[i];
      maxRowIdx = i;
    }
  }
  const K = new Fr(D).mul(maxL1);
  const sumAllRowsL1 = l1.reduce((acc, x) => acc.add(x), new Fr(0n)).mul(new Fr(D));
  console.log(`b=${b}, s=${s}, d=${b**s}, r=${r}: K=${K.toNumber().toFixed(4)} (${K.toString()}), K/d=${(K.toNumber()/(b**s)).toFixed(4)}, D*sum_i ||row_i(H^-1)||_1 = ${sumAllRowsL1.toNumber().toFixed(4)} (${sumAllRowsL1.toString()})`);
  console.log(`  diag(H) = [${H.map((row, i) => Number(row[i])).join(", ")}]`);
  if (r <= 26) {
    H.forEach((row, i) => {
      if (row.filter(v => v !== 0n).length > 1 || s === 1) {
        const nz = [];
        row.forEach((v, j) => { if (v !== 0n) nz.push(`${j}:${v}`); });
        const invNz = [];
        Hinv[i].forEach((v, j) => { if (!v.isZero()) invNz.push(`${j}:${v.mul(new Fr(D)).toString()}`); });
        console.log(`    row ${String(i).padStart(2)}: H=[${nz.join(", ")}] | D*Hinv=[${invNz.join(", ")}]`);
      }
    });
  }
}

