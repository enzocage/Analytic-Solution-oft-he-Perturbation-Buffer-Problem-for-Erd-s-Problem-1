// Test buildHtilde WITH the 5-line off-diagonal reduction (which produces the EXACT canonical HNF H_s!)
// and check how many active rows triangularizeWithTop touches!

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
function reduceOffDiag(H) {
  const r = H.length;
  for (let i = r - 1; i >= 0; i--) {
    for (let j = i + 1; j < r; j++) {
      const q = floorDiv(H[i][j], H[i][i]);
      if (q !== 0n) {
        for (let k = 0; k <= i; k++) H[k][j] -= q * H[k][i];
      }
    }
  }
  return H;
}

function buildH1(b) {
  const r = b - 1;
  const H1 = Array.from({ length: r }, () => Array(r).fill(0n));
  const Mb = (2n ** BigInt(b) + 1n) / 3n;
  H1[0][0] = Mb;
  for (let j = 1; j < r; j++) {
    H1[0][j] = (j % 2 === 1)
      ? (2n ** BigInt(j + 1) - 1n) / 3n
      : (2n ** BigInt(b) - 2n ** BigInt(j + 1)) / 3n;
    H1[j][j] = 1n;
  }
  return H1;
}

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

function triangularizeWithTop(c, vHat, H, topCol0, label = "") {
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
  const activeRows = [];
  const nonTrivialRows = [];
  for (let i = r; i >= 1; i--) {
    const x = M[i][0];
    if (x === 0n) continue;
    const y = M[i][i];
    const [g, alpha, beta] = extGcd(x, y);
    activeRows.push(`i=${i}(x=${x},y=${y},g=${g})`);
    if (g < y) nonTrivialRows.push(`i=${i}(x=${x},y=${y},g=${g},y/g=${y/g})`);
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
  if (label) {
    console.log(`  [${label}] active (${activeRows.length}), nonTrivial g<y (${nonTrivialRows.length}): [${nonTrivialRows.join(", ")}]`);
  }
  return { G: M, Top, steps: activeRows.length };
}

function buildFastHNF(b, s) {
  const H1 = buildH1(b);
  if (s === 1) return H1;
  const Hprev = buildFastHNF(b, s - 1);
  const rPrev = Hprev.length, dPrev = rPrev + 1, r = b ** s - 1;
  const vHat = buildVHat(b, s - 1);
  const vHatBot = vHat.slice(1);
  const Mb = H1[0][0];
  console.log(`\nFastHNF(b=${b}, s=${s}, r=${r}):`);
  const { G: G1, Top: U1 } = triangularizeWithTop(Mb, vHat, Hprev, vHatBot.map(v => -Mb * v), `b=${b},s=${s},Block1`);
  const { G: Gm, Top: TopM } = triangularizeWithTop(1n, vHat, Hprev, [...vHatBot, ...vHat], `b=${b},s=${s},BlockM`);
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

for (const [b, s] of [[3, 4], [5, 3], [7, 3]]) {
  buildFastHNF(b, s);
}
