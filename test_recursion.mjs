import fs from "fs";

const src = fs.readFileSync("analyze_lattice.mjs", "utf8").split("for (const [b, s]")[0];
eval(src.replace(/class Fr/g, "globalThis.Fr = class Fr")
        .replace(/function gcd/g, "globalThis.gcd = function gcd")
        .replace(/function lambda1Basis/g, "globalThis.lambda1Basis = function lambda1Basis")
        .replace(/function lambdaSBasis/g, "globalThis.lambdaSBasis = function lambdaSBasis")
        .replace(/function BMatrix/g, "globalThis.BMatrix = function BMatrix")
        .replace(/function floorDiv/g, "globalThis.floorDiv = function floorDiv")
        .replace(/function hnfUpper/g, "globalThis.hnfUpper = function hnfUpper")
        .replace(/function invertUpperTriangular/g, "globalThis.invertUpperTriangular = function invertUpperTriangular"));

// Inspect how H_s is built recursively from H_{s-1} and H_1!
// Specifically, look at H_s^-1 (or F_s^-1 = D * H_s^-1):
// Why not look at the columns/rows of F_s = H_s / D directly?
for (const [b, s] of [[3, 1], [3, 2], [3, 3], [5, 1], [5, 2]]) {
  const { basis, v } = lambdaSBasis(b, s, new Fr(1n, 2n));
  const B = BMatrix(basis);
  const r = B.length;
  const D = 2n ** BigInt(s);
  const A = B.map(row => row.map(x => x.mul(new Fr(D)).p));
  const Arev = Array.from({ length: r }, (_, i) => A[r - 1 - i]);
  const H = hnfUpper(Arev);
  const Hinv = invertUpperTriangular(H);

  // Print max entry of D * Hinv and max row L1 norm of D * Hinv
  let maxEntry = new Fr(0n);
  for (let i = 0; i < r; i++) {
    for (let j = 0; j < r; j++) {
      const val = Hinv[i][j].abs().mul(new Fr(D));
      if (val.sub(maxEntry).p > 0n) maxEntry = val;
    }
  }
  console.log(`b=${b}, s=${s}: max_{i,j} |D * Hinv[i,j]| = ${maxEntry.toString()} (${maxEntry.toNumber().toFixed(4)})`);
}
