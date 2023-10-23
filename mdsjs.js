/**
 * Created by krause on 2014-10-25.
 */

class Continuation extends Error {
  constructor(f) {
    super();
    this.__continuation = f;
  }
}

export class Matrix {
  constructor(
    /** @type {MDSJs} */ mdsjs,
    /** @type {Float32Array | Float64Array} */ mat,
    /** @type {number} */ rows,
    /** @type {number} */ cols,
  ) {
    /** @type {MDSJs} */ this._mdsjs = mdsjs;
    /** @type {Float32Array | Float64Array} */ this._mat = mat;
    /** @type {number} */ this._rows = rows;
    /** @type {number} */ this._cols = cols;
  }

  mdsjs() {
    return this._mdsjs;
  }

  rows() {
    return this._rows;
  }

  cols() {
    return this._cols;
  }

  isQuadratic() {
    return this._rows === this._cols;
  }

  noNaNs() {
    this._mdsjs.noNaNs(this._mat);
  }

  someRows(
    /** @type {(row: Float32Array | Float64Array, ix: number) => boolean} */ cb,
  ) {
    let pos = 0;
    for (let r = 0; r < this._rows; r += 1) {
      if (cb(this._mat.subarray(pos, pos + this._cols), r)) {
        return true;
      }
      pos += this._cols;
    }
    return false;
  }

  everyRows(
    /** @type {(row: Float32Array | Float64Array, ix: number) => boolean} */ cb,
  ) {
    return !this.someRows((row, ix) => {
      return !cb(row, ix);
    });
  }

  rowsIter(
    /** @type {(row: Float32Array | Float64Array, ix: number) => void} */ cb,
  ) {
    let pos = 0;
    for (let r = 0; r < this._rows; r += 1) {
      cb(this._mat.subarray(pos, pos + this._cols), r);
      pos += this._cols;
    }
  }

  rowIter(
    /** @type {number} */ row,
    /** @type {(val: number, row: number, col: number) => void} */ cb,
  ) {
    let pos = row * this._cols;
    for (let ix = 0; ix < this._cols; ix += 1) {
      cb(this._mat[pos], row, ix);
      pos += 1;
    }
  }

  colIter(
    /** @type {number} */ col,
    /** @type {(val: number, row: number, col: number) => void} */ cb,
  ) {
    let pos = col;
    for (let ix = 0; ix < this._rows; ix += 1) {
      cb(this._mat[pos], ix, col);
      pos += this._cols;
    }
  }

  getUnsafe(/** @type {number} */ pos) {
    return this._mat[pos];
  }

  createArray(/** @type {number} */ rows, /** @type {number} */ cols) {
    const size = rows * cols;
    return this._mat.byteLength > 24
      ? new Float64Array(size)
      : new Float32Array(size);
  }

  toString() {
    let res = '';
    for (let r = 0; r < this.rows(); r += 1) {
      this.rowIter(r, (e, c) => {
        res += ` ${e}`;
      });
      res += '\n';
    }
    return res;
  }

  iter(
    /** @type {Matrix} */ matB,
    /** @type {number} */ row,
    /** @type {number} */ col,
    /** @type {(left: number, right: number, row: number, pos: number, col: number) => void} */ cb,
  ) {
    Matrix.iter(this, matB, row, col, cb);
  }

  mul(/** @type {Matrix} */ matB) {
    return Matrix.mul(this, matB);
  }

  add(/** @type {Matrix} */ matB) {
    return Matrix.add(this, matB);
  }

  neg() {
    const mat = this.createArray(this.rows(), this.cols());
    for (let pos = 0; pos < mat.length; pos += 1) {
      mat[pos] = -this.getUnsafe(pos);
    }
    return new Matrix(this._mdsjs, mat, this.rows(), this.cols());
  }

  scale(/** @type {number} */ scale) {
    const mat = this.createArray(this.rows(), this.cols());
    for (let pos = 0; pos < mat.length; pos += 1) {
      mat[pos] = scale * this.getUnsafe(pos);
    }
    return new Matrix(this._mdsjs, mat, this.rows(), this.cols());
  }

  squareElements() {
    const mat = this.createArray(this.rows(), this.cols());
    for (let pos = 0; pos < mat.length; pos += 1) {
      mat[pos] = this.getUnsafe(pos) * this.getUnsafe(pos);
    }
    return new Matrix(this._mdsjs, mat, this.rows(), this.cols());
  }

  colCenter() {
    const rows = this.rows();
    const cols = this.cols();
    const mat = this.createArray(rows, cols);
    for (let c = 0; c < cols; c += 1) {
      let avg = 0;
      this.colIter(c, (v) => {
        avg += v;
      });
      avg /= rows;
      let pos = c;
      this.colIter(c, (v) => {
        mat[pos] = v - avg;
        pos += cols;
      });
    }
    return new Matrix(this._mdsjs, mat, rows, cols);
  }

  rowCenter() {
    const rows = this.rows();
    const cols = this.cols();
    const mat = this.createArray(rows, cols);
    let pos = 0;
    for (let r = 0; r < rows; r += 1) {
      let avg = 0;
      this.rowIter(r, (v) => {
        avg += v;
      });
      avg /= cols;
      this.rowIter(r, (v) => {
        mat[pos] = v - avg;
        pos += 1;
      });
    }
    return new Matrix(this._mdsjs, mat, rows, cols);
  }

  doubleCenter() {
    const rows = this.rows();
    const cols = this.cols();
    const mat = this.createArray(rows, cols);
    for (let r = 0; r < rows; r += 1) {
      let avg = 0;
      this.rowIter(r, (v) => {
        avg += v;
      });
      avg /= cols;
      let pos = r * cols;
      this.rowIter(r, (v) => {
        mat[pos] = v - avg;
        pos += 1;
      });
    }
    for (let c = 0; c < cols; c += 1) {
      let avg = 0;
      let pos = c;
      for (let r = 0; r < rows; r += 1) {
        avg += mat[pos];
        pos += cols;
      }
      avg /= rows;
      pos = c;
      for (let r = 0; r < rows; r += 1) {
        mat[pos] -= avg;
        pos += cols;
      }
    }
    return new Matrix(this._mdsjs, mat, rows, cols);
  }

  distance(/** @type {number} */ colA, /** @type {number} */ colB) {
    let res = 0;
    let posA = colA;
    let posB = colB;
    for (let r = 0; r < this.rows(); r += 1) {
      const v = this.getUnsafe(posA) - this.getUnsafe(posB);
      res += v * v;
      posA += this.cols();
      posB += this.cols();
    }
    return Math.sqrt(res);
  }

  eigen(/** @type {Float32Array | Float64Array} */ eigenVals) {
    let /** @type {Matrix | undefined} */ res;
    this.eigenAsync(
      eigenVals,
      (mat) => {
        res = mat;
      },
      this._mdsjs.getCallDirect(),
    );
    return res;
  }

  eigenAsync(
    /** @type {Float32Array | Float64Array} */ eigenVals,
    /** @type {(mat: Matrix) => void} */ cb,
    /** @type {(cb: () => void) => void} */ argCall,
  ) {
    const call = argCall ?? this._mdsjs.CALL_ASYNC;
    const d = eigenVals.length;
    const rows = this.rows();
    const cols = this.cols();
    const content = this.createArray(rows, cols);
    let pos = 0;
    for (let r = 0; r < rows; r += 1) {
      this.rowIter(r, (v) => {
        content[pos] = v;
        pos += 1;
      });
    }
    const eigenVecs = this.createArray(d, rows);
    let ePos = -rows;
    let m = 0;
    let r = 0;
    let iter = 0;

    const innerLoop = () => {
      for (let ix = 0; ix < this._mdsjs.EIGEN_ITER_ASYNC; ix += 1) {
        if (
          !(
            Math.abs(1 - r) > this._mdsjs.EIGEN_EPS &&
            iter < this._mdsjs.EIGEN_ITER
          )
        ) {
          m += 1;
          iterate();
          return;
        }
        const q = this.createArray(1, rows);
        pos = 0;
        for (let rix = 0; rix < rows; rix += 1) {
          for (let cix = 0; cix < cols; cix += 1) {
            q[rix] += content[pos] * eigenVecs[ePos + cix];
            pos += 1;
          }
        }
        eigenVals[m] = this._mdsjs.prod(eigenVecs, ePos, q, 0, rows);
        this._mdsjs.normalizeVec(q);
        r = Math.abs(this._mdsjs.prod(eigenVecs, ePos, q, 0, rows));
        this._mdsjs.xcopy(q, 0, eigenVecs, ePos, rows);
        iter += 1;
      }
      call(innerLoop);
    }; // innerLoop

    const iterate = () => {
      if (!(m < d)) {
        cb(new Matrix(this._mdsjs, eigenVecs, d, rows));
        return;
      }

      if (m > 0) {
        pos = 0;
        for (let rix = 0; rix < rows; rix += 1) {
          for (let cix = 0; cix < cols; cix += 1) {
            content[pos] -=
              eigenVals[m - 1] * eigenVecs[ePos + rix] * eigenVecs[ePos + cix];
            pos += 1;
          }
        }
      }
      ePos += rows;
      pos = ePos;
      for (let ix = 0; ix < rows; ix += 1) {
        eigenVecs[pos] = Math.random();
        pos += 1;
      }
      this._mdsjs.normalizeVec(eigenVecs, ePos, ePos + rows);
      r = 0;
      iter = 0;

      call(innerLoop);
    }; // iterate

    iterate();
  }

  powerIter() {
    let /** @type {Float32Array | Float64Array | undefined} */ res;
    this.powerIterAsync((r) => {
      res = r;
    }, this._mdsjs.getCallDirect());
    return res;
  }

  powerIterAsync(
    /** @type {(r: Float32Array | Float64Array) => void} */ cb,
    /** @type {(cb: () => void) => void} */ argCall,
  ) {
    const call = argCall ?? this._mdsjs.CALL_ASYNC;
    const rows = this.rows();
    const cols = this.cols();
    let r = this.createArray(1, cols);
    for (let ix = 0; ix < cols; ix += 1) {
      r[ix] = Math.random();
    }
    let len = Number.POSITIVE_INFINITY;
    let stop = false;
    let iter = 0;

    const iterate = () => {
      for (let ix = 0; ix < this._mdsjs.EIGEN_ITER_ASYNC; ix += 1) {
        if (iter >= this._mdsjs.EIGEN_ITER || stop) {
          cb(r);
          return;
        }
        const s = this.createArray(1, cols);
        for (let row = 0; row < rows; row += 1) {
          let prod = 0;
          this.rowIter(row, (v, row, col) => {
            prod += v * r[col];
          });
          this.rowIter(row, (v, row, col) => {
            s[col] += prod * v;
          });
        }
        const nl = this._mdsjs.lengthSq(s);
        if (Math.abs(len - nl) < this._mdsjs.EIGEN_EPS) {
          stop = true;
        }
        len = nl;
        this._mdsjs.normalizeVec(s);
        r = s;
        iter += 1;
      }
      call(iterate);
    };

    iterate();
  }

  static iter(
    /** @type {Matrix} */ matA,
    /** @type {Matrix} */ matB,
    /** @type {number} */ row,
    /** @type {number} */ col,
    /** @type {(left: number, right: number, row: number, pos: number, col: number) => void} */ cb,
  ) {
    if (matA.cols() !== matB.rows()) {
      console.warn(
        'incompatible dimensions',
        matA.rows() + 'x' + matA.cols(),
        matB.rows() + 'x' + matB.cols(),
      );
      return;
    }
    let posA = row * matA.cols();
    let posB = col;
    for (let ix = 0; ix < matA.cols(); ix += 1) {
      cb(matA.getUnsafe(posA), matB.getUnsafe(posB), row, ix, col);
      posA += 1;
      posB += matB.cols();
    }
  }

  static mul(/** @type {Matrix} */ matA, /** @type {Matrix} */ matB) {
    if (matA.cols() !== matB.rows()) {
      console.warn(
        'incompatible dimensions',
        `${matA.rows()}x${matA.cols()}`,
        `${matB.rows()}x${matB.cols()}`,
      );
      return null;
    }
    // cache friendly iteration (a rows -> a cols/b rows -> b cols)
    // TODO experiment with (a cols -> a rows -> b cols)
    const mat = matA.createArray(matA.rows(), matB.cols());
    for (let r = 0; r < matA.rows(); r += 1) {
      matA.rowIter(r, (a, _, k) => {
        let pos = r * matB.cols();
        matB.rowIter(k, (b, _, __) => {
          mat[pos] += a * b;
          pos += 1;
        });
      });
    }
    return new Matrix(matA.mdsjs(), mat, matA.rows(), matB.cols());
  }

  static add(/** @type {Matrix} */ matA, /** @type {Matrix} */ matB) {
    if (matA.rows() !== matB.rows() || matA.cols() !== matB.cols()) {
      console.warn(
        'incompatible dimensions',
        `${matA.rows()}x${matA.cols()}`,
        `${matB.rows()}x${matB.cols()}`,
      );
      return null;
    }
    const mat = matA.createArray(matA.rows(), matA.cols());
    for (let pos = 0; pos < mat.length; pos += 1) {
      mat[pos] = matA.getUnsafe(pos) + matB.getUnsafe(pos);
    }
    return new Matrix(matA.mdsjs(), mat, matA.rows(), matA.cols());
  }
} // Matrix

export class MDSJs {
  /** @type {boolean} */ DEBUG = false;
  /** @type {number} */ GRAM_SCHMIDT_EPS = 1e-12;
  /** @type {number} */ EIGEN_EPS = 1e-7;
  /** @type {number} */ EIGEN_ITER = 10000;
  /** @type {number} */ EIGEN_ITER_ASYNC = 200;

  constructor() {}

  noNaNs(/** @type {Float32Array | Float64Array | number[]} */ arr) {
    for (let ix = 0; ix < arr.length; ix += 1) {
      if (Number.isNaN(arr[ix])) {
        throw new Error('NaN in array');
      }
    }
  }

  noZeros(/** @type {Float32Array | Float64Array | number[]} */ arr) {
    for (let ix = 0; ix < arr.length; ix += 1) {
      if (Number.isNaN(arr[ix]) || arr[ix] === 0) {
        throw new Error(arr[ix] === 0 ? '0 in array' : 'NaN in array');
      }
    }
  }

  onlyPositive(/** @type {Float32Array | Float64Array | number[]} */ arr) {
    for (let ix = 0; ix < arr.length; ix += 1) {
      if (Number.isNaN(arr[ix]) || !(arr[ix] > 0)) {
        throw new Error(
          !(arr[ix] > 0) ? `${arr[ix]} in array` : 'NaN in array',
        );
      }
    }
  }

  CALL_ASYNC(/** @type {() => void} */ cb) {
    setTimeout(cb, 0);
  }

  getCallDirect() {
    let depth = 0;
    return (cb) => {
      if (depth > 20) {
        // prevent stack-overflows
        throw new Continuation(cb);
      }
      if (!depth) {
        let cc = cb;
        while (cc) {
          try {
            depth += 1;
            cc();
            cc = null;
          } catch (e) {
            if (!e.__continuation) {
              throw e;
            }
            cc = e.__continuation;
          }
          depth = 0;
        }
      } else {
        depth += 1;
        cb();
      }
    };
  }

  pca(/** @type {Matrix} */ positions) {
    let res;
    this.pcaAsync(
      positions,
      (mat) => {
        res = mat;
      },
      this.getCallDirect(),
    );
    return res;
  }

  pcaAsync(
    /** @type {Matrix} */ positions,
    /** @type {(mat: Matrix) => void} */ cb,
    /** @type {(cb: () => void) => void} */ argCall,
  ) {
    const call = arguments.length > 2 ? argCall : this.CALL_ASYNC;
    const centered = positions.colCenter();
    const cols = centered.cols();
    centered.powerIterAsync((pca0) => {
      const mat = this.removeComponent(centered, pca0);
      mat.powerIterAsync((pca1) => {
        const res = centered.createArray(cols, 2);
        for (let ix = 0; ix < cols; ix += 1) {
          res[2 * ix + 0] = pca0[ix];
          res[2 * ix + 1] = pca1[ix];
        }
        cb(new Matrix(this, res, cols, 2));
      }, call);
    }, call);
  }

  removeComponent(
    /** @type {Matrix} */ mat,
    /** @type {Float32Array | Float64Array} */ comp,
  ) {
    // Gram–Schmidt process
    const rows = mat.rows();
    const cols = mat.cols();
    if (comp.length !== cols) {
      console.warn('incompatible size', comp.length, cols);
    }

    const proj = (vec, from, sub, fromSub, len) => {
      const res = mat.createArray(1, len);
      let uv = 0;
      let uu = 0;
      for (let ix = 0; ix < len; ix += 1) {
        uv += sub[fromSub + ix] * vec[from + ix];
        uu += sub[fromSub + ix] * sub[fromSub + ix];
      }
      if (
        Math.abs(uv) < this.GRAM_SCHMIDT_EPS ||
        Math.abs(uu) < this.GRAM_SCHMIDT_EPS ||
        !Number.isFinite(uu) ||
        !Number.isFinite(uv)
      ) {
        for (let ix = 0; ix < len; ix += 1) {
          res[ix] = 0;
        }
      } else {
        for (let ix = 0; ix < len; ix += 1) {
          res[ix] = (uv / uu) * sub[fromSub + ix];
        }
      }
      return res;
    };

    const nextMat = mat.createArray(rows, cols);
    let pos = 0;
    for (let r = 0; r < rows; r += 1) {
      mat.rowIter(r, (v) => {
        nextMat[pos] = v;
        pos += 1;
      });
    }
    for (let r = 0; r < rows; r += 1) {
      const curPos = r * cols;
      this.normalizeVec(nextMat, curPos, curPos + cols);
      for (let ix = r + 1; ix < rows; ix += 1) {
        const mPos = ix * cols;
        const p = proj(nextMat, mPos, nextMat, curPos, cols);
        for (let c = 0; c < cols; c += 1) {
          nextMat[mPos + c] -= p[c];
        }
      }
    }
    return new Matrix(this, nextMat, rows, cols);
  }

  pcaPositions(/** @type {Matrix} */ positions) {
    const pca = this.pca(positions);
    return positions.mul(pca);
  }

  landmarkMDS(/** @type {Matrix} */ dist, /** @type {number} */ dims) {
    let res;
    this.landmarkMDSAsync(
      dist,
      dims,
      (mat) => {
        res = mat;
      },
      this.getCallDirect(),
    );
    return res;
  }

  landmarkMDSAsync(
    /** @type {Matrix} */ dist,
    /** @type {number} */ dims,
    /** @type {(mat: Matrix) => void} */ cb,
    /** @type {(cb: () => void) => void} */ argCall,
  ) {
    const landmarkMatrix = (/** @type {Matrix} */ mat) => {
      if (this.DEBUG) {
        mat.noNaNs();
      }
      const rows = mat.rows();
      const cols = mat.cols();
      const perm = new Uint32Array(rows);
      for (let r = 0; r < rows; r += 1) {
        mat.rowIter(r, (v, r, c) => {
          if (!v) {
            perm[r] = c;
          }
        });
      }
      if (this.DEBUG) {
        this.noNaNs(perm);
      }
      const lm = mat.createArray(rows, rows);
      let pos = 0;
      for (let r = 0; r < rows; r += 1) {
        const mPos = r * cols;
        for (let c = 0; c < cols; c += 1) {
          lm[pos] = mat.getUnsafe(mPos + perm[c]);
          pos += 1;
        }
      }
      if (this.DEBUG) {
        this.noNaNs(lm);
      }
      return new Matrix(this, lm, rows, rows);
    };

    const landmarkResult = (
      /** @type {Matrix} */ dist,
      /** @type {number} */ dims,
      /** @type {Matrix} */ eigenVecs,
      /** @type {Float32Array | Float64Array} */ eigenVals,
    ) => {
      const rows = dist.rows();
      const cols = dist.cols();
      const distSq = dist.squareElements();
      if (this.DEBUG) {
        distSq.noNaNs();
      }

      const mean = dist.createArray(1, cols);
      for (let c = 0; c < cols; c += 1) {
        distSq.colIter(c, (v) => {
          mean[c] += v;
        });
        mean[c] /= rows;
      }
      if (this.DEBUG) {
        this.noNaNs(mean);
        eigenVecs.noNaNs();
        this.noZeros(eigenVals);
      }
      const tmp = eigenVecs.createArray(eigenVecs.rows(), eigenVecs.cols());
      let pos = 0;
      for (let r = 0; r < eigenVecs.rows(); r += 1) {
        const div = Math.sqrt(Math.abs(eigenVals[r])); // TODO not sure how to handle negative values
        eigenVecs.rowIter(r, (v) => {
          tmp[pos] = v / div;
          pos += 1;
        });
      }
      if (this.DEBUG) {
        this.noNaNs(tmp);
      }

      const positions = dist.createArray(cols, dims);
      pos = 0;
      for (let e = 0; e < cols; e += 1) {
        const m = mean[e];
        let tPos = 0;
        for (let d = 0; d < dims; d += 1) {
          let cur = 0;
          distSq.colIter(e, (v) => {
            cur -= 0.5 * (v - m) * tmp[tPos];
            tPos += 1;
          });
          positions[pos] = cur;
          pos += 1;
        }
      }
      if (this.DEBUG) {
        this.noNaNs(positions);
      }
      return new Matrix(this, positions, cols, dims);
    };

    const call = argCall ?? this.CALL_ASYNC;
    const lm = landmarkMatrix(dist);
    const eigenVals = dist.createArray(1, dims);
    lm.squareElements()
      .doubleCenter()
      .scale(-0.5)
      .eigenAsync(
        eigenVals,
        (eigenVecs) => {
          cb(landmarkResult(dist, dims, eigenVecs, eigenVals));
        },
        call,
      );
  }

  normalizeVec(
    /** @type {Float32Array | Float64Array | number[]} */ vec,
    /** @type {number} */ f,
    /** @type {number} */ t,
  ) {
    const from = f ?? 0;
    const to = t ?? vec.length;
    let sum = 0;
    for (let ix = from; ix < to; ix += 1) {
      sum += vec[ix] * vec[ix];
    }
    sum = Math.sqrt(sum);
    if (sum < 1e-30 || !Number.isFinite(sum)) {
      // don't scale when really small
      return;
    }
    for (let ix = from; ix < to; ix += 1) {
      vec[ix] /= sum;
    }
  }

  lengthSq(
    /** @type {Float32Array | Float64Array | number[]} */ vec,
    /** @type {number} */ f,
    /** @type {number} */ t,
  ) {
    const from = f ?? 0;
    const to = t ?? vec.length;
    let sum = 0;
    for (let ix = from; ix < to; ix += 1) {
      sum += vec[ix] * vec[ix];
    }
    return sum;
  }

  prod(
    /** @type {Float32Array | Float64Array | number[]} */ vecA,
    /** @type {number} */ fromA,
    /** @type {Float32Array | Float64Array | number[]} */ vecB,
    /** @type {number} */ fromB,
    /** @type {number} */ len,
  ) {
    let sum = 0;
    let posA = fromA;
    let posB = fromB;
    for (let ix = 0; ix < len; ix += 1) {
      sum += vecA[posA] * vecB[posB];
      posA += 1;
      posB += 1;
    }
    return sum;
  }

  xcopy(
    /** @type {Float32Array | Float64Array | number[]} */ fromVec,
    /** @type {number} */ fromStart,
    /** @type {Float32Array | Float64Array | number[]} */ toVec,
    /** @type {number} */ toStart,
    /** @type {number} */ len,
  ) {
    let fromPos = fromStart;
    let toPos = toStart;
    for (let ix = 0; ix < len; ix += 1) {
      toVec[toPos] = fromVec[fromPos];
      fromPos += 1;
      toPos += 1;
    }
  }

  convertToMatrix(
    /** @type {number[][]} */ arrs,
    /** @type {boolean} */ useFloat32,
  ) {
    const rows = arrs.length;
    if (!rows) {
      console.warn('invalid dimension (rows)', rows);
      return null;
    }
    const cols = arrs[0].length;
    if (!cols) {
      console.warn('invalid dimension (cols)', cols);
      return null;
    }
    const size = rows * cols;
    const mat = useFloat32 ? new Float32Array(size) : new Float64Array(size);
    let pos = 0;
    for (let r = 0; r < rows; r += 1) {
      const row = arrs[r];
      if (row.length !== cols) {
        console.warn('invalid dimension in row ' + r, row.length, cols);
        return null;
      }
      for (let c = 0; c < cols; c += 1) {
        mat[pos] = row[c];
        pos += 1;
      }
    }
    return new Matrix(this, mat, rows, cols);
  }

  eye(
    /** @type {number} */ rows,
    /** @type {number} */ c,
    /** @type {boolean} */ useFloat32,
  ) {
    const cols = c ?? rows;
    const size = rows * cols;
    if (rows <= 0 || cols <= 0) {
      console.warn('invalid dimensions', rows, cols);
      return null;
    }
    const mat = useFloat32 ? new Float32Array(size) : new Float64Array(size);
    let pos = 0;
    for (let ix = 0; ix < Math.min(rows, cols); ix += 1) {
      mat[pos] = 1;
      pos += cols + 1;
    }
    return new Matrix(this, mat, rows, cols);
  }

  pivotRandom(/** @type {Matrix} */ m, /** @type {number} */ k) {
    if (!m.isQuadratic()) {
      console.warn('quadratic matrix needed', m.rows(), m.cols());
      return null;
    }
    if (k < m.rows()) {
      console.warn(
        'requested more pivots than elements',
        k,
        m.rows(),
        m.cols(),
      );
      return null;
    }
    const mat = m.createArray(k, m.cols());
    const pivots = {};
    let pos = 0;
    for (let ix = 0; ix < k; ix += 1) {
      let pivot = 0;
      do {
        pivot = Math.random() * m.cols();
      } while (pivots[pivot]);
      pivots[pivot] = true;
      for (let c = 0; c < m.cols(); c += 1) {
        mat[pos] = m.distance(pivot, c);
        pos += 1;
      }
    }
    return new Matrix(this, mat, k, m.cols());
  }
} // MDSJs

export default new MDSJs(); // create instance
