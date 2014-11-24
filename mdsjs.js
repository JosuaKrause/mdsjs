/**
 * Created by krause on 2014-10-25.
 */

mdsjs = function() {
  var thatMDS = this;

  this.pca = function(positions) {
    var centered = positions.colCenter();
    var rows = centered.rows();
    var cols = centered.cols();
    var pca0 = centered.powerIter();
    var mat = thatMDS.removeComponent(centered, pca0);
    var pca1 = mat.powerIter();
    var res = centered.createArray(cols, 2);
    for(var ix = 0;ix < cols;ix += 1) {
      res[2*ix + 0] = pca0[ix];
      res[2*ix + 1] = pca1[ix];
    }
    return new Matrix(res, cols, 2);
  };
  this.GRAM_SCHMIDT_EPS = 1e-12;
  this.removeComponent = function(mat, comp) {
    // Gram–Schmidt process
    var rows = mat.rows();
    var cols = mat.cols();
    comp.length === cols || console.warn("incompatible size", comp.length, cols);

    function proj(vec, from, sub, fromSub, len) {
      var res = mat.createArray(1, len);
      var uv = 0;
      var uu = 0;
      for(var ix = 0;ix < len;ix += 1) {
        uv += sub[fromSub + ix] * vec[from + ix];
        uu += sub[fromSub + ix] * sub[fromSub + ix];
      }
      if(Math.abs(uv) < thatMDS.GRAM_SCHMIDT_EPS || Math.abs(uu) < thatMDS.GRAM_SCHMIDT_EPS || !Number.isFinite(uu) || !Number.isFinite(uv)) {
        for(var ix = 0;ix < len;ix += 1) {
          res[ix] = 0;
        }
      } else {
        for(var ix = 0;ix < len;ix += 1) {
          res[ix] = uv / uu * sub[fromSub + ix];
        }
      }
      return res;
    }

    var nextMat = mat.createArray(rows, cols);
    var pos = 0;
    for(var r = 0;r < rows;r += 1) {
      mat.rowIter(r, function(v) {
        nextMat[pos] = v;
        pos += 1;
      });
    }
    for(var r = 0;r < rows;r += 1) {
      var curPos = r * cols;
      thatMDS.normalizeVec(nextMat, curPos, curPos + cols);
      for(var ix = r + 1;ix < rows;ix += 1) {
        var pos = ix * cols;
        var p = proj(nextMat, pos, nextMat, curPos, cols);
        for(var c = 0;c < cols;c += 1) {
          nextMat[pos + c] -= p[c];
        }
      }
    }
    return new Matrix(nextMat, rows, cols);
  };
  this.pcaPositions = function(positions) {
    var pca = thatMDS.pca(positions);
    return positions.mul(pca);
  };
  this.landmarkMDS = function(dist, dims) {
    var rows = dist.rows();
    var cols = dist.cols();
    var distSq = dist.squareElements();
    var mean = dist.createArray(1, cols);
    for(var c = 0;c < cols;c += 1) {
      distSq.colIter(c, function(v) {
        mean[c] += v;
      });
      mean[c] /= rows;
    }

    function landmarkMatrix(mat) {
      var rows = mat.rows();
      var cols = mat.cols();
      var perm = new Uint32Array(rows);
      for(var r = 0;r < rows;r += 1) {
        mat.rowIter(r, function(v, r, c) {
          if(!v) {
            perm[r] = c;
          }
        });
      }
      var lm = mat.createArray(rows, rows);
      var pos = 0;
      for(var r = 0;r < rows;r += 1) {
        var mPos = r * cols;
        for(var c = 0;c < cols;c += 1) {
          lm[pos] = mat.getUnsafe(mPos + perm[c]);
          pos += 1;
        }
      }
      return new Matrix(lm, rows, rows);
    }

    var landmarkMatrix = landmarkMatrix(dist);
    var eigenVals = dist.createArray(1, dims);
    var eigenVecs = landmarkMatrix.squareElements().doubleCenter().scale(-0.5).eigen(eigenVals);
    var tmp = eigenVecs.createArray(eigenVecs.rows(), eigenVecs.cols());
    var pos = 0;
    for(var r = 0;r < eigenVecs.rows();r += 1) {
      var div = Math.sqrt(eigenVals[r]);
      eigenVecs.rowIter(r, function(v) {
        tmp[pos] = v / div;
        pos += 1;
      });
    }
    var positions = dist.createArray(cols, dims);
    pos = 0;
    for(var e = 0;e < cols;e += 1) {
      var m = mean[e];
      var tPos = 0;
      for(var d = 0;d < dims;d += 1) {
        var cur = 0;
        distSq.colIter(e, function(v) {
          cur -= 0.5 * (v - m) * tmp[tPos];
          tPos += 1;
        });
        positions[pos] = cur;
        pos += 1;
      }
    }
    return new Matrix(positions, cols, dims);
  };
  this.normalizeVec = function(vec, f, t) {
    var from = arguments.length > 1 ? f : 0;
    var to = arguments.length > 2 ? t : vec.length;
    var sum = 0;
    for(var i = from;i < to;i += 1) {
      sum += vec[i] * vec[i];
    }
    sum = Math.sqrt(sum);
    if(sum < 1e-30 || !Number.isFinite(sum)) { // don't scale when really small
      return;
    }
    for(var i = from;i < to;i += 1) {
      vec[i] /= sum;
    }
  };
  this.lengthSq = function(vec, f, t) {
    var from = arguments.length > 1 ? f : 0;
    var to = arguments.length > 2 ? t : vec.length;
    var sum = 0;
    for(var i = from;i < to;i += 1) {
      sum += vec[i] * vec[i];
    }
    return sum;
  };
  this.prod = function(vecA, fromA, vecB, fromB, len) {
    var sum = 0;
    var posA = fromA;
    var posB = fromB;
    for(var i = 0;i < len;i += 1) {
      sum += vecA[posA] * vecB[posB];
      posA += 1;
      posB += 1;
    }
    return sum;
  };
  this.xcopy = function(fromVec, fromStart, toVec, toStart, len) {
    var fromPos = fromStart;
    var toPos = toStart;
    for(var i = 0;i < len;i += 1) {
      toVec[toPos] = fromVec[fromPos];
      fromPos += 1;
      toPos += 1;
    }
  };
  this.convertToMatrix = function(arrs, useFloat32) {
    var rows = arrs.length;
    if(!rows) {
      console.warn("invalid dimension (rows)", rows);
      return null;
    }
    var cols = arrs[0].length;
    if(!cols) {
      console.warn("invalid dimension (cols)", cols);
      return null;
    }
    var size = rows * cols;
    var mat = useFloat32 ? new Float32Array(size) : new Float64Array(size);
    var pos = 0;
    for(var r = 0;r < rows;r += 1) {
      var row = arrs[r];
      if(row.length !== cols) {
        console.warn("invalid dimension in row " + r, row.length, cols);
        return null;
      }
      for(var c = 0;c < cols;c += 1) {
        mat[pos] = row[c];
        pos += 1;
      }
    }
    return new Matrix(mat, rows, cols);
  };
  this.eye = function(rows, c, useFloat32) {
    var cols = arguments.length < 2 ? rows : c;
    var size = rows * cols;
    if(rows <= 0 || cols <= 0) {
      console.warn("invalid dimensions", rows, cols);
      return null;
    }
    var mat = useFloat32 ? new Float32Array(size) : new Float64Array(size);
    var pos = 0;
    for(var i = 0;i < Math.min(rows, cols);i += 1) {
      mat[pos] = 1;
      pos += cols + 1;
    }
    return new Matrix(mat, rows, cols);
  };
  this.pivotRandom = function(m, k) {
    if(!m.isQuadratic()) {
      console.warn("quadratic matrix needed", m.rows(), m.cols());
      return null;
    }
    if(k < m.rows()) {
      console.warn("requested more pivots than elements", k, m.rows(), m.cols());
      return null;
    }
    var mat = m.createArray(k, m.cols());
    var pivots = {};
    var pos = 0;
    for(var i = 0;i < k;i += 1) {
      var pivot = 0;
      do {
        pivot = Math.random() * m.cols();
      } while(pivots[pivot]);
      pivots[pivot] = true;
      for(var c = 0;c < m.cols();c += 1) {
        mat[pos] = m.distance(pivot, c);
        pos += 1;
      }
    }
    return new Matrix(mat, k, m.cols());
  };

  function Matrix(mat, rows, cols) {

    this.rows = function() {
      return rows;
    };
    this.cols = function() {
      return cols;
    };
    this.isQuadratic = function() {
      return rows === cols;
    };
    this.rowsIter = function(cb) {
      var pos = 0;
      for(var r = 0;r < rows;r += 1) {
        cb(mat.subarray(pos, pos + cols), r);
        pos += cols;
      }
    };
    this.rowIter = function(row, cb) {
      var pos = row * cols;
      for(var i = 0;i < cols;i += 1) {
        cb(mat[pos], row, i);
        pos += 1;
      }
    };
    this.colIter = function(col, cb) {
      var pos = col;
      for(var i = 0;i < rows;i += 1) {
        cb(mat[pos], i, col);
        pos += cols;
      }
    };
    this.getUnsafe = function(pos) {
      return mat[pos];
    };
    this.createArray = function(rows, cols) {
      var size = rows * cols;
      return mat.byteLength > 24 ? new Float64Array(size) : new Float32Array(size);
    };
  } // Matrix
  Matrix.prototype.toString = function() {
    var res = "";
    for(var r = 0;r < this.rows();r += 1) {
      this.rowIter(r, function(e, c) {
        res += " " + e;
      });
      res += "\n";
    }
    return res;
  };
  Matrix.prototype.iter = function(matB, row, col, cb) {
    Matrix.iter(this, matB, row, col, cb);
  };
  Matrix.prototype.mul = function(matB) {
    return Matrix.mul(this, matB);
  };
  Matrix.prototype.add = function(matB) {
    return Matrix.add(this, matB);
  };
  Matrix.prototype.neg = function() {
    var mat = this.createArray(this.rows(), this.cols());
    for(var pos = 0;pos < mat.length;pos += 1) {
      mat[pos] = -this.getUnsafe(pos);
    }
    return new Matrix(mat, this.rows(), this.cols());
  };
  Matrix.prototype.scale = function(scale) {
    var mat = this.createArray(this.rows(), this.cols());
    for(var pos = 0;pos < mat.length;pos += 1) {
      mat[pos] = scale * this.getUnsafe(pos);
    }
    return new Matrix(mat, this.rows(), this.cols());
  };
  Matrix.prototype.squareElements = function() {
    var mat = this.createArray(this.rows(), this.cols());
    for(var pos = 0;pos < mat.length;pos += 1) {
      mat[pos] = this.getUnsafe(pos) * this.getUnsafe(pos);
    }
    return new Matrix(mat, this.rows(), this.cols());
  };
  Matrix.prototype.colCenter = function() {
    var rows = this.rows();
    var cols = this.cols();
    var orig = this;
    var mat = this.createArray(rows, cols);
    for(var c = 0;c < cols;c += 1) {
      var avg = 0;
      orig.colIter(c, function(v) {
        avg += v;
      });
      avg /= rows;
      var pos = c;
      orig.colIter(c, function(v) {
        mat[pos] = v - avg;
        pos += cols;
      });
    }
    return new Matrix(mat, rows, cols);
  };
  Matrix.prototype.rowCenter = function() {
    var rows = this.rows();
    var cols = this.cols();
    var orig = this;
    var mat = this.createArray(rows, cols);
    var pos = 0;
    for(var r = 0;r < rows;r += 1) {
      var avg = 0;
      orig.rowIter(r, function(v) {
        avg += v;
      });
      avg /= cols;
      orig.rowIter(r, function(v) {
        mat[pos] = v - avg;
        pos += 1;
      });
    }
    return new Matrix(mat, rows, cols);
  };
  Matrix.prototype.doubleCenter = function() {
    var rows = this.rows();
    var cols = this.cols();
    var mat = this.createArray(rows, cols);
    for(var r = 0;r < rows;r += 1) {
      var avg = 0;
      this.rowIter(r, function(v) {
        avg += v;
      });
      avg /= cols;
      var pos = r * cols;
      this.rowIter(r, function(v) {
        mat[pos] = v - avg;
        pos += 1;
      });
    }
    for(var c = 0;c < cols;c += 1) {
      var avg = 0;
      var pos = c;
      for(var r = 0;r < rows;r += 1) {
        avg += mat[pos];
        pos += cols;
      }
      avg /= rows;
      pos = c;
      for(var r = 0;r < rows;r += 1) {
        mat[pos] -= avg;
        pos += cols;
      }
    }
    return new Matrix(mat, rows, cols);
  };
  Matrix.prototype.distance = function(colA, colB) {
    var res = 0;
    var posA = colA;
    var posB = colB;
    for(var r = 0;r < this.rows();r += 1) {
      var v = this.getUnsafe(posA) - this.getUnsafe(posB);
      res += v * v;
      posA += this.cols();
      posB += this.cols();
    }
    return Math.sqrt(res);
  };
  this.EIGEN_EPS = 1e-7;
  this.EIGEN_ITER = 10000;
  Matrix.prototype.eigen = function(eigenVals) {
    var mat = this;
    var d = eigenVals.length;
    var rows = mat.rows();
    var cols = mat.cols();
    var content = mat.createArray(rows, cols);
    var pos = 0;
    for(var r = 0;r < rows;r += 1) {
      mat.rowIter(r, function(v) {
        content[pos] = v;
        pos += 1;
      });
    }
    var eigenVecs = mat.createArray(d, rows);
    var ePos = -rows;
    for(var m = 0;m < d;m += 1) {
      if(m > 0) {
        pos = 0;
        for(var r = 0;r < rows;r += 1) {
          for(var c = 0;c < cols;c += 1) {
            content[pos] -= eigenVals[m - 1] * eigenVecs[ePos + r] * eigenVecs[ePos + c];
            pos += 1;
          }
        }
      }
      ePos += rows;
      pos = ePos;
      for(var i = 0;i < rows;i += 1) {
        eigenVecs[pos] = Math.random();
        pos += 1;
      }
      thatMDS.normalizeVec(eigenVecs, ePos, ePos + rows);
      var r = 0;
      for(var iter = 0;Math.abs(1 - r) > thatMDS.EIGEN_EPS && iter < thatMDS.EIGEN_ITER;iter += 1) {
        var q = mat.createArray(1, rows);
        pos = 0;
        for(var r = 0;r < rows;r += 1) {
          for(var c = 0;c < cols;c += 1) {
            q[r] += content[pos] * eigenVecs[ePos + c];
            pos += 1;
          }
        }
        eigenVals[m] = thatMDS.prod(eigenVecs, ePos, q, 0, rows);
        thatMDS.normalizeVec(q);
        r = Math.abs(thatMDS.prod(eigenVecs, ePos, q, 0, rows));
        thatMDS.xcopy(q, 0, eigenVecs, ePos, rows);
      }
    }
    return new Matrix(eigenVecs, d, rows);
  };
  Matrix.prototype.powerIter = function() {
    var mat = this;
    var rows = mat.rows();
    var cols = mat.cols();
    var r = mat.createArray(1, cols);
    for(var i = 0;i < cols;i += 1) {
      r[i] = Math.random();
    }
    var len = Number.POSITIVE_INFINITY;
    var stop = false;
    for(var iter = 0;iter < thatMDS.EIGEN_ITER && !stop;iter += 1) {
      var s = mat.createArray(1, cols);
      for(var row = 0;row < rows;row += 1) {
        var prod = 0;
        mat.rowIter(row, function(v, row, col) {
          prod += v * r[col];
        });
        mat.rowIter(row, function(v, row, col) {
          s[col] += prod * v;
        });
      }
      var nl = thatMDS.lengthSq(s);
      if(Math.abs(len - nl) < thatMDS.EIGEN_EPS) {
        stop = true;
      }
      len = nl;
      thatMDS.normalizeVec(s);
      r = s;
    }
    return r;
  };
  Matrix.iter = function(matA, matB, row, col, cb) {
    if(matA.cols() !== matB.rows()) {
      console.warn("incompatible dimensions", matA.rows() + "x" + matA.cols(), matB.rows() + "x" + matB.cols());
      return;
    }
    var posA = row * matA.cols();
    var posB = col;
    for(var i = 0;i < matA.cols();i += 1) {
      cb(matA.getUnsafe(posA), matB.getUnsafe(posB), row, i, col);
      posA += 1;
      posB += matB.cols();
    }
  };
  Matrix.mul = function(matA, matB) {
    if(matA.cols() !== matB.rows()) {
      console.warn("incompatible dimensions", matA.rows() + "x" + matA.cols(), matB.rows() + "x" + matB.cols());
      return null;
    }
    var mat = matA.createArray(matA.rows(), matB.cols());
    var pos = 0;
    for(var r = 0;r < matA.rows();r += 1) {
      for(var c = 0;c < matB.cols();c += 1) {
        var sum = 0;
        Matrix.iter(matA, matB, r, c, function(a, b, r, i, c) {
          sum += a * b;
        });
        mat[pos] = sum;
        pos += 1;
      }
    }
    return new Matrix(mat, matA.rows(), matB.cols());
  };
  Matrix.add = function(matA, matB) {
    if(matA.rows() !== matB.rows() || matA.cols() !== matB.cols()) {
      console.warn("incompatible dimensions", matA.rows() + "x" + matA.cols(), matB.rows() + "x" + matB.cols());
      return null;
    }
    var mat = matA.createArray(matA.rows(), matA.cols());
    for(var pos = 0;pos < mat.length;pos += 1) {
      mat[pos] = matA.getUnsafe(pos) + matB.getUnsafe(pos);
    }
    return new Matrix(mat, matA.rows(), matA.cols());
  };

}; // mdsjs

mdsjs = new mdsjs(); // create instance
