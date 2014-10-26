/**
 * Created by krause on 2014-10-25.
 */

mdsjs = function() {

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
    if(k < m.rows()) {
      console.warn("requested more pivots than elements", k, m.rows(), m.cols());
      return null;
    }
    var mat = this.createArray(k, m.cols());
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
    this.isSquare = function() {
      return rows === cols;
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
