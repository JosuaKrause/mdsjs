MDS for JavaScript
==================

A multi dimensional scaling library for JavaScript. Only a distance matrix
for elements is required to position them on a 2D plane. This can be useful
for visualizing high-dimensional data, drawing arbitrary graphs, or projecting
points using a custom distance function.

**This is a work in progress and currently only computes MDS without most of the optimizations in the paper!**

Refer to [index.html](index.html) for a small [example](https://josuakrause.github.io/mdsjs/) of how to use the library.
Pull requests are highly appreciated.

The algorithm used for calculating the MDS is published at

```
@incollection{
  year = {2007},
  isbn = {978-3-540-70903-9},
  booktitle = {Graph Drawing},
  volume = {4372},
  series = {Lecture Notes in Computer Science},
  editor = {Kaufmann, Michael and Wagner, Dorothea},
  doi = {10.1007/978-3-540-70904-6_6},
  title = {Eigensolver Methods for Progressive Multidimensional Scaling of Large Data},
  url = {http://dx.doi.org/10.1007/978-3-540-70904-6_6},
  publisher = {Springer Berlin Heidelberg},
  author = {Brandes, Ulrik and Pich, Christian},
  pages = {42-53},
  language = {English}
}
```

A Java implementation of the paper can be found at [mdsj] (http://www.inf.uni-konstanz.de/algo/software/mdsj/).
