# multidimensional median

There are two slightly different generalizations of the median in higher dimensions.

**The centerpoint**, a point that, no matter which projection has at most ⌈n/2⌉ points on the left and at most ⌈n/2⌉ points on its right. There is an algorithm for finding the centerpoint in linear time (Jadhav & Mukhopadhyay (1994)) which seemed to me to be overly complex. I posit that quickselect can be generalized for higher dimensions also. ```fast_centerpoint``` is an implementation of quickselect in higher dimensions.

**The geometric median** is the point that has the smallest sum of absolute distances to all other points. The centerpoint cannot be far off, so using a linear time result to jumpstart an optimization procedure will yield a very fast convergence. Especially because, if the geometric median is not coincidental with one of the points, the sum of absolute distances in higher dimensions is smooth and monotonic. This should give quadratic complexity, and thereby for a required precision ε, the centerpoint of n points in m dimensions will already be about n^(-1/m) close to the geometric median, and Newton's Method will roughly double the precision in each step, n^(-1/m) , n^(-2/m) , n^(-4/m) ... and so on. Therefore the number of steps required is
ε = n^(-s/m) <-> s = -log(ε)*m/log(n)
which is O(log(1/ε)) or thereabouts, similar to Cohen et al. (2016).

Jadhav & Mukhopadhyay (1994) "Computing a centerpoint of a finite planar set of points in linear time", Discrete and Computational Geometry, [doi:10.1007/BF02574382](https://doi.org/10.1007%2FBF02574382).

Cohen, Michael; Lee, Yin Tat; Miller, Gary; Pachocki, Jakub; Sidford, Aaron (2016). "Geometric median in nearly linear time". Proc. 48th Symposium on Theory of Computing (STOC 2016). Association for Computing Machinery. arXiv:1606.05225. [doi:10.1145/2897518.2897647](https://doi.org/10.1145%2F2897518.2897647).
