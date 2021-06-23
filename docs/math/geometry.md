# geometry.hpp

It contains only elementary methods for geometry

## Geomerty

It is a namespace contains following methods

``` cpp
using Point = std::pair<double, double>;
double cross(const Point &op, const Point &sp, const Point &ep) // cross product
bool crossLeft(const Point &op, const Point &sp, const Point &ep) // check if cross product is positive
double dist2(const Point &p, const Point &q) // square of distance
double dist (const Point& p, const Point &q) // dist
std::vector<Point> convexHull(std::vector<Point> p) // convexHull of p
double diameter(std::vector<Point> p) // find greatest distance in between Vertices of p
double minDist(std::vector<Point> a)  // find smallest distance in between Vertices of a
```

__Complexity__:

- cross, crossLeft, dist, dist2:  $O(1)$
- convexHull, diameter: $O(n \log n)$
- minDist: $O(n \log^2 n)$

## partialOrder

partial order of dimension $k$ (optimed by bitset)

__Complexity__:

- $O(\frac{n^2}{w})$ where $w = 32,\;64$

