# point_in_polygon
This project implements and benchmarks two geometric algorithms for testing whether a point lies inside a **convex polygon** in ℝ² (2D Euclidean space), with a strong focus on **performance** and **algorithmic complexity**.

Developed for evaluating and comparing the computational cost of:
- A naive **O(N)** method (sequential side test)
- An optimized **O(log N)** method using binary search and projective geometry
  
---

## Project Structure

- `measure_time.cpp`  
  Compares the execution time of both point-in-polygon algorithms over increasing polygon sizes (from 10 to 1000 vertices). Outputs results to a CSV file.

- `measure_ratio.cpp`  
  Measures how the performance of the `O(log N)` method changes as the number of test points increases. Also outputs results to CSV.

- Internal utilities:
  - Random convex polygon generation
  - Projective geometry primitives (`Point`, `Line`, homogeneous coordinates)
  - Binary search–based point inclusion test
  - Timing via `std::chrono`

---

## Algorithms

### `O(N)` Naive Inclusion Test
Iterates over all polygon edges and uses the **orientation test** (via projective representation) to check whether the point lies on the "inner" side of each edge.

### `O(log N)` Optimized Test
Assumes:
- The polygon is **strictly convex**
- Vertices are ordered counterclockwise

It uses:
- A **binary search** over the polygon wedges
- Centroid-based ray casting
- Orientation via **projective duality**

---

## How to Use

### Compile

```bash
g++ -std=c++11 -O2 measure_time.cpp -o measure_time
g++ -std=c++11 -O2 measure_ratio.cpp -o measure_ratio
