# Vietoris-Rips Algorithm in Julia

A frugal implementation of the Vietoris-Rips algorithm for topological data analysis in Julia.

## Overview

The Vietoris-Rips algorithm is a fundamental tool in topological data analysis (TDA) that constructs a simplicial complex from a point cloud based on pairwise distances between points. This implementation provides:

1. A basic version (`vietoris_rips.jl`) that is simple and easy to understand
2. An optimized version (`vietoris_rips_optimized.jl`) that uses sparse matrices for better memory efficiency
3. An example script (`vietoris_rips_example.jl`) that demonstrates usage
4. Benchmark utilities (`benchmark_comparison.jl`) for performance comparison
5. Tests for larger point clouds (`test_optimized.jl`)

## Features

- Build Vietoris-Rips complexes from point cloud data
- Compute filtrations over a range of distance thresholds
- Support for higher-dimensional simplices
- Visualization of the 1-skeleton of the complex
- Basic persistence diagram computation (for 0-dimensional homology)
- Performance benchmarking tools

## Requirements

- Julia 1.0 or later
- Packages:
  - `LinearAlgebra` (standard library)
  - `SparseArrays` (standard library)
  - `DataStructures` (for the optimized version)
  - `Plots` (for visualization in the example)
  - `Random` (for generating test data)
  - `Printf` (for benchmark output formatting)
  - `Test` (for unit tests)

## Installation

Clone this repository or download the files:

```bash
git clone https://github.com/jmcrook/riplia.git
cd riplia

# Install required packages
julia -e 'using Pkg; Pkg.add(["DataStructures", "Plots", "Test", "Random", "Printf"])'
```

## Usage

### Basic Usage

```julia
include("vietoris_rips.jl")
using .VietorisRips
using LinearAlgebra

# Create a point cloud (each column is a point)
points = randn(2, 10)  # 10 random points in 2D

# Compute Vietoris-Rips complex with epsilon = 0.5 and max dimension = 2
simplices = rips_complex(points, 0.5, 2)

# Compute a filtration with 10 steps, from 0 to 1.0
filtration = rips_filtration(points, 1.0, 10, 2)
```

### Optimized Version

```julia
include("vietoris_rips_optimized.jl")
using .VietorisRipsOptimized
using LinearAlgebra

# Create a point cloud (each column is a point)
points = randn(2, 100)  # 100 random points in 2D

# Compute Vietoris-Rips complex 
simplices = rips_complex(points, 0.5, 2)

# Compute birth times of simplices
birth_times = rips_filtration(points, 1.0, 2)

# Compute 0-dimensional persistence diagram
persistence_intervals = compute_persistence(birth_times)
```

### Run the Examples

```bash
# Basic example with a circle of 20 points
julia vietoris_rips_example.jl

# Larger example with a torus of 200 points
julia test_optimized.jl

# Performance comparison between implementations
julia benchmark_comparison.jl

# Run unit tests
julia test_vietoris_rips.jl
```

## Test Results

We've run several tests to validate the implementation:

1. **Unit Tests**: All 22 tests passed, confirming the correctness of the algorithm on simple geometric examples.

2. **Example with 20 points (Circle)**: 
   - For ε = 0.5 and max_dim = 2, we got 46 simplices:
     - 20 0-dimensional simplices (vertices)
     - 23 1-dimensional simplices (edges)
     - 3 2-dimensional simplices (triangles)

3. **Larger Test (Torus with 200 points)**:
   - For ε = 0.5 and max_dim = 2, we got 894 simplices:
     - 200 0-dimensional simplices (vertices)
     - 415 1-dimensional simplices (edges)
     - 279 2-dimensional simplices (triangles)

## Benchmark Results

Performance comparison between basic and optimized implementations:

| Points | Basic (s) | Optimized (s) | Speedup |
|--------|-----------|---------------|---------|
| 50     | 0.000775  | 0.030727      | 0.03x   |
| 100    | 0.017136  | 0.001031      | 16.6x   |
| 150    | 0.067096  | 0.002688      | 25.0x   |

**Filtration computation** (150 points):
- Basic implementation: 1.268272 seconds
- Optimized implementation: 0.050350 seconds
- **Speedup: 25.2x**

These results demonstrate that:
1. For very small point clouds, the basic implementation may be faster due to lower overhead
2. As the point cloud size increases, the optimized implementation becomes significantly faster
3. For filtration computation, the optimized implementation offers substantial performance improvements

## Algorithm Details

The Vietoris-Rips algorithm works as follows:

1. Construct a graph where vertices are points and edges connect points within distance ε of each other
2. Add all simplices (vertices, edges, triangles, etc.) where all constituent vertices are pairwise connected
3. For filtrations, vary ε and track when simplices appear

## Performance Considerations

- For small point clouds (< 50 points), the basic implementation may be sufficient
- For larger datasets (> 100 points) or higher-dimensional computations, use the optimized version
- Computation time grows exponentially with the maximum dimension of simplices
- Memory usage can be high for large point clouds with many edges
- The optimized implementation shows a 25x speedup for point clouds of 150+ points

## License

MIT

## References

- Edelsbrunner, H., & Harer, J. (2010). Computational Topology: An Introduction.
- Carlsson, G. (2009). Topology and data. Bulletin of the American Mathematical Society, 46(2), 255-308. 