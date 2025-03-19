#!/usr/bin/env julia

"""
Example usage of the VietorisRips module.
This script demonstrates how to use the Vietoris-Rips algorithm for
computing simplicial complexes from point cloud data.
"""

# Include the VietorisRips module
include("vietoris_rips.jl")
using .VietorisRips
using LinearAlgebra
using Plots

function main()
    # Generate a simple point cloud (a circle with some noise)
    n_points = 20
    θ = range(0, 2π, length=n_points)
    points = [cos.(θ)'; sin.(θ)']
    
    # Add some noise
    points .+= 0.05 .* randn(size(points))
    
    # Compute Vietoris-Rips complex at a specific epsilon
    epsilon = 0.5
    max_dim = 2
    simplices = rips_complex(points, epsilon, max_dim)
    
    println("Number of simplices: ", length(simplices))
    for (dim, count) in countmap(length.(simplices) .- 1)
        println("  Dimension $dim: $count simplices")
    end
    
    # Visualize the 1-skeleton of the complex
    p = scatter(points[1,:], points[2,:], 
                label="Points", aspect_ratio=:equal, 
                markersize=6, markershape=:circle,
                title="Vietoris-Rips Complex (ε=$epsilon)")
    
    # Draw the edges (1-simplices)
    for simplex in simplices
        if length(simplex) == 2  # It's an edge
            i, j = simplex
            plot!([points[1,i], points[1,j]], [points[2,i], points[2,j]], 
                  color=:blue, linewidth=1.5, label="")
        end
    end
    
    display(p)
    
    # Compute a filtration
    filtration = rips_filtration(points, 1.0, 5, max_dim)
    
    println("\nFiltration:")
    for (eps, simplices) in filtration
        println("  ε = $eps: $(length(simplices)) simplices")
    end
end

# Helper function to count occurrences
function countmap(x)
    result = Dict{eltype(x),Int}()
    for val in x
        result[val] = get(result, val, 0) + 1
    end
    return result
end

main() 