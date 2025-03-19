#!/usr/bin/env julia

"""
Test script for the optimized Vietoris-Rips implementation.
"""

include("vietoris_rips_optimized.jl")
using .VietorisRipsOptimized
using LinearAlgebra
using Plots
using Random

function main()
    # Set random seed for reproducibility
    Random.seed!(42)
    
    # Generate a larger point cloud (random points on a torus)
    n_points = 200
    points = generate_torus_points(n_points, 2.0, 0.5)
    
    println("Point cloud size: $n_points points in 3D")
    
    # Benchmark time for computing Rips complex
    epsilon = 0.5
    max_dim = 2
    
    println("\nComputing Vietoris-Rips complex with ε = $epsilon, max_dim = $max_dim...")
    @time simplices = rips_complex(points, epsilon, max_dim)
    
    # Count simplices by dimension
    counts = Dict{Int,Int}()
    for simplex in simplices
        dim = length(simplex) - 1
        counts[dim] = get(counts, dim, 0) + 1
    end
    
    println("\nResult:")
    println("  Total simplices: $(length(simplices))")
    for dim in sort(collect(keys(counts)))
        println("  Dimension $dim: $(counts[dim]) simplices")
    end
    
    # Compute birth times for persistence analysis
    println("\nComputing persistence data...")
    @time birth_times = rips_filtration(points, 1.0, max_dim)
    
    # Compute 0-dimensional persistence diagram
    println("\nComputing persistence diagram...")
    @time intervals = compute_persistence(birth_times)
    
    # Display persistence intervals
    println("\nPersistence intervals (birth, death):")
    for (i, (birth, death)) in enumerate(intervals)
        death_str = isinf(death) ? "∞" : round(death, digits=3)
        println("  $i: ($(round(birth, digits=3)), $death_str)")
        if i >= 10  # Just show the first 10 intervals
            println("  ...")
            break
        end
    end
    
    # Plot points in 3D
    plt = scatter(
        points[1, :], points[2, :], points[3, :],
        title="Point Cloud (Torus)",
        label="Points",
        marker=:circle,
        markersize=2,
        legend=:none,
        camera=(30, 30)
    )
    
    # Save or display the plot
    savefig(plt, "torus_points.png")
    println("\nPlot saved as 'torus_points.png'")
end

"""
Generate points on a torus in 3D space.
"""
function generate_torus_points(n, R, r)
    points = zeros(3, n)
    
    for i in 1:n
        θ = 2π * rand()  # Angle around the center of the tube
        φ = 2π * rand()  # Angle around the center of the torus
        
        # Parametric equations for a torus
        points[1, i] = (R + r * cos(θ)) * cos(φ)
        points[2, i] = (R + r * cos(θ)) * sin(φ)
        points[3, i] = r * sin(θ)
    end
    
    return points
end

main() 