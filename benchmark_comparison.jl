#!/usr/bin/env julia

"""
Benchmark comparison between the basic and optimized Vietoris-Rips implementations.
This script compares the performance of both implementations on point clouds of different sizes.
"""

include("vietoris_rips.jl")
include("vietoris_rips_optimized.jl")
using .VietorisRips
using .VietorisRipsOptimized
using LinearAlgebra
using Random
using Printf

function main()
    # Set random seed for reproducibility
    Random.seed!(42)
    
    # Parameters
    epsilon = 0.3
    max_dim = 2
    
    # Test different point cloud sizes
    point_counts = [50, 100, 150]
    
    println("Benchmarking Vietoris-Rips implementations")
    println("==========================================")
    println("Parameters:")
    println("  Epsilon: $epsilon")
    println("  Max dimension: $max_dim")
    println("==========================================")
    
    # Print header
    @printf("%-10s %-15s %-15s %-15s %-15s\n", 
            "Points", "Basic (s)", "Optimized (s)", "Basic count", "Optimized count")
    println("-" ^ 70)
    
    for n_points in point_counts
        # Generate random points in a unit cube
        points = rand(3, n_points)
        
        # Benchmark basic implementation
        basic_time = @elapsed begin
            basic_simplices = VietorisRips.rips_complex(points, epsilon, max_dim)
        end
        
        # Benchmark optimized implementation
        optimized_time = @elapsed begin
            optimized_simplices = VietorisRipsOptimized.rips_complex(points, epsilon, max_dim)
        end
        
        # Count simplices
        basic_count = length(basic_simplices)
        optimized_count = length(optimized_simplices)
        
        # Print results
        @printf("%-10d %-15.6f %-15.6f %-15d %-15d\n", 
                n_points, basic_time, optimized_time, basic_count, optimized_count)
                
        # Verify results match
        if basic_count != optimized_count
            println("Warning: Simplex counts don't match for n=$n_points!")
        end
    end
    
    println("\nBenchmarking filtration computation")
    println("==========================================")
    
    # Test filtration on the largest point cloud
    n_points = point_counts[end]
    points = rand(3, n_points)
    max_epsilon = 0.5
    
    # Benchmark basic filtration
    basic_time = @elapsed begin
        basic_filtration = VietorisRips.rips_filtration(points, max_epsilon, 5, max_dim)
    end
    
    # Benchmark optimized filtration
    optimized_time = @elapsed begin
        optimized_filtration = VietorisRipsOptimized.rips_filtration(points, max_epsilon, max_dim)
    end
    
    # Print results
    println("Filtration with $n_points points:")
    @printf("  Basic implementation: %.6f seconds\n", basic_time)
    @printf("  Optimized implementation: %.6f seconds\n", optimized_time)
    @printf("  Speedup: %.2fx\n", basic_time / optimized_time)
end

main() 