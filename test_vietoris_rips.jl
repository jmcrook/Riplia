#!/usr/bin/env julia

"""
Unit tests for the Vietoris-Rips implementation.
"""

include("vietoris_rips.jl")
using .VietorisRips
using Test
using LinearAlgebra

@testset "VietorisRips" begin
    # Test 1: Simple square
    @testset "Square" begin
        # Four corners of a unit square
        points = [0.0 1.0 0.0 1.0; 
                  0.0 0.0 1.0 1.0]
        
        # With epsilon = 1.0, we should get 4 vertices, 4 edges, and no triangles
        # (diagonals have length √2 > 1)
        simplices = rips_complex(points, 1.0, 2)
        
        # Count simplices by dimension
        counts = Dict{Int,Int}()
        for simplex in simplices
            dim = length(simplex) - 1
            counts[dim] = get(counts, dim, 0) + 1
        end
        
        @test counts[0] == 4  # 4 vertices
        @test counts[1] == 4  # 4 edges (sides of the square)
        @test get(counts, 2, 0) == 0  # No triangles
        
        # With epsilon = 1.5, we should get triangles as well
        simplices = rips_complex(points, 1.5, 2)
        
        # Count simplices by dimension
        counts = Dict{Int,Int}()
        for simplex in simplices
            dim = length(simplex) - 1
            counts[dim] = get(counts, dim, 0) + 1
        end
        
        @test counts[0] == 4  # 4 vertices
        @test counts[1] == 6  # 6 edges (sides + diagonals)
        @test counts[2] == 4  # 4 triangles
    end
    
    # Test 2: Circle
    @testset "Circle" begin
        # Generate points on a circle
        n_points = 8
        θ = range(0, 2π, length=n_points+1)[1:end-1]
        points = [cos.(θ)'; sin.(θ)']
        
        # With a small epsilon, we should only get edges between adjacent points
        small_eps = 0.8
        simplices = rips_complex(points, small_eps, 2)
        
        # Count simplices by dimension
        counts = Dict{Int,Int}()
        for simplex in simplices
            dim = length(simplex) - 1
            counts[dim] = get(counts, dim, 0) + 1
        end
        
        @test counts[0] == n_points  # n_points vertices
        @test counts[1] == n_points  # n_points edges (the circle)
        @test get(counts, 2, 0) == 0  # No triangles
        
        # With a larger epsilon, we should get more edges and triangles
        large_eps = 1.5
        simplices = rips_complex(points, large_eps, 2)
        
        # Count simplices by dimension
        counts = Dict{Int,Int}()
        for simplex in simplices
            dim = length(simplex) - 1
            counts[dim] = get(counts, dim, 0) + 1
        end
        
        @test counts[0] == n_points  # n_points vertices
        @test counts[1] > n_points   # More than n_points edges
        @test get(counts, 2, 0) > 0  # Some triangles
    end
    
    # Test 3: Filtration test
    @testset "Filtration" begin
        # Four corners of a unit square
        points = [0.0 1.0 0.0 1.0; 
                  0.0 0.0 1.0 1.0]
        
        # Test filtration with multiple steps
        filtration = rips_filtration(points, 2.0, 5, 2)
        
        # Epsilon values should be evenly spaced
        epsilons = [eps for (eps, _) in filtration]
        @test length(epsilons) == 5
        @test epsilons[1] ≈ 0.0
        @test epsilons[5] ≈ 2.0
        
        # Number of simplices should be non-decreasing
        simplex_counts = [length(simplices) for (_, simplices) in filtration]
        for i in 1:length(simplex_counts)-1
            @test simplex_counts[i] <= simplex_counts[i+1]
        end
        
        # Last step should have all simplices (tetrahedron)
        _, last_simplices = filtration[end]
        counts = Dict{Int,Int}()
        for simplex in last_simplices
            dim = length(simplex) - 1
            counts[dim] = get(counts, dim, 0) + 1
        end
        
        @test counts[0] == 4  # 4 vertices
        @test counts[1] == 6  # 6 edges (complete graph K4)
        @test counts[2] == 4  # 4 triangles
    end
end

println("All tests passed!") 