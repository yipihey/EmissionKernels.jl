# EmissionKernels test suite. Run standalone (no macOS C-grackle oracle needed):
#   <julia> --project=lib/EmissionKernels/test lib/EmissionKernels/test/runtests.jl
using EmissionKernels
using Test

try
    @eval using Metal
catch
end

@testset "EmissionKernels" begin
    include("test_emission.jl")
end

@testset "sub-5K f32: H2 low/LTE interpolation is 0/0-guarded" begin
    # Below ~5 K BOTH H2LTE(T) and the Galli–Palla low-density collision sum
    # underflow Float32 → the unguarded h2lte/galdl was 0/0 = NaN (hit by the
    # first BlockAMR CICASS run to pass z<15; cold-void cells reach T≈4 K).
    for T in (1.5f0, 2.0f0, 3.0f0, 4.0f0, 4.9f0)
        v = EmissionKernels.cooling_rate_total(3.4f-4, 5.1f-8, 2.6f-5, 5.1f-8,
                                               7.7f-10, 0.0f0, T, 14.9f0)
        @test isfinite(v)
    end
end
