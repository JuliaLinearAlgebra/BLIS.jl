using Test
using Random
using LinearAlgebra
using LinearAlgebra: BLAS
using InteractiveUtils: @which
using DelimitedFiles
using Statistics
using BLIS

Random.seed!(1234)

rtype(::Type{Complex{T}}) where {T} = T

zrtest(val, atol, label) = begin
    iszr = ≈(val, 0.0, atol=atol)
    if !iszr
        @info "`$label` test failed (err=$val). Consider adding it to ~/.blis_jlbla_blacklist."
    end
    return iszr
end

"Run & evict this method."
macro run_evict(func, io, largs)
    return quote
        # Invokation.
        $func($(esc(largs))..., $(esc(io)))

        # Evict.
        method_blis = @which $func($(esc(largs))..., $(esc(io)))
        Base.delete_method(method_blis)
    end
end

include("test_level1.jl")
include("test_level2.jl")
include("test_level3.jl")

