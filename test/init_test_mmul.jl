# Simple test on matrix multiplication.
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\vvvvvv/!!!!!!!!!!!!!!!!!!!!
# !CAUTION: THIS TEST MUST OCCUR >BEFORE< BLIS WAS IMPORTED.!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/^^^^^^\!!!!!!!!!!!!!!!!!!!!
#
using Test
using Random
using LinearAlgebra
using LinearAlgebra: BLAS
using Statistics
Random.seed!(1234)
rtype(::Type{Complex{T}}) where {T} = T
zrtest(val, atol, label) = begin
    iszr = ≈(val, 0.0, atol=atol)
    if !iszr
        @info "`$label` test failed. Consider adding it to ~/.blis_jlbla_blacklist."
    end
    return iszr
end

@testset "BLAS level-3 LinearAlgebra interface" begin
αr = 1.1
βr = 1.1
αc = 1.1 + 0.2im
βc = 1.1 + 0.3im

χlarge = 500
χsmall = 20

Alarge_base = rand(ComplexF64, χlarge, χlarge)
Blarge_base = rand(ComplexF64, χlarge, χlarge)
Clarge_base = rand(ComplexF64, χlarge, χlarge)

Asmall_base = rand(ComplexF64, χsmall, χsmall)
Bsmall_base = rand(ComplexF64, χsmall, χsmall)
Csmall_base = rand(ComplexF64, χsmall, χsmall)

global Clarge_gemm = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_hemm = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_symm = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_her2k= [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_syr2k= [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_herk = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Clarge_syrk = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Clarge_trmm = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Clarge_trsm = [zeros(T, χlarge, χlarge) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_gemm = [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_hemm = [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_symm = [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_her2k= [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_syr2k= [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_herk = [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Csmall_syrk = [zeros(T, χsmall, χsmall) for T=(Float32, Float64, ComplexF32, ComplexF64)]

global Cst_lg_gemm = [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Cst_lg_hemm = [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Cst_lg_symm = [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Cst_lg_her2k= [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Cst_lg_syr2k= [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Cst_lg_herk = [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
# global Cst_lg_syrk = [zeros(T, χlarge÷2, χlarge÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Cst_sm_gemm = [zeros(T, χsmall÷2, χsmall÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Cst_sm_hemm = [zeros(T, χsmall÷2, χsmall÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]
global Cst_sm_symm = [zeros(T, χsmall÷2, χsmall÷2) for T=(Float32, Float64, ComplexF32, ComplexF64)]

for (i, T)=zip(1:4, [Float32, Float64, ComplexF32, ComplexF64])
    Alarge = zeros(T, χlarge, χlarge)
    Blarge = zeros(T, χlarge, χlarge)
    Asmall = zeros(T, χsmall, χsmall)
    Bsmall = zeros(T, χsmall, χsmall)

    local elconv, αu, βu, αR, βR
    local locl_hemm!, locl_her2k!, locl_herk!
    if eltype(Alarge)<:Complex
        elconv = x -> x
        αu = T(αc)
        βu = T(βc)
        αR = rtype(T)(αr)
        βR = rtype(T)(βr)
        locl_hemm!  = BLAS.hemm!
        locl_herk!  = BLAS.herk!
        locl_her2k! = BLAS.her2k!
    else
        elconv = x -> real(x)
        αu = T(αr)
        βu = T(βr)
        αR = αu
        βR = βu
        locl_hemm!  = BLAS.symm!
        locl_herk!  = BLAS.syrk!
        locl_her2k! = BLAS.syr2k!
    end
    Alarge .= elconv.(Alarge_base)
    Blarge .= elconv.(Blarge_base)
    Asmall .= elconv.(Asmall_base)
    Bsmall .= elconv.(Bsmall_base)
    Clarge_gemm[i]  .= elconv.(Clarge_base)
    Clarge_hemm[i]  .= elconv.(Clarge_base)
    Clarge_symm[i]  .= elconv.(Clarge_base)
    Clarge_her2k[i] .= elconv.(Clarge_base)
    Clarge_syr2k[i] .= elconv.(Clarge_base)
    Clarge_herk[i]  .= elconv.(Clarge_base)
    Clarge_syrk[i]  .= elconv.(Clarge_base)
    Csmall_gemm[i]  .= elconv.(Csmall_base)
    Csmall_hemm[i]  .= elconv.(Csmall_base)
    Csmall_symm[i]  .= elconv.(Csmall_base)
    Csmall_her2k[i] .= elconv.(Csmall_base)
    Csmall_syr2k[i] .= elconv.(Csmall_base)
    Csmall_herk[i]  .= elconv.(Csmall_base)
    Csmall_syrk[i]  .= elconv.(Csmall_base)

    # Strided.
    Ast_lg = view(Alarge, 1:2:χlarge, 1:2:χlarge)
    Bst_lg = view(Blarge, 1:2:χlarge, 1:2:χlarge)
    Ast_sm = view(Asmall, 1:2:χsmall, 1:2:χsmall)
    Bst_sm = view(Bsmall, 1:2:χsmall, 1:2:χsmall)

    # Execute: column-major.
    BLAS.gemm!('N', 'N', αu, Alarge, Blarge, βu, Clarge_gemm[i])
    BLAS.gemm!('N', 'N', αu, Asmall, Bsmall, βu, Csmall_gemm[i])
    locl_hemm!('L', 'U', αu, Alarge, Blarge, βu, Clarge_hemm[i])
    locl_hemm!('R', 'U', αu, Asmall, Bsmall, βu, Csmall_hemm[i])
    BLAS.symm!('L', 'L', αu, Alarge, Blarge, βu, Clarge_symm[i])
    BLAS.symm!('L', 'L', αu, Asmall, Bsmall, βu, Csmall_symm[i])
    locl_her2k!('U', 'N', αu, Alarge, Blarge, βR, Clarge_her2k[i])
    locl_her2k!('U', 'N', αu, Asmall, Bsmall, βR, Csmall_her2k[i])
    BLAS.syr2k!('U', 'N', αu, Alarge, Blarge, βu, Clarge_syr2k[i])
    BLAS.syr2k!('U', 'N', αu, Asmall, Bsmall, βu, Csmall_syr2k[i])
    locl_herk!('U', 'N', αR, Alarge, βR, Clarge_herk[i])
    locl_herk!('U', 'N', αR, Asmall, βR, Csmall_herk[i])
    BLAS.syrk!('U', 'N', αu, Alarge, βu, Clarge_syrk[i])
    BLAS.syrk!('U', 'N', αu, Asmall, βu, Csmall_syrk[i])

    # Execute: generic-strided.
    Cst_lg_gemm[i] .= Ast_lg * Bst_lg
    Cst_sm_gemm[i] .= Ast_sm * Bst_sm
    Cst_lg_hemm[i] .= Hermitian(Array(Ast_lg)) * Array(Bst_lg)
    Cst_sm_hemm[i] .= Hermitian(Array(Ast_sm)) * Array(Bst_sm)
    Cst_lg_symm[i] .= Symmetric(Array(Ast_lg)) * Array(Bst_lg)
    Cst_sm_symm[i] .= Symmetric(Array(Ast_sm)) * Array(Bst_sm)
end

# Import BLIS for testing.
using BLIS
if length(BLIS.BLASInterface.blacklist) > 0
    @info "Blacklisted methods: $(BLIS.BLASInterface.blacklist)."
end

for (i, T)=zip(1:4, [Float32, Float64, ComplexF32, ComplexF64])
    Alarge = zeros(T, χlarge, χlarge)
    Blarge = zeros(T, χlarge, χlarge)
    Asmall = zeros(T, χsmall, χsmall)
    Bsmall = zeros(T, χsmall, χsmall)

    local elconv, αu, βu, αR, βR
    local locl_hemm!, locl_her2k!, locl_herk!
    if eltype(Alarge)<:Complex
        elconv = x -> x
        αu = T(αc)
        βu = T(βc)
        αR = rtype(T)(αr)
        βR = rtype(T)(βr)
        locl_hemm!  = BLAS.hemm!
        locl_herk!  = BLAS.herk!
        locl_her2k! = BLAS.her2k!
    else
        elconv = x -> real(x)
        αu = T(αr)
        βu = T(βr)
        αR = αu
        βR = βu
        locl_hemm!  = BLAS.symm!
        locl_herk!  = BLAS.syrk!
        locl_her2k! = BLAS.syr2k!
    end
    Alarge .= elconv.(Alarge_base)
    Blarge .= elconv.(Blarge_base)
    Asmall .= elconv.(Asmall_base)
    Bsmall .= elconv.(Bsmall_base)
    Clarge_gemm_cur  = T.(elconv.(Clarge_base))
    Clarge_hemm_cur  = T.(elconv.(Clarge_base))
    Clarge_symm_cur  = T.(elconv.(Clarge_base))
    Clarge_her2k_cur = T.(elconv.(Clarge_base))
    Clarge_syr2k_cur = T.(elconv.(Clarge_base))
    Clarge_herk_cur  = T.(elconv.(Clarge_base))
    Clarge_syrk_cur  = T.(elconv.(Clarge_base))
    Csmall_gemm_cur  = T.(elconv.(Csmall_base))
    Csmall_hemm_cur  = T.(elconv.(Csmall_base))
    Csmall_symm_cur  = T.(elconv.(Csmall_base))
    Csmall_her2k_cur = T.(elconv.(Csmall_base))
    Csmall_syr2k_cur = T.(elconv.(Csmall_base))
    Csmall_herk_cur  = T.(elconv.(Csmall_base))
    Csmall_syrk_cur  = T.(elconv.(Csmall_base))

    # Strided.
    Ast_lg = view(Alarge, 1:2:χlarge, 1:2:χlarge)
    Bst_lg = view(Blarge, 1:2:χlarge, 1:2:χlarge)
    Ast_sm = view(Asmall, 1:2:χsmall, 1:2:χsmall)
    Bst_sm = view(Bsmall, 1:2:χsmall, 1:2:χsmall)

    # Execute: column-major.
    BLAS.gemm!('N', 'N', αu, Alarge, Blarge, βu, Clarge_gemm_cur)
    BLAS.gemm!('N', 'N', αu, Asmall, Bsmall, βu, Csmall_gemm_cur)
    locl_hemm!('L', 'U', αu, Alarge, Blarge, βu, Clarge_hemm_cur)
    locl_hemm!('R', 'U', αu, Asmall, Bsmall, βu, Csmall_hemm_cur)
    BLAS.symm!('L', 'L', αu, Alarge, Blarge, βu, Clarge_symm_cur)
    BLAS.symm!('L', 'L', αu, Asmall, Bsmall, βu, Csmall_symm_cur)
    locl_her2k!('U', 'N', αu, Alarge, Blarge, βR, Clarge_her2k_cur)
    locl_her2k!('U', 'N', αu, Asmall, Bsmall, βR, Csmall_her2k_cur)
    BLAS.syr2k!('U', 'N', αu, Alarge, Blarge, βu, Clarge_syr2k_cur)
    BLAS.syr2k!('U', 'N', αu, Asmall, Bsmall, βu, Csmall_syr2k_cur)
    locl_herk!('U', 'N', αR, Alarge, βR, Clarge_herk_cur)
    locl_herk!('U', 'N', αR, Asmall, βR, Csmall_herk_cur)
    BLAS.syrk!('U', 'N', αu, Alarge, βu, Clarge_syrk_cur)
    BLAS.syrk!('U', 'N', αu, Asmall, βu, Csmall_syrk_cur)

    # Execute: generic-strided.
    Cst_lg_gemm_cur = Ast_lg * Bst_lg
    Cst_sm_gemm_cur = Ast_sm * Bst_sm
    Cst_lg_hemm_cur = Hermitian(Ast_lg) * Bst_lg
    Cst_sm_hemm_cur = Hermitian(Ast_sm) * Bst_sm
    Cst_lg_symm_cur = Symmetric(Ast_lg) * Bst_lg
    Cst_sm_symm_cur = Symmetric(Ast_sm) * Bst_sm

    # Check.
    @test zrtest(mean(abs.(Clarge_gemm_cur  - Clarge_gemm[i] )), 1e-6*χlarge^1.2, "gemm")
    @test zrtest(mean(abs.(Clarge_hemm_cur  - Clarge_hemm[i] )), 1e-6*χlarge^1.2, "hemm")
    @test zrtest(mean(abs.(Clarge_symm_cur  - Clarge_symm[i] )), 1e-6*χlarge^1.2, "symm")
    @test zrtest(mean(abs.(Clarge_her2k_cur - Clarge_her2k[i])), 1e-6*χlarge^1.2, "her2k")
    @test zrtest(mean(abs.(Clarge_syr2k_cur - Clarge_syr2k[i])), 1e-6*χlarge^1.2, "syr2k")
    @test zrtest(mean(abs.(Clarge_herk_cur  - Clarge_herk[i] )), 1e-6*χlarge^1.2, "herk")
    @test zrtest(mean(abs.(Clarge_syrk_cur  - Clarge_syrk[i] )), 1e-6*χlarge^1.2, "syrk")
    @test zrtest(mean(abs.(Csmall_gemm_cur  - Csmall_gemm[i] )), 1e-6*χsmall^1.2, "gemm")
    @test zrtest(mean(abs.(Csmall_hemm_cur  - Csmall_hemm[i] )), 1e-6*χsmall^1.2, "hemm")
    @test zrtest(mean(abs.(Csmall_symm_cur  - Csmall_symm[i] )), 1e-6*χsmall^1.2, "symm")
    @test zrtest(mean(abs.(Csmall_her2k_cur - Csmall_her2k[i])), 1e-6*χsmall^1.2, "her2k")
    @test zrtest(mean(abs.(Csmall_syr2k_cur - Csmall_syr2k[i])), 1e-6*χsmall^1.2, "syr2k")
    @test zrtest(mean(abs.(Csmall_herk_cur  - Csmall_herk[i] )), 1e-6*χsmall^1.2, "herk")
    @test zrtest(mean(abs.(Csmall_syrk_cur  - Csmall_syrk[i] )), 1e-6*χsmall^1.2, "syrk")

    # Check - strided.
    @test zrtest(mean(abs.(Cst_lg_gemm_cur - Cst_lg_gemm[i])), 1e-6*χlarge^1.2, "gemm")
    @test zrtest(mean(abs.(Cst_sm_gemm_cur - Cst_sm_gemm[i])), 1e-6*χsmall^1.2, "gemm")
    @test zrtest(mean(abs.(Cst_lg_hemm_cur - Cst_lg_hemm[i])), 1e-6*χlarge^1.2, "hemm")
    @test zrtest(mean(abs.(Cst_sm_hemm_cur - Cst_sm_hemm[i])), 1e-6*χsmall^1.2, "hemm")
    @test zrtest(mean(abs.(Cst_lg_symm_cur - Cst_lg_symm[i])), 1e-6*χlarge^1.2, "symm")
    @test zrtest(mean(abs.(Cst_sm_symm_cur - Cst_sm_symm[i])), 1e-6*χsmall^1.2, "symm")
end
end

