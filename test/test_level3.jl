# Level-3 BLAS: Test & evict.
#

macro l3_test_evict_(func, χ, αtype, βtype, atype, btype, ctype, largs)
    return quote
        local χ  = $χ
        local αr = 1.1
        local βr = 1.1
        local αc = 1.1 + 0.2im
        local βc = 1.1 + 0.3im
        αr = $αtype(αr)
        βr = $βtype(βr)
        if $αtype <: Complex
            αc = $αtype(αc)
        end
        if $βtype <: Complex
            βc = $βtype(βc)
        end

        local A = rand($atype, χ, χ)
        local B = rand($btype, χ, χ)
        local C = rand($ctype, χ, χ)
        local C_= copy(C)

        @run_evict $func C $(largs)

        A = $ctype.(A)
        B = $ctype.(B)

        # Check.
        $func($(largs)..., C_)
        abs.(C - C_)
    end
end

macro l3_test_evict(func, χ, atype, btype, ctype, largs)
    return quote
        @l3_test_evict_ $func $χ $ctype $ctype $atype $btype $ctype $largs
    end
end

@testset "BLAS level-3 LinearAlgebra interface" begin

    @test zrtest(reduce(max, @l3_test_evict BLAS.gemm!  100    Float32 ComplexF32 ComplexF64 ('N', 'N', αr, A, B, βr) ), 1e-6, "xgemm_nn_100")
    @test zrtest(reduce(max, @l3_test_evict BLAS.gemm!   20    Float32    Float32    Float32 ('N', 'N', αr, A, B, βr) ), 1e-3, "sgemm_nn_20")
    @test zrtest(reduce(max, @l3_test_evict BLAS.gemm! 1200    Float64    Float64    Float64 ('T', 'N', αr, A, B, βr) ), 1e-6, "dgemm_tn_1200")
    @test zrtest(reduce(max, @l3_test_evict BLAS.gemm!  240 ComplexF32 ComplexF32 ComplexF32 ('T', 'T', αc, A, B, βc) ), 1e-3, "cgemm_tt_240")
    @test zrtest(reduce(max, @l3_test_evict BLAS.gemm! 1200 ComplexF64 ComplexF64 ComplexF64 ('N', 'T', αc, A, B, βc) ), 1e-6, "zgemm_nt_1200")

    @test zrtest(reduce(max, @l3_test_evict BLAS.hemm!  240 ComplexF32 ComplexF32 ComplexF32 ('R', 'L', αc, A, B, βc) ), 1e-3, "chemm_rl_240")
    @test zrtest(reduce(max, @l3_test_evict BLAS.hemm! 1200 ComplexF64 ComplexF64 ComplexF64 ('R', 'L', αc, A, B, βc) ), 1e-6, "zhemm_rl_1200")

    @test zrtest(reduce(max, @l3_test_evict BLAS.symm! 2000    Float32    Float32    Float32 ('L', 'U', αr, A, B, βr) ), 1e-3, "ssymm_lu_2000")
    @test zrtest(reduce(max, @l3_test_evict BLAS.symm!   32    Float64    Float64    Float64 ('R', 'U', αr, A, B, βr) ), 1e-6, "dsymm_ru_32")
    @test zrtest(reduce(max, @l3_test_evict BLAS.symm! 2400 ComplexF32 ComplexF32 ComplexF32 ('R', 'L', αc, A, B, βc) ), 1e-2, "csymm_rl_2400")
    @test zrtest(reduce(max, @l3_test_evict BLAS.symm!   50 ComplexF64 ComplexF64 ComplexF64 ('R', 'L', αc, A, B, βc) ), 1e-6, "zsymm_rl_50")

    @test zrtest(reduce(max, @l3_test_evict_ BLAS.her2k!  240 ComplexF32 Float32 ComplexF32 ComplexF32 ComplexF32 ('U', 'N', αr, A, B, βr) ), 1e-3, "cher2k_un_240")
    @test zrtest(reduce(max, @l3_test_evict_ BLAS.her2k! 1200 ComplexF64 Float64 ComplexF64 ComplexF64 ComplexF64 ('L', 'N', αr, A, B, βr) ), 1e-6, "zher2k_lt_1200")

    @test zrtest(reduce(max, @l3_test_evict BLAS.syr2k! 2000    Float32    Float32    Float32 ('U', 'N', αr, A, B, βr) ), 1e-3, "ssyr2k_un_2000")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syr2k!   32    Float64    Float64    Float64 ('U', 'T', αr, A, B, βr) ), 1e-6, "dsyr2k_ut_32")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syr2k! 2400 ComplexF32 ComplexF32 ComplexF32 ('L', 'N', αc, A, B, βc) ), 1e-2, "csyr2k_ln_2400")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syr2k!   50 ComplexF64 ComplexF64 ComplexF64 ('L', 'T', αc, A, B, βc) ), 1e-6, "zsyr2k_lt_50")

    @test zrtest(reduce(max, @l3_test_evict_ BLAS.herk!  240 Float32 Float32 ComplexF32 ComplexF32 ComplexF32 ('U', 'N', αr, A, βr) ), 1e-3, "cherk_un_240")
    @test zrtest(reduce(max, @l3_test_evict_ BLAS.herk! 1200 Float64 Float64 ComplexF64 ComplexF64 ComplexF64 ('L', 'N', αr, A, βr) ), 1e-6, "zherk_lt_1200")

    @test zrtest(reduce(max, @l3_test_evict BLAS.syrk! 2000    Float32    Float32    Float32 ('U', 'N', αr, A, βr) ), 1e-3, "ssyrk_un_2000")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syrk!   32    Float64    Float64    Float64 ('U', 'T', αr, A, βr) ), 1e-6, "dsyrk_ut_32")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syrk! 2400 ComplexF32 ComplexF32 ComplexF32 ('L', 'N', αc, A, βc) ), 1e-2, "csyrk_ln_2400")
    @test zrtest(reduce(max, @l3_test_evict BLAS.syrk!   50 ComplexF64 ComplexF64 ComplexF64 ('L', 'T', αc, A, βc) ), 1e-6, "zsyrk_lt_50")

    @test zrtest(reduce(max, @l3_test_evict BLAS.trmm!   20    Float32    Float32    Float32 ('L', 'U', 'N', 'U', αr, A) ), 1e-3, "strmm_lunu_20")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trmm! 1200    Float64    Float64    Float64 ('L', 'L', 'N', 'U', αr, A) ), 1e-6, "dtrmm_llnu_1200")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trmm!  240 ComplexF32 ComplexF32 ComplexF32 ('R', 'U', 'T', 'N', αc, A) ), 1e-3, "ctrmm_rutn_240")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trmm! 1200 ComplexF64 ComplexF64 ComplexF64 ('R', 'L', 'T', 'N', αc, A) ), 1e-6, "ztrmm_rltn_1200")

    @test zrtest(reduce(max, @l3_test_evict BLAS.trsm!  100    Float32    Float32    Float32 ('L', 'U', 'N', 'U', αr, A) ), 1e-3, "strsm_un_100")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trsm!  100    Float64    Float64    Float64 ('L', 'L', 'N', 'U', αr, A) ), 1e-4, "dtrsm_ut_100")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trsm!   30 ComplexF32 ComplexF32 ComplexF32 ('R', 'U', 'T', 'N', αc, A) ), 1e-2, "ctrsm_ln_30")
    @test zrtest(reduce(max, @l3_test_evict BLAS.trsm!   30 ComplexF64 ComplexF64 ComplexF64 ('R', 'L', 'T', 'N', αc, A) ), 1e-4, "ztrsm_lt_30")
end

