# Level-2 BLAS: Test & evict.
#

macro l2_test_evict(func, χ, type, largs)
    return quote
        local χ  = $χ
        local αr = 1.1
        local βr = 1.1
        local αc = 1.1 + 0.2im
        local βc = 1.1 + 0.3im
        αr = $type(αr)
        βr = $type(βr)
        if $type <: Complex
            αc = $type(αc)
        end
        if $type <: Complex
            βc = $type(βc)
        end

        local A = rand($type, χ, χ)
        local x = rand($type, χ)
        local y = rand($type, χ)
        local y_= copy(y)

        @run_evict $func y $(largs)

        A = $type.(A)
        x = $type.(x)

        # Check.
        $func($(largs)..., y_)
        abs.(y - y_)
    end
end

macro l2r_test_evict(func, χ, αtype, type, largs)
    return quote
        local χ  = $χ
        local αr = 1.1
        local βr = 1.1
        local αc = 1.1 + 0.2im
        local βc = 1.1 + 0.3im
        αr = $αtype(αr)
        βr =  $type(βr)
        if $type <: Complex
            αc = $type(αc)
        end
        if $type <: Complex
            βc = $type(βc)
        end

        local x = rand($type, χ)
        local y = rand($type, χ)
        local C = rand($type, χ, χ)
        local C_= copy(C)

        @run_evict $func C $(largs)

        x = $type.(x)
        y = $type.(y)

        # Check.
        $func($(largs)..., C_)
        abs.(C - C_)
    end
end

@testset "BLAS level-2 LinearAlgebra interface" begin

    @test zrtest(reduce(max, @l2_test_evict BLAS.gemv!  100    Float32 ('N', αr, A, x, βr) ), 1e-4, "sgemv_n_100")
    @test zrtest(reduce(max, @l2_test_evict BLAS.gemv!   30    Float64 ('T', αr, A, x, βr) ), 1e-6, "dgemv_t_30")
    @test zrtest(reduce(max, @l2_test_evict BLAS.gemv! 1100 ComplexF32 ('N', αc, A, x, βc) ), 4e-3, "cgemv_n_1100")
    @test zrtest(reduce(max, @l2_test_evict BLAS.gemv!  300 ComplexF64 ('N', αc, A, x, βc) ), 4e-7, "zgemv_n_300")

    @test zrtest(reduce(max, @l2_test_evict BLAS.hemv! 1100 ComplexF32 ('U', αc, A, x, βc) ), 4e-3, "chemv_u_1100")
    @test zrtest(reduce(max, @l2_test_evict BLAS.hemv!  300 ComplexF64 ('L', αc, A, x, βc) ), 4e-7, "zhemv_l_300")

    @test zrtest(reduce(max, @l2_test_evict BLAS.symv!   32    Float32 ('U', αr, A, x, βr) ), 1e-4, "ssymv_u_32")
    @test zrtest(reduce(max, @l2_test_evict BLAS.symv! 2000    Float64 ('U', αr, A, x, βr) ), 1e-3, "dsymv_u_2000")
    @test zrtest(reduce(max, @l2_test_evict BLAS.symv!   24 ComplexF32 ('L', αc, A, x, βc) ), 4e-6, "csymv_l_24")
    @test zrtest(reduce(max, @l2_test_evict BLAS.symv! 1200 ComplexF64 ('L', αc, A, x, βc) ), 4e-7, "zsymv_l_1200")

    @test zrtest(reduce(max, @l2_test_evict BLAS.trmv!   32    Float32 ('U', 'N', 'U', A) ), 1e-4, "strmv_u_32")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trmv! 2000    Float64 ('U', 'T', 'N', A) ), 1e-3, "dtrmv_u_2000")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trmv!   24 ComplexF32 ('L', 'N', 'U', A) ), 4e-3, "ctrmv_l_24")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trmv! 1200 ComplexF64 ('L', 'T', 'N', A) ), 4e-7, "ztrmv_l_1200")

    @test zrtest(reduce(max, @l2_test_evict BLAS.trsv!  100    Float32 ('U', 'N', 'U', A) ), 1e-2, "strsv_u_100")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trsv!   30    Float64 ('U', 'T', 'N', A) ), 1e-4, "dtrsv_u_30")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trsv!  100 ComplexF32 ('L', 'N', 'U', A) ), 1e-1, "ctrsv_l_100")
    @test zrtest(reduce(max, @l2_test_evict BLAS.trsv!  100 ComplexF64 ('L', 'T', 'N', A) ), 1e-2, "ztrsv_l_100")

    @test zrtest(reduce(max, @l2r_test_evict BLAS.ger!  100    Float32    Float32 (αr, x, y) ), 1e-4, "sger_100")
    @test zrtest(reduce(max, @l2r_test_evict BLAS.ger!   30    Float64    Float64 (αr, x, y) ), 1e-6, "dger_30")
    # @test zrtest(reduce(max, @l2r_test_evict BLAS.ger! 1100 ComplexF32 ComplexF32 (αc, x, y) ), 4e-3, "cger_1100")
    # @test zrtest(reduce(max, @l2r_test_evict BLAS.ger!  300 ComplexF64 ComplexF64 (αc, x, y) ), 4e-7, "zger_300")

    @test zrtest(reduce(max, @l2r_test_evict BLAS.her!  200    Float32 ComplexF32 ('U', αr, x) ), 4e-4, "csyr_u_200")
    @test zrtest(reduce(max, @l2r_test_evict BLAS.her!  800    Float64 ComplexF64 ('L', αr, x) ), 4e-5, "zsyr_l_800")

    @test zrtest(reduce(max, @l2r_test_evict BLAS.syr!  100    Float32    Float32 ('U', αr, x) ), 1e-4, "ssyr_u_100")
    @test zrtest(reduce(max, @l2r_test_evict BLAS.syr!   30    Float64    Float64 ('L', αr, x) ), 1e-6, "dsyr_l_30")
    @test zrtest(reduce(max, @l2r_test_evict BLAS.syr! 1100 ComplexF32 ComplexF32 ('U', αc, x) ), 4e-3, "csyr_u_1100")
    @test zrtest(reduce(max, @l2r_test_evict BLAS.syr!  300 ComplexF64 ComplexF64 ('L', αc, x) ), 4e-7, "zsyr_l_300")

end
