# Level-1 BLAS: Test & evict.
#

macro l1_test_evict(func, χ, type, ctype, largs)
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

        local x = rand($type, χ)
        local y = rand($type, χ)
        local y_= copy(y)

        @run_evict $func y $(largs)

        x = $type.(x)

        # Check.
        $func($(largs)..., y_)
        abs.(y - y_)
    end
end


@testset "BLAS level-1 LinearAlgebra interface" begin

    @test zrtest(reduce(max, @l1_test_evict BLAS.axpy!  100    Float32    Float64 (αr, x) ), 1e-6, "xaxpy_100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpy!  100    Float32    Float32 (αr, x) ), 1e-6, "saxpy_100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpy!   30    Float64    Float64 (αr, x) ), 1e-6, "daxpy_30")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpy! 1100 ComplexF32 ComplexF32 (αc, x) ), 4e-3, "caxpy_1100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpy!  300 ComplexF64 ComplexF64 (αc, x) ), 4e-7, "zaxpy_300")

    @test zrtest(reduce(max, @l1_test_evict BLAS.axpby!  100    Float32    Float64 (αr, x, βr) ), 1e-6, "xaxpby_100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpby!  100    Float32    Float32 (αr, x, βr) ), 1e-6, "saxpby_100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpby!   30    Float64    Float64 (αr, x, βr) ), 1e-6, "daxpby_30")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpby! 1100 ComplexF32 ComplexF32 (αc, x, βc) ), 4e-3, "caxpby_1100")
    @test zrtest(reduce(max, @l1_test_evict BLAS.axpby!  300 ComplexF64 ComplexF64 (αc, x, βc) ), 4e-7, "zaxpby_300")

end
