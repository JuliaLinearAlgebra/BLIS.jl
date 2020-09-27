# HEMM / SYMM
#

macro blis_linalg_hemm(typename, targetfunc, bliname)

    blifunc = Symbol("bli_", cblas_typechar[ eval(typename) ], bliname)

    return quote
        """
        BLIS-based HEMM/SYMM with approximately the same interface as `LinearAlgebra.BLAS`.
        Strided arrays are directly supported without further Julia-level dispatching.
        """
        $(esc(targetfunc))(side::AbstractChar,
                           uplo::AbstractChar,
                           alpha::$typename,
                           A::StridedMatrix{$typename},
                           B::StridedMatrix{$typename},
                           beta::$typename,
                           C::StridedMatrix{$typename}) = begin

            bli_sidea = char_to_side[side]
            bli_uploa = char_to_uplo[uplo]

            k , k_= size(A)
            m , n = size(B)
            m_, n_= size(C)

            (k == k_) || throw(DimensionMismatch("HE(SY)MM A storage is not square."))
            (m == m_ && 
             n == n_) || throw(DimensionMismatch("HE(SY)MM target buffer size mismatch."))
            (bli_sidea != BLIS_LEFT  || 
             k == m ) || throw(DimensionMismatch("HE(SY)MM contracted size mismatch."))
            (bli_sidea != BLIS_RIGHT || 
             k == n ) || throw(DimensionMismatch("HE(SY)MM contracted size mismatch."))

            $(esc(blifunc))(bli_sidea,
                            bli_uploa,
                            BLIS_NO_CONJUGATE,
                            BLIS_NO_TRANSPOSE,
                            m, n,
                            [alpha],
                            A, strides(A)...,
                            B, strides(B)...,
                            [beta],
                            C, strides(C)...)
            C

        end
    end
end

@blis_linalg_hemm(Float32, hemm!, hemm!)
@blis_linalg_hemm(Float64, hemm!, hemm!)
@blis_linalg_hemm(ComplexF32, hemm!, hemm!)
@blis_linalg_hemm(ComplexF64, hemm!, hemm!)

@blis_linalg_hemm(Float32, symm!, symm!)
@blis_linalg_hemm(Float64, symm!, symm!)
@blis_linalg_hemm(ComplexF32, symm!, symm!)
@blis_linalg_hemm(ComplexF64, symm!, symm!)

