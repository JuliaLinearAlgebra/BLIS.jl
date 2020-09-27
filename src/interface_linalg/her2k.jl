# HER2K / SKR2K
#

macro blis_linalg_her2k(typename, btypename, targetfunc, bliname)

    blifunc = Symbol("bli_", cblas_typechar[ eval(typename) ], bliname)

    return quote
        """
        BLIS-based HER2K/SYR2K with approximately the same interface as `LinearAlgebra.BLAS`.
        Strided arrays are directly supported without further Julia-level dispatching.
        """
        $(esc(targetfunc))(uplo::AbstractChar,
                           tAB::AbstractChar,
                           alpha::$typename,
                           A::StridedMatrix{$typename},
                           B::StridedMatrix{$typename},
                           beta::$btypename,
                           C::StridedMatrix{$typename}) = begin

            bli_uploc = char_to_uplo[uplo]
            bli_tAB = char_to_trans[tAB]

            n , k = size(A)
            n_, k_= size(B)
            m , m_= size(C)

            if bli_tAB == BLIS_TRANSPOSE || bli_tAB == BLIS_CONJ_TRANSPOSE
                n , k  = k , n
                n_, k_ = k_, n_
            end

            (m == m_) || throw(DimensionMismatch("HE(SY)R2K C storage not square."))
            (k == k_) || throw(DimensionMismatch("HE(SY)R2K contracted size mismatch."))
            (m == n &&
             n == n_) || throw(DimensionMismatch("HE(SY)R2K external size mismatch."))

            $(esc(blifunc))(bli_uploc,
                            bli_tAB,
                            m, k,
                            [alpha],
                            A, strides(A)...,
                            B, strides(B)...,
                            [beta],
                            C, strides(C)...)
            C

        end
    end
end

@blis_linalg_her2k(Float32, Float32, her2k!, her2k!)
@blis_linalg_her2k(Float64, Float64, her2k!, her2k!)
@blis_linalg_her2k(ComplexF32, Float32, her2k!, her2k!)
@blis_linalg_her2k(ComplexF64, Float64, her2k!, her2k!)

@blis_linalg_her2k(Float32, Float32, syr2k!, syr2k!)
@blis_linalg_her2k(Float64, Float64, syr2k!, syr2k!)
@blis_linalg_her2k(ComplexF32, ComplexF32, syr2k!, syr2k!)
@blis_linalg_her2k(ComplexF64, ComplexF64, syr2k!, syr2k!)

