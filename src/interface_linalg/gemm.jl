# GEMM
#

macro blis_linalg_gemm(typename, targetfunc, bliname)

    blifunc = Symbol("bli_", cblas_typechar[ eval(typename) ], bliname)

    return quote
        """
        BLIS-based GEMM with approximately the same interface as `LinearAlgebra.BLAS.gemm!`.
        Strided arrays are directly supported without further Julia-level dispatching.
        """
        $(esc(targetfunc))(tA::AbstractChar,
                           tB::AbstractChar,
                           alpha::$typename,
                           A::StridedMatrix{$typename},
                           B::StridedMatrix{$typename},
                           beta::$typename,
                           C::StridedMatrix{$typename}) = begin

            bli_tA = char_to_trans[tA]
            bli_tB = char_to_trans[tB]

            m,  k = size(A)
            k_, n = size(B)
            m_, n_= size(C)

            # Transpose sizes.
            if bli_tA == BLIS_TRANSPOSE || bli_tA == BLIS_CONJ_TRANSPOSE
                m, k = k, m
            end
            if bli_tB == BLIS_TRANSPOSE || bli_tB == BLIS_CONJ_TRANSPOSE
                k_, n = n, k_
            end

            (k == k_) || throw(DimensionMismatch("GEMM contracted size mismatch."))
            (m == m_ && 
             n == n_) || throw(DimensionMismatch("GEMM target buffer size mismatch."))

            $(esc(blifunc))(bli_tA,
                            bli_tB,
                            m, n, k,
                            [alpha],
                            A, strides(A)...,
                            B, strides(B)...,
                            [beta],
                            C, strides(C)...)
            C

        end
    end
end

@blis_linalg_gemm(Float32, gemm!, gemm!)
@blis_linalg_gemm(Float64, gemm!, gemm!)
@blis_linalg_gemm(ComplexF32, gemm!, gemm!)
@blis_linalg_gemm(ComplexF64, gemm!, gemm!)

