# TRMM
#

macro blis_linalg_trmm(typename, targetfunc, bliname)

    # Get method for the typed API backend.
    blifuncname = Symbol("bli_", cblas_typechar[ eval(typename) ], bliname)
    blifunc = getproperty(TypedBackend, blifuncname)

    return quote
        """
        BLIS-based TRMM with approximately the same interface as `LinearAlgebra.BLAS`.
        Strided arrays are directly supported without further Julia-level dispatching.
        """
        $(esc(targetfunc))(side::AbstractChar,
                           uplo::AbstractChar,
                           tA::AbstractChar,
                           dA::AbstractChar,
                           alpha::$typename,
                           A::StridedMatrix{$typename},
                           B::StridedMatrix{$typename}) = begin

            bli_sidea = char_to_side[side]
            bli_uploa = char_to_uplo[uplo]
            bli_tA = char_to_trans[tA]
            bli_dA = char_to_diag[dA]

            k , k_= size(A)
            m , n = size(B)

            (k == k_) || throw(DimensionMismatch("TRMM A storage is not square."))
            (bli_sidea != BLIS_LEFT  ||
             k == m ) || throw(DimensionMismatch("TRMM contracted size mismatch."))
            (bli_sidea != BLIS_RIGHT ||
             k == n ) || throw(DimensionMismatch("TRMM contracted size mismatch."))

            $blifunc(bli_sidea,
                     bli_uploa,
                     bli_tA,
                     bli_dA,
                     m, n,
                     [alpha],
                     A, strides(A)...,
                     B, strides(B)...)
            B

        end
    end
end

@blis_linalg_trmm(Float32, trmm!, trmm!)
@blis_linalg_trmm(Float64, trmm!, trmm!)
@blis_linalg_trmm(ComplexF32, trmm!, trmm!)
@blis_linalg_trmm(ComplexF64, trmm!, trmm!)

