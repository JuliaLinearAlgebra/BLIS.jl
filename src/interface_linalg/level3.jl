# Level-3 LinearAlgebra.BLAS interface.
#

" Level3 BLAS size checker. "
bli_check_lv3(transa::BliTrans,
              transb::BliTrans,
              transc::BliTrans,
              m ::Integer, k ::Integer,
              k_::Integer, n ::Integer,
              m_::Integer, n_::Integer) = begin
    # Assumption:
    # m, k * k_, n -> m_, n_

    ((transa.enum & BLIS_TRANS_BIT) != 0) && ((m , k ) = (k , m ))
    ((transb.enum & BLIS_TRANS_BIT) != 0) && ((k_, n ) = (n , k_))
    ((transc.enum & BLIS_TRANS_BIT) != 0) && ((m_, n_) = (n_, m_))

    (m == m_ && 
     n == n_) || throw(DimensionMismatch("Target buffer size mismatch."))
    (k == k_) || throw(DimensionMismatch("Contracted size mismatch."))

    nothing
end

macro blis_interface_linalg_lv3_gemm(Tc1, T1, T2, Tc2, T3, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            gemm!(tA, tB, α, A, B, β, C)
        BLIS-based GEMM with generic strides & mixed precision directly supported.
        """
        $(esc(targetfunc))(tA::AbstractChar,
                           tB::AbstractChar,
                           α::$Tc1,
                           A::StridedVecOrMat{<:$T1},
                           B::StridedVecOrMat{<:$T2},
                           β::$Tc2,
                           C::StridedVecOrMat{<:$T3}) = begin

            bli_tA = char_to_trans[tA]
            bli_tB = char_to_trans[tB]

            bli_check_lv3(bli_tA, bli_tB, BLIS_NO_TRANSPOSE,
                          size(A)...,
                          size(B)...,
                          size(C)...)

            oα = BliObj(α)
            oA = BliObj(A)
            oB = BliObj(B)
            oβ = BliObj(β)
            oC = BliObj(C)

            ObjectBackend.bli_obj_set_onlytrans!(bli_tA, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tB, oB.obj)
            $blifunc(oα, oA, oB, oβ, oC)
            C

        end
    end
end

@blis_interface_linalg_lv3_gemm(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                gemm!, gemm!)
# To avoid obfuscating with LinearAlgebra.gemm!,
#  instantiate more specifically.
@blis_interface_linalg_lv3_gemm Float32    Float32    Float32    Float32    Float32    gemm! gemm!
@blis_interface_linalg_lv3_gemm Float64    Float64    Float64    Float64    Float64    gemm! gemm!
@blis_interface_linalg_lv3_gemm ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 gemm! gemm!
@blis_interface_linalg_lv3_gemm ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 gemm! gemm!

# Do this for gemm_wrapper! as well.
# Direct overriding triggers warning.
@blis_interface_linalg_lv3_gemm_wrapper Float32    Float32    Float32
@blis_interface_linalg_lv3_gemm_wrapper Float64    Float64    Float64
@blis_interface_linalg_lv3_gemm_wrapper ComplexF32 ComplexF32 ComplexF32
@blis_interface_linalg_lv3_gemm_wrapper ComplexF64 ComplexF64 ComplexF64

macro blis_interface_linalg_lv3_hemm(Tc1, T1, T2, Tc2, T3, targetfunc, bliname, bli_struc)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            ??mm(side, uplo, α, A, B, β, C)
        BLIS-based HEMM/SYMM with generic strides & mixed precision directly supported.
        """
        $(esc(targetfunc))(si::AbstractChar,
                           ul::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           B::StridedMatrix{<:$T2},
                           β::$Tc2,
                           C::StridedMatrix{<:$T3}) = begin

            bli_si = char_to_side[si]
            bli_ul = char_to_uplo[ul]

            if bli_si == BLIS_LEFT
                bli_check_lv3(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
                              size(A)...,
                              size(B)...,
                              size(C)...)
            else
                bli_check_lv3(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
                              size(B)...,
                              size(A)...,
                              size(C)...)
            end

            oα = BliObj(α)
            oA = BliObj(A)
            oB = BliObj(B)
            oβ = BliObj(β)
            oC = BliObj(C)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oA.obj)
            ObjectBackend.bli_obj_set_struc!($bli_struc, oA.obj)
            $blifunc(bli_si, oα, oA, oB, oβ, oC)
            C

        end
    end
end

@blis_interface_linalg_lv3_hemm(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                hemm!, hemm!,
                                BLIS_HERMITIAN)
@blis_interface_linalg_lv3_hemm Float32    Float32    Float32    Float32    Float32    hemm! hemm! BLIS_HERMITIAN
@blis_interface_linalg_lv3_hemm Float64    Float64    Float64    Float64    Float64    hemm! hemm! BLIS_HERMITIAN
@blis_interface_linalg_lv3_hemm ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 hemm! hemm! BLIS_HERMITIAN
@blis_interface_linalg_lv3_hemm ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 hemm! hemm! BLIS_HERMITIAN

@blis_interface_linalg_lv3_hemm(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                symm!, symm!,
                                BLIS_SYMMETRIC)
@blis_interface_linalg_lv3_hemm Float32    Float32    Float32    Float32    Float32    symm! symm! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_hemm Float64    Float64    Float64    Float64    Float64    symm! symm! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_hemm ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 symm! symm! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_hemm ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 symm! symm! BLIS_SYMMETRIC

macro blis_interface_linalg_lv3_her2k(Tc1, T1, T2, Tc2, T3, targetfunc, bliname, bli_struc)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            ??r2k(uplo, tAB, α, A, B, β, C)
        BLIS-based HER2K/SYR2K with generic strides & mixed precision directly supported.
        """
        $(esc(targetfunc))(ul::AbstractChar,
                           tAB::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           B::StridedMatrix{<:$T2},
                           β::$Tc2,
                           C::StridedMatrix{<:$T3}) = begin

            bli_ul = char_to_uplo[ul]
            bli_tAB = char_to_side[tAB]

            bli_check_lv3(bli_tAB, bli_tAB, BLIS_NO_TRANSPOSE,
                          size(A)...,
                          size(B')...,
                          size(C)...)

            oα = BliObj(α)
            oA = BliObj(A)
            oB = BliObj(B)
            oβ = BliObj(β)
            oC = BliObj(C)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oC.obj)
            ObjectBackend.bli_obj_set_struc!($bli_struc, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tAB, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tAB, oB.obj)
            $blifunc(oα, oA, oB, oβ, oC)
            C

        end
    end
end

@blis_interface_linalg_lv3_her2k(BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 her2k!, her2k!,
                                 BLIS_HERMITIAN)
@blis_interface_linalg_lv3_her2k Float32    Float32    Float32    Float32    Float32    her2k! her2k! BLIS_HERMITIAN
@blis_interface_linalg_lv3_her2k Float64    Float64    Float64    Float64    Float64    her2k! her2k! BLIS_HERMITIAN
@blis_interface_linalg_lv3_her2k ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 her2k! her2k! BLIS_HERMITIAN
@blis_interface_linalg_lv3_her2k ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 her2k! her2k! BLIS_HERMITIAN

@blis_interface_linalg_lv3_her2k(BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 BliCompatibleType,
                                 syr2k!, syr2k!,
                                 BLIS_SYMMETRIC)
@blis_interface_linalg_lv3_her2k Float32    Float32    Float32    Float32    Float32    syr2k! syr2k! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_her2k Float64    Float64    Float64    Float64    Float64    syr2k! syr2k! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_her2k ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 syr2k! syr2k! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_her2k ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 syr2k! syr2k! BLIS_SYMMETRIC

macro blis_interface_linalg_lv3_herk(Tc1, T1, Tc2, T2, targetfunc, bliname, bli_struc)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            ??rk(uplo, tA, α, A, β, C)
        BLIS-based HER2K/SYR2K with generic strides & mixed precision directly supported.
        """
        $(esc(targetfunc))(ul::AbstractChar,
                           tA::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           β::$Tc2,
                           C::StridedMatrix{<:$T2}) = begin

            bli_ul = char_to_uplo[ul]
            bli_tA = char_to_side[tA]

            bli_check_lv3(bli_tA, bli_tA, BLIS_NO_TRANSPOSE,
                          size(A)...,
                          size(A')...,
                          size(C)...)

            oα = BliObj(α)
            oA = BliObj(A)
            oβ = BliObj(β)
            oC = BliObj(C)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oC.obj)
            ObjectBackend.bli_obj_set_struc!($bli_struc, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tA, oA.obj)
            $blifunc(oα, oA, oβ, oC)
            C

        end
    end
end

@blis_interface_linalg_lv3_herk(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                herk!, herk!,
                                BLIS_HERMITIAN)
@blis_interface_linalg_lv3_herk Float32    Float32    Float32    Float32    herk! herk! BLIS_HERMITIAN
@blis_interface_linalg_lv3_herk Float64    Float64    Float64    Float64    herk! herk! BLIS_HERMITIAN
@blis_interface_linalg_lv3_herk ComplexF32 ComplexF32 ComplexF32 ComplexF32 herk! herk! BLIS_HERMITIAN
@blis_interface_linalg_lv3_herk ComplexF64 ComplexF64 ComplexF64 ComplexF64 herk! herk! BLIS_HERMITIAN

@blis_interface_linalg_lv3_herk(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                syrk!, syrk!,
                                BLIS_SYMMETRIC)
@blis_interface_linalg_lv3_herk Float32    Float32    Float32    Float32    syrk! syrk! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_herk Float64    Float64    Float64    Float64    syrk! syrk! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_herk ComplexF32 ComplexF32 ComplexF32 ComplexF32 syrk! syrk! BLIS_SYMMETRIC
@blis_interface_linalg_lv3_herk ComplexF64 ComplexF64 ComplexF64 ComplexF64 syrk! syrk! BLIS_SYMMETRIC

macro blis_interface_linalg_lv3_trmm(Tc1, T1, T2, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            tr?m(side, uplo, tA, dA, α, A, B)
        BLIS-based TRMM/TRSM with generic strides & mixed precision directly supported.
        """
        $(esc(targetfunc))(si::AbstractChar,
                           ul::AbstractChar,
                           tA::AbstractChar,
                           dA::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           B::StridedMatrix{<:$T2}) = begin

            bli_si = char_to_side[si]
            bli_ul = char_to_uplo[ul]
            bli_dA = char_to_diag[dA]
            bli_tA = char_to_trans[tA]

            if bli_si == BLIS_LEFT
                bli_check_lv3(bli_tA, BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
                              size(A)...,
                              size(B)...,
                              size(B)...)
            else
                bli_check_lv3(BLIS_NO_TRANSPOSE, bli_tA, BLIS_NO_TRANSPOSE,
                              size(B)...,
                              size(A)...,
                              size(B)...)
            end

            oα = BliObj(α)
            oA = BliObj(A)
            oB = BliObj(B)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oA.obj)
            ObjectBackend.bli_obj_set_diag!(bli_dA, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tA, oA.obj)
            $blifunc(bli_si, oα, oA, oB)
            B

        end
    end
end

@blis_interface_linalg_lv3_trmm(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                trmm!, trmm!)
@blis_interface_linalg_lv3_trmm Float32    Float32    Float32    trmm! trmm!
@blis_interface_linalg_lv3_trmm Float64    Float64    Float64    trmm! trmm!
@blis_interface_linalg_lv3_trmm ComplexF32 ComplexF32 ComplexF32 trmm! trmm!
@blis_interface_linalg_lv3_trmm ComplexF64 ComplexF64 ComplexF64 trmm! trmm!

@blis_interface_linalg_lv3_trmm(BliCompatibleType,
                                BliCompatibleType,
                                BliCompatibleType,
                                trsm!, trsm!)
@blis_interface_linalg_lv3_trmm Float32    Float32    Float32    trsm! trsm!
@blis_interface_linalg_lv3_trmm Float64    Float64    Float64    trsm! trsm!
@blis_interface_linalg_lv3_trmm ComplexF32 ComplexF32 ComplexF32 trsm! trsm!
@blis_interface_linalg_lv3_trmm ComplexF64 ComplexF64 ComplexF64 trsm! trsm!

