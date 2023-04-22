# Level-2 LinearAlgebra.BLAS interface.
#
# TODO: support gbmv!, sbmv! and hbmv!.
#       hpmv! and spmv! is not possible within current frame.

# Level2 BLAS size checker.
bli_check_lv2(trans::BliTrans,
              m ::Integer, n ::Integer,
              n_::Integer, m_::Integer) = begin

    ((trans.enum & BLIS_TRANS_BIT) != 0) && ((m, n) = (n, m))
    (m == m_) || throw(DimensionMismatch("Target buffer size mismatch."))
    (n == n_) || throw(DimensionMismatch("Contracted size mismatch."))

    nothing
end

macro blis_interface_linalg_lv2_gemv(Tc1, T1, T2, Tc2, T3, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        $(esc(targetfunc))(tA::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           x::StridedVector{<:$T2},
                           β::$Tc2,
                           y::StridedVector{<:$T3}) = begin

            bli_tA = char_to_trans[tA]

            bli_check_lv2(bli_tA,
                          size(A)...,
                          length(x),
                          length(y))

            oα = BliObj(α)
            oA = BliObj(A)
            ox = BliObj(x)
            oβ = BliObj(β)
            oy = BliObj(y)

            ObjectBackend.bli_obj_set_onlytrans!(bli_tA, oA.obj)
            $blifunc(oα, oA, ox, oβ, oy)
            y

        end
    end
end

@doc """
        gemv!(tA, α, A, x, β, y)
    BLIS-based GEMV with strides support.
    """ gemv!
@blis_interface_linalg_lv2_gemv Float32    Float32    Float32    Float32    Float32    gemv! gemv!
@blis_interface_linalg_lv2_gemv Float64    Float64    Float64    Float64    Float64    gemv! gemv!
@blis_interface_linalg_lv2_gemv ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 gemv! gemv!
@blis_interface_linalg_lv2_gemv ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 gemv! gemv!

macro blis_interface_linalg_lv2_hemv(Tc1, T1, T2, Tc2, T3, targetfunc, bliname, bli_struc)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        $(esc(targetfunc))(ul::AbstractChar,
                           α::$Tc1,
                           A::StridedMatrix{<:$T1},
                           x::StridedVector{<:$T2},
                           β::$Tc2,
                           y::StridedVector{<:$T3}) = begin

            bli_ul = char_to_uplo[ul]

            bli_check_lv2(BLIS_NO_TRANSPOSE,
                          size(A)...,
                          length(x),
                          length(y))

            oα = BliObj(α)
            oA = BliObj(A)
            ox = BliObj(x)
            oβ = BliObj(β)
            oy = BliObj(y)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oA.obj)
            ObjectBackend.bli_obj_set_struc!($bli_struc, oA.obj)
            $blifunc(oα, oA, ox, oβ, oy)
            y

        end
    end
end

@doc """
        hemv!(ul, α, A, x, β, y)
    BLIS-based HEMV with strides support.
    """ hemv!
@blis_interface_linalg_lv2_hemv Float32    Float32    Float32    Float32    Float32    hemv! hemv! BLIS_HERMITIAN
@blis_interface_linalg_lv2_hemv Float64    Float64    Float64    Float64    Float64    hemv! hemv! BLIS_HERMITIAN
@blis_interface_linalg_lv2_hemv ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 hemv! hemv! BLIS_HERMITIAN
@blis_interface_linalg_lv2_hemv ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 hemv! hemv! BLIS_HERMITIAN

@doc """
        symv!(ul, α, A, x, β, y)
    BLIS-based SYMV with strides support.
    """ symv!
@blis_interface_linalg_lv2_hemv Float32    Float32    Float32    Float32    Float32    symv! symv! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_hemv Float64    Float64    Float64    Float64    Float64    symv! symv! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_hemv ComplexF32 ComplexF32 ComplexF32 ComplexF32 ComplexF32 symv! symv! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_hemv ComplexF64 ComplexF64 ComplexF64 ComplexF64 ComplexF64 symv! symv! BLIS_SYMMETRIC

macro blis_interface_linalg_lv2_trmv(T1, T2, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        $(esc(targetfunc))(ul::AbstractChar,
                           tA::AbstractChar,
                           dA::AbstractChar,
                           A::StridedMatrix{<:$T1},
                           b::StridedVector{<:$T2}) = begin

            bli_ul = char_to_uplo[ul]
            bli_dA = char_to_diag[dA]
            bli_tA = char_to_trans[tA]

            bli_check_lv2(bli_tA,
                          size(A)...,
                          length(b),
                          length(b))

            oα = BliObj($T1(1.0))
            oA = BliObj(A)
            ob = BliObj(b)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oA.obj)
            ObjectBackend.bli_obj_set_diag!(bli_dA, oA.obj)
            ObjectBackend.bli_obj_set_onlytrans!(bli_tA, oA.obj)
            ObjectBackend.bli_obj_set_struc!(BLIS_TRIANGULAR, oA.obj)
            $blifunc(oα, oA, ob)
            b

        end
    end
end

@doc """
        trmv!(ul, tA, dA, A, b)
    BLIS-based TRMV with strides support.
    """
@blis_interface_linalg_lv2_trmv Float32    Float32    trmv! trmv!
@blis_interface_linalg_lv2_trmv Float64    Float64    trmv! trmv!
@blis_interface_linalg_lv2_trmv ComplexF32 ComplexF32 trmv! trmv!
@blis_interface_linalg_lv2_trmv ComplexF64 ComplexF64 trmv! trmv!

@doc """
        trsv!(ul, tA, dA, A, b)
    BLIS-based TRSV with strides support.
    """
@blis_interface_linalg_lv2_trmv Float32    Float32    trsv! trsv!
@blis_interface_linalg_lv2_trmv Float64    Float64    trsv! trsv!
@blis_interface_linalg_lv2_trmv ComplexF32 ComplexF32 trsv! trsv!
@blis_interface_linalg_lv2_trmv ComplexF64 ComplexF64 trsv! trsv!

macro blis_interface_linalg_lv2_ger(Tc1, T1, T2, T3, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        $(esc(targetfunc))(α::$Tc1,
                           x::StridedVector{<:$T1},
                           y::StridedVector{<:$T2},
                           A::StridedMatrix{<:$T3}) = begin

            bli_check_lv2(BLIS_NO_TRANSPOSE,
                          size(A)...,
                          length(y),
                          length(x))

            oα = BliObj(α)
            ox = BliObj(x)
            oy = BliObj(y)
            oA = BliObj(A)

            $blifunc(oα, ox, oy, oA)
            A

        end
    end
end

@doc """
        ger!(α, x, y, A)
    BLIS-based GER with strides support.
    """ ger!
@blis_interface_linalg_lv2_ger Float32    Float32    Float32    Float32    ger! ger!
@blis_interface_linalg_lv2_ger Float64    Float64    Float64    Float64    ger! ger!
@blis_interface_linalg_lv2_ger ComplexF32 ComplexF32 ComplexF32 ComplexF32 ger! ger!
@blis_interface_linalg_lv2_ger ComplexF64 ComplexF64 ComplexF64 ComplexF64 ger! ger!

macro blis_interface_linalg_lv2_her(Tc1, T1, T2, targetfunc, bliname, bli_struc)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        $(esc(targetfunc))(ul::AbstractChar,
                           α::$Tc1,
                           x::StridedVector{<:$T1},
                           A::StridedMatrix{<:$T2}) = begin
            
            bli_ul = char_to_uplo[ul]

            bli_check_lv2(BLIS_NO_TRANSPOSE,
                          size(A)...,
                          length(x),
                          length(x))

            oα = BliObj(α)
            ox = BliObj(x)
            oA = BliObj(A)

            ObjectBackend.bli_obj_set_uplo!(bli_ul, oA.obj)
            ObjectBackend.bli_obj_set_struc!($bli_struc, oA.obj)
            $blifunc(oα, ox, oA)
            A

        end
    end
end

@doc """
        her!(uplo, α, x, A)
    BLIS-based HER with strides support.
    """ her!
@blis_interface_linalg_lv2_her Float32    Float32    Float32    her! her! BLIS_HERMITIAN
@blis_interface_linalg_lv2_her Float64    Float64    Float64    her! her! BLIS_HERMITIAN
@blis_interface_linalg_lv2_her Float32    ComplexF32 ComplexF32 her! her! BLIS_HERMITIAN
@blis_interface_linalg_lv2_her Float64    ComplexF64 ComplexF64 her! her! BLIS_HERMITIAN

@doc """
        syr!(uplo, α, x, A)
    BLIS-based SYR with strides support.
    """ syr!
@blis_interface_linalg_lv2_her Float32    Float32    Float32    syr! syr! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_her Float64    Float64    Float64    syr! syr! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_her ComplexF32 ComplexF32 ComplexF32 syr! syr! BLIS_SYMMETRIC
@blis_interface_linalg_lv2_her ComplexF64 ComplexF64 ComplexF64 syr! syr! BLIS_SYMMETRIC

