# Level-1 LinearAlgebra.BLAS interface.
#

# NOTE: scal! and blascopy! have incx in its parmeter.
#       No further StridedVector interface is provided.
macro blis_interface_linalg_lv1_scal(T1, targetfunc, bliname)

    # Get method for the typed API backend.
    blifuncname = Symbol("bli_", cblas_typechar[ eval(T1) ], bliname)
    blifunc = getproperty(TypedBackend, blifuncname)

    return quote

        """
            scal!(n, α, x, incx)
        BLIS-based SCAL. This is the same as BLAS' `?scal`.
        """
        $(esc(targetfunc))(n::Integer,
                           α::$T1,
                           x::Union{Ptr{$T1}, Vector{$T1}},
                           incx::Integer) = begin

            $blifunc(BLIS_NO_CONJUGATE, n, α, x, incx)
            x

        end
    end
end

@blis_interface_linalg_lv1_scal Float32    scal! scalv!
@blis_interface_linalg_lv1_scal Float64    scal! scalv!
@blis_interface_linalg_lv1_scal ComplexF32 scal! scalv!
@blis_interface_linalg_lv1_scal ComplexF64 scal! scalv!

macro blis_interface_linalg_lv1_copy(T1, targetfunc, bliname)

    # Get method for the typed API backend.
    blifuncname = Symbol("bli_", cblas_typechar[ eval(T1) ], bliname)
    blifunc = getproperty(TypedBackend, blifuncname)

    return quote

        """
            copy!(n, x, incx, y, incy)
        BLIS-based COPY. This is the same as BLAS' `?copy`.
        """
        $(esc(targetfunc))(n::Integer,
                           x::Union{Ptr{$T1}, Vector{$T1}},
                           incx::Integer,
                           y::Union{Ptr{$T1}, Vector{$T1}},
                           incy::Integer) = begin

            $blifunc(BLIS_NO_CONJUGATE, n, x, incx, y, incy)
            y

        end
    end
end

@blis_interface_linalg_lv1_copy Float32    blascopy! copyv!
@blis_interface_linalg_lv1_copy Float64    blascopy! copyv!
@blis_interface_linalg_lv1_copy ComplexF32 blascopy! copyv!
@blis_interface_linalg_lv1_copy ComplexF64 blascopy! copyv!

macro blis_interface_linalg_lv1_axpy(Tc1, T1, T2, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            axpy!(α, x, y)
        BLIS-based AXPY.
        """
        $(esc(targetfunc))(α::$Tc1,
                           x::StridedVector{$T1},
                           y::StridedVector{$T2}) = begin

            (length(x) == length(y)) || throw(DimensionMismatch("Vector size mismatch."))

            oα = BliObj(α)
            ox = BliObj(x)
            oy = BliObj(y)

            $blifunc(oα, ox, oy)
            y

        end
    end
end

@blis_interface_linalg_lv1_axpy Float32    Float32    Float32    axpy! axpyv!
@blis_interface_linalg_lv1_axpy Float64    Float64    Float64    axpy! axpyv!
@blis_interface_linalg_lv1_axpy Float32    Float32    Float64    axpy! axpyv!
@blis_interface_linalg_lv1_axpy ComplexF32 ComplexF32 ComplexF32 axpy! axpyv!
@blis_interface_linalg_lv1_axpy ComplexF64 ComplexF64 ComplexF64 axpy! axpyv!
@blis_interface_linalg_lv1_axpy Float64    ComplexF64 ComplexF64 axpy! axpyv!
@blis_interface_linalg_lv1_axpy ComplexF64 Float64    ComplexF64 axpy! axpyv!

macro blis_interface_linalg_lv1_axpby(Tc1, T1, Tc2, T2, targetfunc, bliname)

    # Get method for the object API backend.
    blifuncname = Symbol("bli_", bliname)
    blifunc = getproperty(ObjectBackend, blifuncname)

    return quote

        """
            axpby!(α, x, β, y)
        BLIS-based AXPBY.
        """
        $(esc(targetfunc))(α::$Tc1,
                           x::StridedVector{$T1},
                           β::$Tc2,
                           y::StridedVector{$T2}) = begin

            (length(x) == length(y)) || throw(DimensionMismatch("Vector size mismatch."))

            oα = BliObj(α)
            ox = BliObj(x)
            oβ = BliObj(β)
            oy = BliObj(y)

            $blifunc(oα, ox, oβ, oy)
            y

        end
    end
end

@blis_interface_linalg_lv1_axpby Float32    Float32    Float32    Float32    axpby! axpbyv!
@blis_interface_linalg_lv1_axpby Float64    Float64    Float64    Float64    axpby! axpbyv!
@blis_interface_linalg_lv1_axpby Float32    Float32    Float64    Float64    axpby! axpbyv!
@blis_interface_linalg_lv1_axpby ComplexF32 ComplexF32 ComplexF32 ComplexF32 axpby! axpbyv!
@blis_interface_linalg_lv1_axpby ComplexF64 ComplexF64 ComplexF64 ComplexF64 axpby! axpbyv!
@blis_interface_linalg_lv1_axpby Float64    Float64    ComplexF64 ComplexF64 axpby! axpbyv!
@blis_interface_linalg_lv1_axpby Float64    ComplexF64 ComplexF64 ComplexF64 axpby! axpbyv!
@blis_interface_linalg_lv1_axpby ComplexF64 Float64    ComplexF64 ComplexF64 axpby! axpbyv!

