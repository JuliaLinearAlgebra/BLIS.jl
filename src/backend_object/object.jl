# Defines BLIS Objects and related methods.
#

export BliObj, BliObjBase

" BliObjBase resembles obj_t of BLIS' object-based interface. "
mutable struct BliObjBase
    # Basic fields.
    root::Ptr{BliObjBase} 
    off::NTuple{2, BliDim}
    dim::NTuple{2, BliDim}
    diag_off::BliDoff
    info1::BliObjBits
    info2::BliObjBits
    enum_size::BliSiz
    buffer::Ptr{Nothing}
    rs::BliInc
    cs::BliInc
    is::BliInc

    # Bufferless scalar.
    scalar::BliAtomic

    # Pack-related fields.
    m_padded::BliDim
    n_padded::BliDim
    ps::BliInc
    pd::BliInc
    m_panel::BliDim
    n_panel::BliDim

    # Empty inner-constructor only,
    #  resembling a plain C struct.
    BliObjBase() = begin
        obj = new(Ptr{Nothing}(), 
                  (0, 0), 
                  (0, 0), 
                  0, 0, 0, 0,
                  Ptr{Nothing}(), 0, 0, 0,
                  0,
                  0, 0, 0, 0,
                  0, 0)
    end
end

"""
Julia object BliObj packs obj_t with a reference to buffer array.
This is to preserve the base array from being recycled by GC module.
"""
struct BliObj
    obj::BliObjBase
    # Referee of obj.buffer, which itself is also a reference.
    # Purpose is to keep data alive.
    data::RefValue
end


" Object from strided matrix (view). "
BliObj(A::StridedMatrix{T}) where {T<:BliCompatibleType} = begin
    # Query data format.
    mA, nA = size(A)
    rsA, csA = strides(A)
    dtA = ctype_to_bli_num[T]

    # Initialize object with reference.
    objA = BliObj(BliObjBase(), Ref(A))
    # Pass address of A to object initializer.
    bli_obj_create_with_attached_buffer!(dtA, mA, nA,
                                         unsafe_convert(Ptr{Nothing}, A), rsA, csA,
                                         objA.obj)

    objA
end

" Object from strided vector. "
BliObj(x::StridedVector{T}) where {T<:BliCompatibleType} = begin
    mx, = size(x)
    incx, = strides(x)
    dt = ctype_to_bli_num[T]

    objx = BliObj(BliObjBase(), Ref(x))
    # Creation of vector buffer is roughly the same as matrix.
    bli_obj_create_with_attached_buffer!(dt, mx, 1,
                                         unsafe_convert(Ptr{Nothing}, x), incx, mx*incx,
                                         objx.obj)

    objx
end

" Object from scalar number. "
BliObj(γ::T) where {T<:BliCompatibleType} = begin
    dtγ = ctype_to_bli_num[T]

    objγ = BliObj(BliObjBase(), Ref(nothing))
    # bli_obj_create_1x1_with_attached_buffer!(dtγ,
    #                                          unsafe_convert(Ptr{Nothing}, objγ.data),
    #                                          objγ.obj)
    # Instead of something like above, use internal scalar storage.
    bli_obj_scalar_init_detached!(dtγ, objγ.obj)
    bli_setsc!(real(γ), imag(γ), objγ.obj)

    objγ
end


#----------------------------------------------------------
# Object-creation methods are exposed for the in-principle
#  consistency with BLIS C API and for advanced usages.
# Usually object constructors should be enough for Julia.
#
# NOTE: bli_obj_create, bli_obj_free, bli_obj_alloc_buffer,
#       bli_obj_create_1x1 and bli_obj_create_conf_to.
#       are not interfaced as Julia is a GC language.
#

bli_obj_create_without_buffer!(dt::BliNum,
                               m::BliDim,
                               n::BliDim,
                               obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_create_without_buffer),
                                                        Nothing,
                                                        (BliNum,
                                                         BliDim,
                                                         BliDim,
                                                         Ptr{BliObjBase}),
                                                        dt, m, n, Ref(obj))

bli_obj_attach_buffer!(p::Ptr{Nothing},
                       rs::BliInc,
                       cs::BliInc,
                       is::BliInc,
                       obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_attach_buffer),
                                                Nothing,
                                                (Ptr{Nothing},
                                                 BliInc,
                                                 BliInc,
                                                 BliInc,
                                                 Ptr{BliObjBase}),
                                                p, rs, cs, is, Ref(obj))

bli_obj_create_with_attached_buffer!(dt::BliNum,
                                     m::BliDim,
                                     n::BliDim,
                                     p::Ptr{Nothing},
                                     rs::BliInc,
                                     cs::BliInc,
                                     obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_create_with_attached_buffer),
                                                              Nothing,
                                                              (BliNum,
                                                               BliDim,
                                                               BliDim,
                                                               Ptr{Nothing},
                                                               BliInc,
                                                               BliInc,
                                                               Ptr{BliObjBase}),
                                                              dt, m, n, p, rs, cs, Ref(obj))

bli_obj_create_1x1_with_attached_buffer!(dt::BliNum,
                                         p::Ptr{Nothing},
                                         obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_create_1x1_with_attached_buffer),
                                                                  Nothing,
                                                                  (BliNum,
                                                                   Ptr{Nothing},
                                                                   Ptr{BliObjBase}),
                                                                  dt, p, Ref(obj))

bli_obj_scalar_init_detached!(dt::BliNum,
                              obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_scalar_init_detached),
                                                       Nothing,
                                                       (BliNum,
                                                        Ptr{BliObjBase}),
                                                       dt, Ref(obj))

bli_setsc!(zeta_r::Real,
           zeta_i::Real,
           obj::BliObjBase) = ccall(dlsym(libblis, :bli_setsc),
                                    Nothing,
                                    (Float64,
                                     Float64,
                                     Ptr{BliObjBase}),
                                    Float64(zeta_r), Float64(zeta_i), Ref(obj))

#----------------------------------------------------------------
# These methods replicates BLIS' C marcos for setting obj_t.info.
# They are exposed in a minimum fashion.

bli_obj_set_info_bits!(dest::BliObjBits,
                       src::BliObjBits,
                       obj::BliObjBase) = begin
    obj.info1 = UInt32( ( obj.info1 & ~dest ) | src )

    nothing
end

bli_obj_set_conjtrans!(t::BliTrans,
                       obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_CONJTRANS_BITS,
                                                                 BliObjBits(t.enum),
                                                                 obj)

bli_obj_set_onlytrans!(t::BliTrans,
                       obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_TRANS_BIT,
                                                                 BliObjBits(t.enum),
                                                                 obj)

bli_obj_set_conj!(c::BliConj,
                  obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_CONJ_BIT,
                                                            BliObjBits(c.enum),
                                                            obj)

bli_obj_set_uplo!(ul::BliUpLo,
                  obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_UPLO_BITS,
                                                            BliObjBits(ul.enum),
                                                            obj)

bli_obj_set_diag!(d::BliDiag,
                  obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_DIAG_BIT,
                                                            BliObjBits(d.enum),
                                                            obj)

bli_obj_set_invert_diag!(id::BliInvDiag,
                         obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_INVERT_DIAG_BIT,
                                                                   BliObjBits(id.enum),
                                                                   obj)

bli_obj_set_struc!(st::BliStruc,
                   obj::BliObjBase) = bli_obj_set_info_bits!(BLIS_STRUC_BITS,
                                                             BliObjBits(st.enum),
                                                             obj)

