# Defines BLIS Objects and related methods.
#

" BliObjBase resembles obj_t of BLIS' object-based interface. "
struct BliObjBase
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
    # Referee of obj.buffer, which itself can also be a reference.
    data::Any
end


" Object from strided matrix (view). "
BliObj(A::StridedMatrix) = begin
    # Query data format.
    mA, nA = size(A)
    rsA, csA = strides(A)
    dtA = ctype_to_bli_num[eltype(A)]

    # Initialize object with reference.
    objA = new(BliObjBase(), A)
    # NOTE: though A could not be immutable, pass always objA.data for safety.
    bli_obj_create_with_attached_buffer!(dtA, mA, nA, objA.data, rsA, csA, objA, obj)

    objA
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
                                                        dt, m, n, obj)

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
                                                p, rs, cs, is, obj)

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
                                                              dt, m, n, p, rs, cs, obj)

bli_obj_create_1x1_with_attached_buffer!(dt::BliNum,
                                         p::Ptr{Nothing},
                                         obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_create_1x1_with_attached_buffer),
                                                                  Nothing,
                                                                  (BliNum,
                                                                   Ptr{Nothing},
                                                                   Ptr{BliObjBase}),
                                                                  dt, p, obj)

bli_obj_scalar_init_detached!(dt::BliNum,
                              obj::BliObjBase) = ccall(dlsym(libblis, :bli_obj_scalar_init_detached),
                                                       Nothing,
                                                       (BliNum,
                                                        Ptr{BliObjBase}),
                                                       dt, obj)

