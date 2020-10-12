# Override mul! method.
#

"""
    mul!(C, A, B, α, β)
BLIS-based matrix multiply.
Generic strides are directly supported.
"""
mul!(C::StridedMatrix{TC},
     A::StridedMatrix{TA},
     B::StridedMatrix{TB},
     α::Tα,
     β::Tβ) where {TA<:BliCompatibleType,
                   TB<:BliCompatibleType,
                   TC<:BliCompatibleType,
                   Tα<:BliCompatibleType,
                   Tβ<:BliCompatibleType} = begin

    oA = BliObj(A)
    oB = BliObj(B)
    oC = BliObj(C)
    oα = BliObj(α)
    oβ = BliObj(β)

    ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
    C
end

"""
    mul!(C, A, B)
BLIS-based matrix multiply.
Generic strides are directly supported.
"""
mul!(C::StridedMatrix{TC},
     A::StridedMatrix{TA},
     B::StridedMatrix{TB}
    ) where {TA<:BliCompatibleType,
             TB<:BliCompatibleType,
             TC<:BliCompatibleType} = mul!(C, A, B, TA(1.0), TC(0.0))

#--------------------------------------------------
# Instantiate mul! separately to avoid conflicting
#  with Julia's intrinsic.

mul!(C::StridedMatrix{TC},
     A::StridedMatrix{TA},
     B::Adjoint{TB, <:StridedMatrix{TB}},
     α::Tα,
     β::Tβ) where {TA<:BliCompatibleType,
                   TB<:BliCompatibleType,
                   TC<:BliCompatibleType,
                   Tα<:BliCompatibleType,
                   Tβ<:BliCompatibleType} = begin

    oA = BliObj(A)
    oB = BliObj(B')
    oC = BliObj(C)
    oα = BliObj(α)
    oβ = BliObj(β)

    ObjectBackend.bli_obj_set_conjtrans!(BLIS_CONJ_TRANSPOSE, oB.obj)
    ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
    C
end

# mul!(C::StridedMatrix{TC},
#      A::Adjoint{TA, <:StridedMatrix{TA}},
#      B::StridedMatrix{TB},
#      α::Tα,
#      β::Tβ) where {TA<:BliCompatibleType,
#                    TB<:BliCompatibleType,
#                    TC<:BliCompatibleType,
#                    Tα<:BliCompatibleType,
#                    Tβ<:BliCompatibleType} = begin
#
#     oA = BliObj(A')
#     oB = BliObj(B)
#     oC = BliObj(C)
#     oα = BliObj(α)
#     oβ = BliObj(β)
#     @show oα
#     
#     ObjectBackend.bli_obj_set_conjtrans!(BLIS_CONJ_TRANSPOSE, oA.obj)
#     ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
#     C
# end

mul!(C::StridedMatrix{TC},
     A::Adjoint{TA, <:StridedMatrix{TA}},
     B::Adjoint{TB, <:StridedMatrix{TB}},
     α::Tα,
     β::Tβ) where {TA<:BliCompatibleType,
                   TB<:BliCompatibleType,
                   TC<:BliCompatibleType,
                   Tα<:BliCompatibleType,
                   Tβ<:BliCompatibleType} = begin

    oA = BliObj(A')
    oB = BliObj(B')
    oC = BliObj(C)
    oα = BliObj(α)
    oβ = BliObj(β)

    ObjectBackend.bli_obj_set_conjtrans!(BLIS_CONJ_TRANSPOSE, oA.obj)
    ObjectBackend.bli_obj_set_conjtrans!(BLIS_CONJ_TRANSPOSE, oB.obj)
    ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
    C
end

mul!(C::StridedMatrix{TC},
     A::Adjoint{TA, <:StridedMatrix{TA}},
     B::Adjoint{TB, <:StridedMatrix{TB}}
    ) where {TA<:BliCompatibleType,
             TB<:BliCompatibleType,
             TC<:BliCompatibleType} = mul!(C, A, B, TA(1.0), TC(0.0))

mul!(C::StridedMatrix{TC},
     A::StridedMatrix{TA},
     B::Adjoint{TB, <:StridedMatrix{TB}}
    ) where {TA<:BliCompatibleType,
             TB<:BliCompatibleType,
             TC<:BliCompatibleType} = mul!(C, A, B, TA(1.0), TC(0.0))

# A'B instantiated here as mul!(C, A', B, α, β)
#  has additional ambiguity.
mul!(C::StridedMatrix{TC},
     A::Adjoint{TA, <:StridedMatrix{TA}},
     B::StridedMatrix{TB}
    ) where {TA<:BliCompatibleType,
             TB<:BliCompatibleType,
             TC<:BliCompatibleType} = begin
    oA = BliObj(A')
    oB = BliObj(B)
    oC = BliObj(C)
    oα = BliObj(1.0)
    oβ = BliObj(0.0)

    ObjectBackend.bli_obj_set_conjtrans!(BLIS_CONJ_TRANSPOSE, oA.obj)
    ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
    C
end

