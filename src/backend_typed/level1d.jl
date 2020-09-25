# Wrapper for level-1 BLIS matrix-diagonal routines.
#

@blis_ccall_group("addd",  Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType}, BliInc, BliInc)) 
@blis_ccall_group("axpyd", Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("copyd", Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("invertd",Cvoid,(BliDoff,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("scald", Cvoid, (BliConj,
                                   BliDoff,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("scal2d",Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("setd",  Cvoid, (BliConj,
                                   BliDoff,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("setid", Cvoid, (BliDoff,
                                   BliDim,
                                   BliDim,
                                   Ptr{rType},
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("shiftd",Cvoid, (BliDoff,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("subd",  Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("xpbyd", Cvoid, (BliDoff,
                                   BliDiag,
                                   BliTrans,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc, BliInc))

