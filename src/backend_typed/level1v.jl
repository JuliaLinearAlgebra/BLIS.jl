# Wrapper for level-1 BLIS vector routines.
#

@blis_ccall_group("addv",  Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc)) 
@blis_ccall_group("amaxv", Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{Clong}))
@blis_ccall_group("axpyv", Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("axpbyv",Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc))
@blis_ccall_group("copyv", Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("dotv",  Cvoid, (BliConj,
                                   BliConj,
                                   BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}))
@blis_ccall_group("dotxv", Cvoid, (BliConj,
                                   BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType},
                                   Ptr{xType}))
@blis_ccall_group("invertv",Cvoid,(BliDim,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("scalv", Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc))
@blis_ccall_group("scal2v",Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("setv",  Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc))
@blis_ccall_group("subv",  Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("swapv", Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("xpbyv", Cvoid, (BliConj,
                                   BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{xType},
                                   Ptr{xType}, BliInc))
