# Wrapper for some BLIS utility routines.
#

@blis_ccall_group("asumv", Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("norm1m",Cvoid, (BliDoff,
                                   BliDoff,
                                   BliUpLo,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("normfm",Cvoid, (BliDoff,
                                   BliDoff,
                                   BliUpLo,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("normim",Cvoid, (BliDoff,
                                   BliDoff,
                                   BliUpLo,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("norm1v",Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("normfv",Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("normiv",Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{rType}))
@blis_ccall_group("mkherm",Cvoid, (BliUpLo,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("mksymm",Cvoid, (BliUpLo,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("mktrim",Cvoid, (BliUpLo,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc))
# TODO: printv functions.
@blis_ccall_group("randv", Cvoid, (BliDim,
                                   Ptr{xType}, BliInc))
@blis_ccall_group("randm", Cvoid, (BliDoff,
                                   BliUpLo,
                                   BliDim,
                                   BliDim,
                                   Ptr{xType}, BliInc, BliInc))
@blis_ccall_group("sumsqv",Cvoid, (BliDim,
                                   Ptr{xType}, BliInc,
                                   Ptr{rType},
                                   Ptr{rType}))

