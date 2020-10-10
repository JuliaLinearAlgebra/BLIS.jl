# Object-based wrapper for level-1 BLIS fused-vector routines.
#

@blis_ccall_object(axpy2v,
                   Cvoid,
                   αx, xObj,
                   αy, xObj,
                   x,  xObj,
                   y,  xObj,
                   z,  xObj)
@blis_ccall_object(axpyf,
                   Cvoid,
                   α, xObj,
                   a, xObj,
                   x, xObj,
                   y, xObj)
@blis_ccall_object(dotxf,
                   Cvoid,
                   α, xObj,
                   a, xObj,
                   x, xObj,
                   β, xObj,
                   y, xObj)
# TODO: Document of dotaxpyv and dotxaxpyf is not consistent
#       with header files. Prospone.

