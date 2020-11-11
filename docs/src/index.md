# BLIS.jl

This package provides:

- Wrapper for typed and object-based interface of [BLIS](https://github.com/flame/blis).
- Overwrite of `LinearAlgebra.BLAS` functions so that matrix
  multiplications can also be redirected to the BLIS backend.

Julia code in this repository aims to maximize usage of Julia's metaprogramming features to minimize hard coding. Some of the features might be at the cutting edge.

## Table of contents

```@contents
Pages = [
  "index.md",
  "interface_linalg.md",
  "backend_object.md",
  "backend_typed.md",
  "performance.md"
]
Depth = 4
```

## Installation
```julia
using Pkg
Pkg.add("BLIS")
# or simply `]add BLIS`.
```
