BLIS.jl
=======

This repository to provides:

- Wrapper for typed interface of [BLIS](https://github.com/flame/blis).
- As BLIS itself is using actively templates and macros, 
  this package also aims to maximize usage of Julia's
  metaprogramming features.
  
**Short-term Roadmaps**

- Define interface methods like `Base.*` and `Base.mul!` for basic
  Julia types like `Array` and `StridedArray`.

**Long-term Roadmaps**

- Improve template instantiation style for defining BLIS typed interface.
- Provide support also for object-like interface.
