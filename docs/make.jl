using Documenter
using BLIS

makedocs(modules=[BLIS],
            sitename = "BLIS.jl",
            authors = "RuQing G. Xu (This Wrapper)",
            pages = [
                "Home" => [
                    "index.md",
                    "interface_linalg.md",
                    "backend_object.md",
                    "backend_typed.md",
                    "performance.md"
                ]
            ]
        )

deploydocs(
    repo = "github.com/JuliaLinearAlgebra/BLIS.jl.git",
    devbranch = "master"
)

