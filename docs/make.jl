using CellMetabolism
using Documenter

DocMeta.setdocmeta!(CellMetabolism, :DocTestSetup, :(using CellMetabolism); recursive=true)

makedocs(;
    modules=[CellMetabolism],
    authors="Denis Titov <titov@berkeley.edu>  and contributors",
    sitename="CellMetabolism.jl",
    format=Documenter.HTML(;
        canonical="https://Denis-Titov.github.io/CellMetabolism.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Denis-Titov/CellMetabolism.jl",
    devbranch="main",
)
