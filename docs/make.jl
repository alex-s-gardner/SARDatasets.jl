using SLCDatasets
using Documenter

DocMeta.setdocmeta!(SLCDatasets, :DocTestSetup, :(using SLCDatasets); recursive=true)

makedocs(;
    modules=[SLCDatasets],
    authors="Alex S. Gardner, JPL/NASA",
    sitename="SLCDatasets.jl",
    format=Documenter.HTML(;
        canonical="https://alex-s-gardner.github.io/SLCDatasets.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alex-s-gardner/SLCDatasets.jl",
    devbranch="main",
)
