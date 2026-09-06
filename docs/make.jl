using SARDatasets
using Documenter

DocMeta.setdocmeta!(SARDatasets, :DocTestSetup, :(using SARDatasets); recursive=true)

makedocs(;
    modules=[SARDatasets],
    authors="Alex S. Gardner, JPL/NASA",
    sitename="SARDatasets.jl",
    format=Documenter.HTML(;
        canonical="https://alex-s-gardner.github.io/SARDatasets.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alex-s-gardner/SARDatasets.jl",
    devbranch="main",
)
