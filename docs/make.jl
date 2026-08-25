using SAR
using Documenter

DocMeta.setdocmeta!(SAR, :DocTestSetup, :(using SAR); recursive=true)

makedocs(;
    modules=[SAR],
    authors="Alex S. Gardner, JPL/NASA",
    sitename="SAR.jl",
    format=Documenter.HTML(;
        canonical="https://alex-s-gardner.github.io/SAR.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/alex-s-gardner/SAR.jl",
    devbranch="main",
)
