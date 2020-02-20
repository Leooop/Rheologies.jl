using Documenter, Rheologies

makedocs(;
    modules=[Rheologies],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/Leooop/Rheologies.jl/blob/{commit}{path}#L{line}",
    sitename="Rheologies.jl",
    authors="Léo Petit",
    assets=String[],
)

deploydocs(;
    repo="github.com/Leooop/Rheologies.jl",
)
