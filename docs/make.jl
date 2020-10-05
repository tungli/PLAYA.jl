using Documenter, PLAYA

makedocs(modules=[PLAYA],
         format=Documenter.HTML(),
         sitename="PLAYA.jl",
         pages = Any["Home" => "index.md",
                    ],
         doctest = true,
        )
