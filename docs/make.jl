using Documenter, Monstera

DocMeta.setdocmeta!(Monstera, :DocTestSetup, :(using Monstera); recursive=true)

makedocs(;
    sitename = "Monstera.jl",
    # modules = [Monstera],
    doctest = false
    )
