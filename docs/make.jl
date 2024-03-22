using Documenter, Monstera

DocMeta.setdocmeta!(Monstera, :DocTestSetup, :(using Monstera); recursive=true)

makedocs(;
    sitename = "Monstera",
    modules = [Monstera],
    )
