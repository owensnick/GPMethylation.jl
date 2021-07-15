using Documenter
using GPMethylation

makedocs(
    sitename = "GPMethylation",
    format = Documenter.HTML(),
    modules = [GPMethylation]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/owensnick/GPMethylation.jl.git",
)
