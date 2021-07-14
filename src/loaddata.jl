
"""
    loaddata(file)

    Loads methylation data.
    `file` is a `RData` file containing a matrix of size `(n x m)` 
    for `n` samples and `m` probes.
"""
function loaddata(file)
    
    
    
end

"""
    loadsamplemeta(file)

    Loads sample meta - that must contain the following columns (in any order):

    | Sample    | Sex   | Age   |
    |-----------|-------|-------|
    |  Sample_1 | F     | 1.0   |
    |  Sample_2 | M     | 2.0   |
    |  Sample_3 | M     | 10.0  |
    

"""
function loadsamplemeta(file, mandatorycols = [:Sample, :Sex, :Age])
    meta = CSV.read(file, DataFrame)
    
    if length(intersect(mandatorycols, propertynames(meta))) != length(mandatorycols)
        missingcols = setdist(mandatorycols, propertynames(meta))
        colstring = join(names(meta), ", ")
        mcolstring = join(missingcols, ", ")
        error("Loading sample meta. Found cols: $colstring\nMissing: $mcolstring")
    end
    meta
end


"""
    loadprobemeta(file)

    Loads probe meta - a tab or comma separated file with columns:
"""
function loadprobemeta(file)

end

