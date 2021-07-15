
"""
    loaddata(file)

Loads methylation data.

  
`file` is a `RData` file containing a matrix named `betas.dasen` of size `(n x m)` for `n` samples and `m` probes.


Sample names and probe names taken from colnames and rownames respectively.


"""
function loaddata(file)
    
    R"""
        load($file)
        samples <- colnames(betas.dasen)
        probes  <- rownames(betas.dasen)
        beta    <- betas.dasen
    """
    
    @rget samples probes beta
    samples, probes, beta
    
end


"""
    loadsamplemeta(file)

Loads sample meta, a tab or commma separated file that must contain the following columns (in any order):

| Sample    | Sex   | Age   |
|-----------|-------|-------|
|  Sample_1 | F     | 1.0   |
|  Sample_2 | M     | 2.0   |
|  Sample_3 | M     | 10.0  |

- `Sample` unique string identifier
- `Sex`  {M, F} or missing
- `Age` floating point age



"""
loadsamplemeta(file) = loadmeta(file, [:Sample, :Sex, :Age], "sample")


"""
    loadprobemeta(file)

Loads probe meta, a tab or commma separated file that must contain the following columns (in any order):

| Probe      |
|------------| 
| cg00000029 |
| cg00000103 |
| cg00000109 |

"""
loadprobemeta(file) = loadmeta(file, [:Probe], "probe")
    


"""
    loadmeta(file, mandatorycols, label)

    Helper function used to load metafiles
"""
function loadmeta(file, mandatorycols, label)
    meta = CSV.read(file, DataFrame)
    
    if length(intersect(mandatorycols, propertynames(meta))) != length(mandatorycols)
        missingcols = setdist(mandatorycols, propertynames(meta))
        colstring = join(names(meta), ", ")
        mcolstring = join(missingcols, ", ")
        error("Loading $label meta. Found cols: $colstring\nMissing: $mcolstring")
    end
    meta
end
