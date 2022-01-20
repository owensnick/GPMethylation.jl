

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

Loads sample meta, a tab or commma separated file that must contain the following first threecolumns (in any order):

| Sample    | Sex   | Age   |  Filter  |
|-----------|-------|-------|----------|
|  Sample_1 | F     | 1.0   |   1      |
|  Sample_2 | M     | 2.0   |   0      |
|  Sample_3 | M     | 10.0  |   1      |

- `Sample` unique string identifier
- `Sex`  {M, F} or missing
- `Age` floating point age
- `Filter` optional column of `{1, 0}`s to exclude samples, when `Filter = 0` sample is excluded.

`Sample` should match `colnames` of RData file. 


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


Can be used to override the `rownames` in `RData` file.

"""
loadprobemeta(file) = loadmeta(file, [:Probe], "probe")
    


"""
    loadmeta(file, mandatorycols, label)

    Helper function used to load metafiles
"""
function loadmeta(file, mandatorycols, label)
    meta = CSV.read(file, DataFrame)
    
    if length(intersect(mandatorycols, propertynames(meta))) != length(mandatorycols)
        missingcols = setdiff(mandatorycols, propertynames(meta))
        colstring = join(names(meta), ", ")
        mcolstring = join(missingcols, ", ")
        error("Loading $label meta. Found cols: $colstring\nMissing: $mcolstring")
    end
    meta
end

"""
    loadall(samplemetafile, betafile)

Loads samplefmetafile and `RData` file `betafile` ensures samples details correspond and filters samples
"""
function loadall(samplemetafile, betafile, verbose=true)

    verbose && println("[GPM]\tLoaded sample meta file ", samplemetafile)
    verbose && println("[GPM]\tLoaded probe data file ", betafile)
    meta = loadsamplemeta(samplemetafile)
    samples, probes, beta = loaddata(betafile)



    if meta.Sample != samples
        error("Sample meta does not correspond to beta colnames")
    end

    verbose && println("[GPM]\tLoaded meta information for ", size(meta, 1), " samples.")
    verbose && println("[GPM]\tLoaded probe data on  ", size(beta, 1), " probes and ", size(beta, 2), " samples.")

    ### Filter samples
    if :Filter ∈ propertynames(meta)
        ind = meta.Filter .== 1

        if verbose
            println("[GPM]\tExcluding ", sum(.!ind), " samples:")
            display(meta[.!ind, :])
        end
        meta = meta[ind, :]
        beta = beta[:, ind]


    end

    ### Converts row major R to column major Julia
    meta, probes, Matrix(beta')
end