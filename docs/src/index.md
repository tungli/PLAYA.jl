# PLAYA

# Documentation

To add Package use:
```
julia> Pkg.add("https://github.com/tungli/PLAYA.jl.git")
```

A try to replicate a study of air kinetics from [^1] is in `demo/air-kinetics.jl`.


```@meta
CurrentModule = PLAYA
```

```@docs
parse_reactions
safe_replace
parse_bolos_output
keep_fixed!
```


[^1]: A. Obrusník, P. Bílek, T. Hoder, M. Šimek, Z. Bonaventura, Plasma Sources Sci. Technol., 2018, DOI:10.1088/1361-6595/aad663. 
