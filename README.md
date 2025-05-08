First, we need to load the data we need.
Currently:
1. eki_file.jld2
2. training_locations
3. priors.jl
4. variable_list

Ideally, these should all be in eki_file.jld2. 
Locations (coordinates), variables (e.g., "lhf", "swu", ...), and seasons ("DJF", "MAM", ...) should be in eki metadata (or names)

Now that we have the data, the goal is to be able to access specifically what we want to plot
For example, seasonal_g_data[iteration][member][variable][season] (e.g., [2][3]["lhf"]["DJF"])
To get there, currently, we generate a nested Dict, via functions with indexing. It is specific to ClimaLand flattening of diagnostics.

Ideally, we want to do that directly from eki object, so it can work for ClimaLand, ClimaAtmos, others...
To be able to do this, the eki object needs to contain all info (variable, coordinate, season).
Currently it doesn't, but it could via the new EKP Observation metadata (and or via EKP.Observation.names)
