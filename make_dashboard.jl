import EnsembleKalmanProcesses as EKP
using WGLMakie
using Bonito
using JLD2

########### Load raw data ##################################

# Load eki object, locations, prior
eki = JLD2.load_object("eki_file.jld2")
prior = include("priors.jl")
@load "all_locations.jld2" locations

########### Get what we need in memory #####################

# Get n_ensembles, n_iterations, errors
n_ensembles = EKP.get_N_ens(eki)
n_iterations = EKP.get_N_iterations(eki)
errors = eki.error
normalized_errors = errors ./ errors[1] .* 100

# Get all g
g_all = [EKP.get_g(eki, i) for i in 1:n_iterations]

# Get all y
obs_series = EKP.get_observation_series(eki)
y_obs = obs_series.observations
y_all = [EKP.get_obs(y_obs[i]) for i in 1:n_iterations]

# Get all gamma
#get_obs_noise

# Get all contrained parameters
#get_ϕ

# to consider: Δ_t, rng, ... what else might be insightful?

########## Process data for plotting ######################
"""
    n: 1=lhf, 2=shf, 3=swu, 4=lwu
"""
function slicevar(v, n)
    idx = vcat([(i+(n-1)*4:i+(n-1)*4+3) for i in 1:16:length(v)]...)
    return v[idx]
end

y_data = Dict()
g_data = Dict()

for iteration_n in 1:length(y_all)
    y_data[iteration_n] = Dict(
        "lhf" => slicevar(y_all[iteration_n], 1),
        "shf" => slicevar(y_all[iteration_n], 2),
        "swu" => slicevar(y_all[iteration_n], 3),
        "lwu" => slicevar(y_all[iteration_n], 4),
    )
end

for iteration_n in 1:length(g_all)
    g = g_all[iteration_n]
    idxs = Dict(
        "lhf" => slicevar(1:size(g, 1), 1),
        "shf" => slicevar(1:size(g, 1), 2),
        "swu" => slicevar(1:size(g, 1), 3),
        "lwu" => slicevar(1:size(g, 1), 4),
    )
    g_data[iteration_n] = Dict()
    for ensemble in 1:size(g, 2)
        g_data[iteration_n][ensemble] = Dict(
            var => g[idxs[var], ensemble] for var in keys(idxs)
        )
    end
end

function add_seasonal_access(data_dict)
    result = Dict{Any, Any}()

    for (site_id, site_data) in data_dict
        result[site_id] = Dict{String, Any}()

        for (var_key, values) in site_data
            # Create a new dictionary for each variable with seasonal slices
            result[site_id][var_key] = Dict{String, Vector{Float64}}(
                "winter" => values[1:4:end],
                "spring" => values[2:4:end],
                "summer" => values[3:4:end],
                "fall" => values[4:4:end]
            )
        end
    end

    return result
end

# Create new structure with seasonal access
seasonal_y_data = add_seasonal_access(y_data)

function add_seasonal_access_g(data_dict)
    result = Dict{Any, Any}()

    for (outer_key, outer_val) in data_dict
        result[outer_key] = Dict{Any, Any}()

        for (middle_key, middle_val) in outer_val
            result[outer_key][middle_key] = Dict{String, Any}()

            for (var_key, values) in middle_val
                # Create seasonal slices for each variable
                result[outer_key][middle_key][var_key] = Dict{String, Vector{Float64}}(
                    "winter" => values[1:4:end],
                    "spring" => values[2:4:end],
                    "summer" => values[3:4:end],
                    "fall" => values[4:4:end]
                )
            end
        end
    end

    return result
end

# Create new structure with seasonal access
seasonal_g_data = add_seasonal_access_g(g_data)

# example y: seasonal_y_data[iteration]["swu"]["fall"]
# example g: seasonal_g_data[iteration][ensemble]["lhf"]["winter"]

lons = map(x -> x[1], locations)
lats = map(x -> x[2], locations)

########## Web dashboard using data in memory #############

fig = Figure(size = (1600, 900), fontsize = 14)
ax_y = Axis(fig[1,1],
            title="y (era5 target)",
    xlabel="Longitude",
    ylabel="Latitude")
ax_g = Axis(fig[1,2],
            title="g (model output)",
    xlabel="Longitude",
    ylabel="Latitude")

function update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, seasonal_g_data, seasonal_y_data, lons, lats)
    m_v = menu_var.value
    m_i = menu_iter.value
    m_m = menu_m.value
    m_s = menu_season.value

    g = @lift(seasonal_g_data[$m_i][$m_m][$m_v][$m_s])
    y = @lift(seasonal_y_data[$m_i][$m_v][$m_s])

    min_p = @lift(minimum(vcat($g, $y)))
    max_p = @lift(maximum(vcat($g, $y)))
    limits_p = @lift(($min_p, $max_p))

    p_g = heatmap!(ax_g, lons, lats, g, colorrange = limits_p)
    p_y = heatmap!(ax_y, lons, lats, y, colorrange = limits_p)

    cb_g = Colorbar(fig[1,3], colorrange = limits_p)
    return fig
end

app = App() do
    menu_var = Dropdown(["lhf", "shf", "swu", "lwu"])
    menu_iter = Dropdown(1:n_iterations)
    menu_m = Dropdown(1:n_ensembles)
    menu_season = Dropdown(["winter", "spring", "summer", "fall"])
    maps = update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, seasonal_g_data, seasonal_y_data, lons, lats)
    return DOM.div(menu_var, menu_iter, menu_m, menu_season, maps)
end
# http://localhost:9384/browser-display

# To do
# y - g (bias) in 3rd axis
# gamma (noise) in 4th axis

# NOTE:
# all info from eki object can be dynamically displayed - parameters, config, methods...
# infinity of dashboards can be hosted on server and publicly served - can explore different calibration (eki objects)
# this is scalable, expandable and generalizable, dashboard made from eki and served, could be automated in CI
# not esthetic right now, but can be as pretty as Manuscript figures (just not a prio)

