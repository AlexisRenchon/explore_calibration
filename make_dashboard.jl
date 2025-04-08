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

# example y: y_data[iteration]["swu"]
# example g: g_data[iteration][ensemble]["lhf"]

lons = map(x -> x[1], locations)
lats = map(x -> x[2], locations)

########## Web dashboard using data in memory #############

fig = Figure()
ax_y = Axis(fig[1,1])
ax_g = Axis(fig[1,2])

function update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, g_data, y_data, lons, lats)
    m_v = menu_var.value
    m_i = menu_iter.value
    m_m = menu_m.value
    m_s = menu_season.value

    g = @lift(g_data[$m_i][$m_m][$m_v][1:4:end])
    y = @lift(y_data[$m_i][$m_v][1:4:end])

    p_g = contourf!(ax_g, lons, lats, g)
    p_y = contourf!(ax_y, lons, lats, y)
    return fig
end

app = App() do
    menu_var = Dropdown(["lhf", "shf", "swu", "lwu"])
    menu_iter = Dropdown(1:n_iterations)
    menu_m = Dropdown(1:n_ensembles)
    menu_season = Dropdown(["winter", "spring", "summer", "fall"])
    maps = update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, g_data, y_data, lons, lats)
    return DOM.div(menu_var, menu_iter, menu_m, menu_season, maps)
end

# To do
# y - g (bias) in 3rd axis
# gamma (noise) in 4th axis
