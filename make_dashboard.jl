import EnsembleKalmanProcesses as EKP
import GeoMakie as GM
using WGLMakie
using Bonito
using JLD2
using Statistics
using Printf

########### Load raw data ##################################

# Load eki object, locations, prior
eki = JLD2.load_object("ekifiles/eki_file_all_asnow_zenith.jld2")
prior = include("ekifiles/priors_all_asnow_zenith.jl")
@load "all_locations.jld2" locations

########## Some useful functions ###########################

function RMSE(x, z)
    return sqrt(mean((x.-z).^2))
end

function error_abs(x, z)
    return mean(abs.(x.-z))
end

# TODO: change function below so it works with length(n) = 1, 2, 3...
# probably just replace the 4 with length(n)
"""
    n: 1=lhf, 2=shf, 3=swu, 4=lwu
"""
function slicevar(v, n)
    idx = vcat([(i+(n-1)*4:i+(n-1)*4+3) for i in 1:16:length(v)]...)
    return v[idx]
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

# Calculate mean lat weight averaged
function cosine_weighted_global_mean(values::Vector{Float64}, lats::Vector{Float64})
    @assert length(values) == length(lats) "Length mismatch between values and latitudes"

    weights = cosd.(lats)  # cosd for degrees
    numerator = sum(values .* weights)
    denominator = sum(weights)

    return numerator / denominator
end

########### Get what we need in memory #####################

# TODO: put all that below in a function and call in inside App()

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
params = EKP.get_ϕ(prior, eki)
param_dict = Dict(i => [params[i][:, j] for j in 1:size(params[i], 2)] for i in eachindex(params))
params_name = prior.name

# to consider: Δ_t, rng, ... what else might be insightful?

########## Process data for plotting ######################

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

# Create new structure with seasonal access
seasonal_y_data = add_seasonal_access(y_data)

# Create new structure with seasonal access
seasonal_g_data = add_seasonal_access_g(g_data)

# example y: seasonal_y_data[iteration]["swu"]["fall"]
# example g: seasonal_g_data[iteration][ensemble]["lhf"]["winter"]

lons = map(x -> x[1], locations)
lats = map(x -> x[2], locations)

########## Web dashboard using data in memory #############

function update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, ax_anomalies, ax_sm, seasonal_g_data, seasonal_y_data, lons, lats)
    m_v = menu_var.value
    m_i = menu_iter.value
    m_m = menu_m.value
    m_s = menu_season.value

    g = @lift(seasonal_g_data[$m_i][$m_m][$m_v][$m_s])
    y = @lift(seasonal_y_data[$m_i][$m_v][$m_s])
    anomalies = @lift($g .- $y)
    rmse_y_g = @lift(string("RMSE = ", round(RMSE($g, $y), digits=1), " W m⁻²"))

    min_p = @lift(minimum(vcat($g, $y)))
    max_p = @lift(maximum(vcat($g, $y)))
    limits_p = @lift(($min_p, $max_p))

    min_ano = @lift(minimum($anomalies))
    max_ano = @lift(maximum($anomalies))
    # limits_ano = @lift(($min_ano, $max_ano))
    limits_ano = (-30, 30)

    p_g = heatmap!(ax_g, lons, lats, g, colorrange = limits_p)
    p_y = heatmap!(ax_y, lons, lats, y, colorrange = limits_p)
    p_ano = heatmap!(ax_anomalies, lons, lats, anomalies, colorrange = limits_ano, colormap = cgrad(:Spectral, 3, categorical = true), highclip = :cyan, lowclip = :red)

    cl = @lift($m_v * " (W m⁻²)")
    cb = Colorbar(fig[1, 3], colorrange = limits_p, label = cl, height = 300, tellheight = false)

    cl_ano = @lift($m_v * " (W m⁻²)")
    cb_ano = Colorbar(fig[2, 3], colorrange = limits_ano, label = cl_ano, height = 300, tellheight = false, colormap = cgrad(:Spectral, 3, categorical = true),highclip = :cyan, lowclip = :red)

    y_seasonal_means = @lift([cosine_weighted_global_mean(seasonal_y_data[$m_i][$m_v][season], lats) for season in ["winter", "spring", "summer", "fall"]])
    y_seasonal_means_1 = @lift([cosine_weighted_global_mean(seasonal_y_data[1][$m_v][season], lats) for season in ["winter", "spring", "summer", "fall"]])

    g_seasonal_means = @lift([cosine_weighted_global_mean(seasonal_g_data[$m_i][$m_m][$m_v][season], lats) for season in ["winter", "spring", "summer", "fall"]])
     g_seasonal_means_1 = @lift([cosine_weighted_global_mean(seasonal_g_data[1][1][$m_v][season], lats) for season in ["winter", "spring", "summer", "fall"]])

    min_sm = 0 # @lift(minimum(vcat($y_seasonal_means, $g_seasonal_means)))
    max_sm = @lift(maximum(vcat($y_seasonal_means, $g_seasonal_means)) + 10)
    limits_sm = @lift(($min_sm, $max_sm))
    line_y_1 = lines!(ax_sm, 1:4, y_seasonal_means_1, color= (:green, 0.3), linestyle = :dash)
    lines_g_1 = lines!(ax_sm, 1:4, g_seasonal_means_1, color= (:black, 0.3), linestyle = :dash)
    lines_g = lines!(ax_sm, 1:4, g_seasonal_means, color= :black)
    lines_y = lines!(ax_sm, 1:4, y_seasonal_means, color= :green)
    text!(ax_sm, 0.1, 0.1, text = rmse_y_g, align = (:left, :top), space = :relative)

    seasons = ["winter", "spring", "summer", "fall"]
    current_s = @lift(findfirst(==($m_s), seasons))
    current_s = @lift([$current_s, $current_s])
    lines!(ax_sm, current_s, [0, 1000], color= :red, linewidth = 3)
    @lift(ylims!(ax_sm, $min_sm, $max_sm))
    axislegend(ax_sm, [lines_y, lines_g], ["era5", "ClimaLand"])

    return fig
end

app = App(; title="Explore Land Calibration v0.2") do
    fig = Figure(size = (1600, 900), fontsize = 22)
    menu_var = Dropdown(["lhf", "shf", "swu", "lwu"])
    menu_iter = Dropdown(1:n_iterations)
    menu_m = Dropdown(1:n_ensembles)
    menu_season = Dropdown(["winter", "spring", "summer", "fall"])


    rmse_all = @lift(RMSE(g_all[$(menu_iter.value)][:,1], y_all[$(menu_iter.value)]))
    rmse_all_var = @lift(RMSE(g_data[$(menu_iter.value)][$(menu_m.value)][$(menu_var.value)], y_data[$(menu_iter.value)][$(menu_var.value)]))
    rmse_clm = Dict("lhf" => 20, "shf" => 15, "swu" => 25, "lwu" => 10)
    year_x = @lift(2008+$(menu_iter.value))
    title_ym = @lift("$($(menu_season.value)) $($(menu_var.value)), iteration $($(menu_iter.value)), year $($(year_x)), Y (era5 obs)")
    title_gm = @lift("$($(menu_season.value)) $($(menu_var.value)), iteration $($(menu_iter.value)), ensemble $($(menu_m.value)), year $($(year_x)), G (ClimaLand)")
    ax_y = GM.GeoAxis(
                      fig[1, 1];
                      dest = "+proj=wintri",
                      title = title_ym,
                     )
    lines!(ax_y, GM.coastlines())
    ax_g = GM.GeoAxis(
                      fig[1, 2];
                      dest = "+proj=wintri",
                      title = title_gm,
                     )
    lines!(ax_g, GM.coastlines())
        ax_anomalies = GM.GeoAxis(
                      fig[2, 2];
                      dest = "+proj=wintri",
                      title = "g - y",
                     )
    lines!(ax_anomalies, GM.coastlines())
    ylabel_sm = @lift($(menu_var.value) * " (W m⁻²)")
    ax_sm = Axis(fig[2, 1],
                 title= "Global mean by season",
                 limits=(0.99, 4.01, 0, 400),
                 ylabel= ylabel_sm,
                 xticks = (1:4, ["winter", "spring", "summer", "fall"]),
    xlabel = "Season",
                 # xticklabel = ["winter", "spring", "summer", "fall"],
                )
    maps = update_fig(menu_var, menu_iter, menu_m, menu_season, fig, ax_y, ax_g, ax_anomalies, ax_sm, seasonal_g_data, seasonal_y_data, lons, lats)
    params_current = @lift(param_dict[$(menu_iter.value)][$(menu_m.value)])
    params_initial = param_dict[1][1]
    params_relative = @lift($params_current ./ params_initial)
        return DOM.div(
        style="display: flex; flex-direction: row; flex-wrap: wrap; gap: 20px; max-width: 1600px; margin: 0 auto 0 0;",
        # Left column
        DOM.div(
            style="flex: 0 0 300px; display: flex; flex-direction: column; gap: 20px;",
            Card(
                title="Parameters",
                DOM.div(
                    style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                    DOM.div("Variable", menu_var),
                    DOM.div("Iteration", menu_iter),
                    DOM.div("Ensemble", menu_m),
                    DOM.div("Season", menu_season)
                )
            ),
            Card(
                title="Error Metrics",
                DOM.div(
                    style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                    @lift(begin
                        current_iter = $(menu_iter.value)
                        DOM.div(
                            style="display: flex; flex-direction: column; gap: 5px;",
                            DOM.div(
                                style="display: grid; grid-template-columns: 0.8fr 1.1fr 1.1fr; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                DOM.span("Iteration"),
                                DOM.span("Absolute"),
                                DOM.span("Normalized (%)")
                            ),
                            [DOM.div(
                                style=string(
                                    "display: grid; grid-template-columns: 0.8fr 1.1fr 1.1fr; ",
                                    "padding: 5px; border-bottom: 1px solid #eee; ",
                                    i == current_iter ? "background-color: #f0f8ff; font-weight: bold;" : ""
                                ),
                                DOM.span("$(i)"),
                                DOM.span(string(round(errors[i]/1e6, digits=2), " × 10⁶")),
                                DOM.div(
                                    style=string(
                                        "display: flex; align-items: center; ",
                                        "color: ", i > 1 && normalized_errors[i] < normalized_errors[i-1] ? "green" : "inherit"
                                    ),
                                    DOM.span(round(normalized_errors[i], digits=1)),
                                    i > 1 ? DOM.span(
                                        style=string(
                                            "margin-left: 5px; font-size: 0.8em; ",
                                            "color: ", normalized_errors[i] < normalized_errors[i-1] ? "green" : "red"
                                        ),
                                        normalized_errors[i] < normalized_errors[i-1] ? "↓" : "↑"
                                    ) : DOM.span("")
                                )
                            ) for i in 1:length(errors)]
                        )
                    end)
                )
            )
        ),
        # Middle column
        DOM.div(
            style="flex: 0 0 380px; display: flex; flex-direction: column; gap: 20px;",
            Card(
                title="Parameter Values",
                DOM.div(
                    style="display: flex; flex-direction: column; gap: 10px; padding: 10px;",
                    @lift(begin
                        param_values = $(params_current)
                        relative_values = $(params_relative)
                        DOM.div(
                            style="display: flex; flex-direction: column; gap: 5px;",
                            DOM.div(
                                style="display: grid; grid-template-columns: 1.2fr 1fr 1fr; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                DOM.span("Parameter"),
                                DOM.span("Value"),
                                DOM.span("Relative")
                            ),
                            [DOM.div(
                                style="display: grid; grid-template-columns: 1.2fr 1fr 1fr; padding: 5px; border-bottom: 1px solid #eee;",
                                DOM.span(style="font-weight: bold;", params_name[i]),
                                DOM.span(
                                    # Format based on magnitude
                                    let val = param_values[i]
                                        if abs(val) < 0.001 && val != 0
                                            @sprintf("%.3e", val)  # Scientific for very small numbers
                                        elseif abs(val) > 10000
                                            @sprintf("%.3e", val)  # Scientific for very large numbers
                                        else
                                            @sprintf("%.4g", val)  # General format with 4 significant digits
                                        end
                                    end
                                ),
                                DOM.div(
                                    style=string(
                                        "display: flex; align-items: center; ",
                                        "color: ", relative_values[i] > 1 ? "green" : (relative_values[i] < 1 ? "red" : "black")
                                    ),
                                    DOM.span(@sprintf("%.2f", relative_values[i])),  # Always 2 decimal places
                                    DOM.span(
                                        style="margin-left: 5px; font-size: 0.8em;",
                                        relative_values[i] > 1 ? "↑" : (relative_values[i] < 1 ? "↓" : "")
                                    )
                                )
                            ) for i in 1:length(params_name)]
                        )
                    end)
                )
            ),
            Card(
                title="RMSE Metrics",
                DOM.div(
                    style="display: flex; flex-direction: column; gap: 15px; padding: 10px;",
                    @lift(begin
                        current_iter = $(menu_iter.value)
                        current_var = $(menu_var.value)
                        current_m = $(menu_m.value)
                        rmse_clm_value = rmse_clm[current_var]

                        DOM.div(
                            style="display: flex; flex-direction: column; gap: 15px;",
                            # Section headers
                            DOM.div(
                                style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px;",
                                DOM.h4(style="margin: 0; color: #333; grid-column: 1;", "Overall RMSE"),
                                DOM.h4(style="margin: 0; color: #333; grid-column: 2;", "$(current_var) RMSE (CLM: $(rmse_clm_value) W m⁻²)")
                            ),
                            # Table headers
                            DOM.div(
                                style="display: grid; grid-template-columns: 0.4fr 0.6fr 0.4fr 0.6fr; gap: 15px; border-bottom: 1px solid #999; padding-bottom: 5px; font-weight: bold;",
                                DOM.span("Iter."),
                                DOM.span("RMSE (W m⁻²)"),
                                DOM.span("Iter."),
                                DOM.span("RMSE (W m⁻²)")
                            ),
                            # Table rows
                            [DOM.div(
                                style="display: grid; grid-template-columns: 0.4fr 0.6fr 0.4fr 0.6fr; gap: 15px; padding: 5px; border-bottom: 1px solid #eee;",
                                # Overall RMSE column
                                DOM.span(
                                    style=i == current_iter ? "font-weight: bold;" : "",
                                    "$(i)"
                                ),
                                DOM.div(
                                    style=string(
                                        "display: flex; align-items: center; ",
                                        i == current_iter ? "font-weight: bold;" : "",
                                        "color: ", i > 1 && RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "green" : "inherit"
                                    ),
                                    DOM.span(@sprintf("%.3f", RMSE(g_all[i][:,1], y_all[i]))),
                                    i > 1 ? DOM.span(
                                        style=string(
                                            "margin-left: 5px; font-size: 0.8em; ",
                                            "color: ", RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "green" : "red"
                                        ),
                                        RMSE(g_all[i][:,1], y_all[i]) < RMSE(g_all[i-1][:,1], y_all[i-1]) ? "↓" : "↑"
                                    ) : DOM.span("")
                                ),
                                # Variable RMSE column
                                DOM.span(
                                    style=i == current_iter ? "font-weight: bold;" : "",
                                    "$(i)"
                                ),
                                DOM.div(
                                    style=string(
                                        "display: flex; align-items: center; ",
                                        i == current_iter ? "font-weight: bold;" : "",
                                        "color: ", RMSE(g_data[i][current_m][current_var], y_data[i][current_var]) < rmse_clm_value ? "green" : "red"
                                    ),
                                    DOM.span(@sprintf("%.3f", RMSE(g_data[i][current_m][current_var], y_data[i][current_var]))),
                                    DOM.span(
                                        style="margin-left: 5px; font-size: 0.8em;",
                                        i > 1 && RMSE(g_data[i][current_m][current_var], y_data[i][current_var]) <
                                        RMSE(g_data[i-1][current_m][current_var], y_data[i-1][current_var]) ? "↓" :
                                        (i > 1 ? "↑" : "")
                                    )
                                )
                            ) for i in 1:n_iterations]
                        )
                    end)
                )
            )
        ),
        # Right column (maps)
        DOM.div(
            style="flex: 1; min-width: 800px;",
            Card(
                title="Map Visualization",
                maps
            )
        )
    )
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
# the dashboard allows access to 4 (variables) * 4 (plots) * 8 (iterations) * 19 (ensemble) *  4 (seasons) = 10,000 plots...
