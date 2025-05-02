using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaComms
using ClimaCore
ClimaComms.@import_required_backends

"""
    diagnostics_lat_lon(nelements)

Return latitudes and longitude for the diagnostics of the land
model run at a given resolution (`nelements`).
"""
function diagnostics_lat_lon(nelements)
    radius = 0.1 # These don't matter
    depth = 0.1 # These don't matter

    domain = ClimaLand.Domains.SphericalShell(; radius, depth, nelements)

    num_long, num_lat, _ =
        ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
    longs = collect(range(-180.0, 180.0, length = num_long))
    lats = collect(range(-90.0, 90.0, length = num_lat))

    mask = ClimaLand.landsea_mask(domain)
    return lats, longs, mask
end

"""
    make_training_locations(nelements)

Create a list of geographic training locations (longitude, latitude) for model calibration.

It assumes that the output is on the default grid, as determined by
`ClimaLand.default_num_points`. It also assumes that the domain is the full
globe.

# Notes
- The function applies a land mask to identify suitable training locations
"""
function make_training_locations(nelements)
    lats, longs, mask = diagnostics_lat_lon(nelements)

    # Use the mask ClimaLand.landsea_mask
    target_hcoords = [ClimaCore.Geometry.LatLongPoint(lat, lon) for lat in lats,  lon in longs]
    interpolated_mask = Array(ClimaCore.Remapping.interpolate(mask; target_hcoords))

    training_locations = [
        (lon, lat) for
        (j, lon) in enumerate(longs) for
        (i, lat) in enumerate(lats) if
        !iszero(interpolated_mask[i, j])
    ]

    return training_locations
end
