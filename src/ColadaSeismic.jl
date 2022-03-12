module ColadaSeismic

using Logging
using JUDI
using PyCall
using LinearAlgebra

# h5geo includes
include("h5geo/h5geo.jl")

# Colada includes
include("TimeModeling/inner/util.jl")
include("TimeModeling/inner/options.jl")
include("TimeModeling/inner/h5judi_geometry.jl")
include("TimeModeling/inner/h5judi_physical_parameter.jl")
include("TimeModeling/h5judi_forward_modeling.jl")
include("TimeModeling/h5judi_modeling.jl")



# ONLY WHEN COLADA IS NOT AVAILABLE
# h5geo = pyimport("h5geopy._h5geo")
# h5geo.sr.setSpatialReferenceFromUserInput("EPSG:32056")


end # module
