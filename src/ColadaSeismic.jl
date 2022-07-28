module ColadaSeismic

using Logging
using JUDI
using SegyIO
using PyCall
using LinearAlgebra
using SlimOptim
using Random
using Statistics
using NLopt

# h5geo includes
include("h5geo/h5geo.jl")

# Colada includes
include("TimeModeling/util.jl")
include("TimeModeling/h5judi_physical_parameter.jl")
include("TimeModeling/h5judi_forward_modeling.jl")
include("TimeModeling/h5judi_rtm.jl")
include("TimeModeling/h5judi_fwi.jl")
include("TimeModeling/h5judi_twri.jl")
include("TimeModeling/h5judi_modeling.jl")

# global variables
# Julia 1.6.7: type declarations on global variables are not yet supported
model_origin_x = 0.0
model_origin_y = 0.0
model_orientation = 0.0
save_as = SEGY
spatial_reference = ""
h5geo = nothing
survey_type = nothing
seis_out_cnt = nothing
seis_out = nothing

export setGlobals

# Types are hidden so that it was possible to pass 'nothing' as argument
function setGlobals(;
  model_origin_x=nothing, 
  model_origin_y=nothing, 
  model_orientation=nothing,
  save_as=nothing,
  spatial_reference=nothing,
  h5geo=nothing,
  survey_type=nothing,
  seis_out_cnt=nothing,
  seis_out=nothing)
  if !isnothing(model_origin_x)
    global model_origin_x = model_origin_x
  end
  if !isnothing(model_origin_y)
    global model_origin_y = model_origin_y
  end
  if !isnothing(model_orientation)
    global model_orientation = model_orientation
  end
  if !isnothing(save_as)
    global save_as = save_as
  end
  if !isnothing(spatial_reference)
    global spatial_reference = spatial_reference
  end
  if !isnothing(h5geo)
    global h5geo = h5geo
  end
  if !isnothing(survey_type)
    global survey_type = survey_type
  end
  if !isnothing(seis_out_cnt)
    global seis_out_cnt = seis_out_cnt
  end
  if !isnothing(seis_out)
    global seis_out = seis_out
  end
end


end # module
