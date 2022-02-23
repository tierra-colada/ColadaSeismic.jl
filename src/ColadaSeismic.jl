module ColadaSeismic

using Logging
using JUDI
using PyCall
using LinearAlgebra

include("TimeModeling/inner/util.jl")
include("TimeModeling/inner/options.jl")
include("TimeModeling/inner/prepare_JUDI_model.jl")
include("TimeModeling/inner/prepare_JUDI_geometry.jl")
include("TimeModeling/JUDI_modeling.jl")



# ONLY WHEN COLADA IS NOT AVAILABLE
h5geo = pyimport("h5geopy._h5geo")
h5geo.sr.setSpatialReferenceFromUserInput("EPSG:32056")



function h5wavemodeling2d()
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  model_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/models.h5"
  model_name = "model2d"

  geom_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/wavefields.h5"
  geom_name = "model2d wf"

  modelCnt = h5geo.openSeisContainerByName(model_filename)
  if isnothing(modelCnt)
    @error "Unable to open model container: $model_filename"
    return
  end

  model = modelCnt.openSeis(model_name)
  if isnothing(model)
    @error "Unable to open model: $model_name"
    return
  end

  geomCnt = h5geo.openSeisContainerByName(geom_filename)
  if isnothing(geomCnt)
    @error "Unable to open geom container: $geom_filename"
    return
  end

  geom = geomCnt.openSeis(geom_name)
  if isnothing(geom)
    @error "Unable to open geom: $geom_name"
    return
  end

  opt = JUDI.Options()
  h5wavemodeling2d(model, geom, opt=opt)
end


function h5wavemodeling3d()
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  # MUST BE DELETED
  h5geo.sr.setSpatialReferenceFromUserInput("EPSG:32056")

  model_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/models.h5"
  model_name = "model3d"

  geom_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/wavefields.h5"
  geom_name = "model3d wf"

  modelCnt = h5geo.openSeisContainerByName(model_filename)
  if isnothing(modelCnt)
    @error "Unable to open model container: $model_filename"
    return
  end

  model = modelCnt.openSeis(model_name)
  if isnothing(model)
    @error "Unable to open model: $model_name"
    return
  end

  geomCnt = h5geo.openSeisContainerByName(geom_filename)
  if isnothing(geomCnt)
    @error "Unable to open geom container: $geom_filename"
    return
  end

  geom = geomCnt.openSeis(geom_name)
  if isnothing(geom)
    @error "Unable to open geom: $geom_name"
    return
  end

  opt = JUDI.Options()
  # judi_model, orientation = prepare_JUDI_model_3D(model)
  # srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_3D(geom)

  # a = 0
  h5wavemodeling3d_segy(model, geom, opt=opt)
end

h5wavemodeling2d()
# h5wavemodeling3d()

end # module
