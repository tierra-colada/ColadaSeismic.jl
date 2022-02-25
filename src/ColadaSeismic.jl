module ColadaSeismic

using Logging
using JUDI
using PyCall
using LinearAlgebra

include("TimeModeling/inner/utils.jl")
include("TimeModeling/inner/options.jl")
include("TimeModeling/inner/h5judi_geometry.jl")
# include("TimeModeling/inner/h5judi_model.jl")
include("TimeModeling/inner/h5judi_physical_parameter.jl")
include("TimeModeling/h5judi_modeling.jl")



# ONLY WHEN COLADA IS NOT AVAILABLE
h5geo = pyimport("h5geopy._h5geo")
h5geo.sr.setSpatialReferenceFromUserInput("EPSG:32056")



function H5Modeling2D()
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
  # H5Modeling2D(model, geom, opt=opt)

  # judi_model, orientation = prepare_JUDI_model_3D(model)
  # srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_3D(geom)

  # a = 0

  opt.save_data_to_disk = true
  opt.file_path = "/home/kerim/Documents/Colada_prj/my JUDI test 2D"
  opt.file_name = "test"
  H5Modeling2D_segy(model, geom, opt=opt)
end


function H5Modeling3D()
  h5geo = pyimport("h5geopy._h5geo")

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

  opt.save_data_to_disk = true
  opt.file_path = "/home/kerim/Documents/Colada_prj/my JUDI test 3D"
  opt.file_name = "test"
  H5Modeling3D_segy(model, geom, opt=opt)
end

H5Modeling2D()
# H5Modeling3D()

end # module
