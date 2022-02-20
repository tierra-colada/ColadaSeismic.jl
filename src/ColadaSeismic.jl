module ColadaSeismic

using Logging
using JUDI
using PyCall


include("TimeModeling/JUDI_modeling_2D.jl")


function h5wavemodeling2d()
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  model_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/models.h5"
  model_name = "model2d"

  geom_filename = "/home/kerim/Documents/Colada_prj/default/DATA/seismic/wavefields.h5"
  geom_name = "model2d wf"

  modelCnt = h5geo.openSeisContainerByName(model_filename)
  if isnothing(modelCnt)
    print("Unable to open model container: $model_filename")
    return
  end

  model = modelCnt.openSeis(model_name)
  if isnothing(model)
    print("Unable to open model: $model_name")
    return
  end

  geomCnt = h5geo.openSeisContainerByName(geom_filename)
  if isnothing(geomCnt)
    print("Unable to open geom container: $geom_filename")
    return
  end

  geom = geomCnt.openSeis(geom_name)
  if isnothing(geom)
    print("Unable to open geom: $geom_name")
    return
  end

  opt = JUDI.Options()
  h5wavemodeling2d(model, geom, opt=opt)
end

h5wavemodeling2d()

end # module
