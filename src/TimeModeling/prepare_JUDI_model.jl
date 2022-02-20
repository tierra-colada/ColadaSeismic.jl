function prepare_JUDI_model_2D(
  h5model::PyCall.PyObject;
  model_xkey::String="CDP_X")

  h5geo = pyimport("h5geopy._h5geo")

  v = transpose(h5model.getTrace(0, typemax(Int), 0, typemax(Int), "km/s"))
  if isempty(v)
    @error "Unable to read model. Probably data units incorrect/missing"
    return
  end

  if ndims(v) != 2
    @error "Model must be two dimensional array"
    return
  end
  
  model_x = h5model.getTraceHeader(model_xkey, 0, typemax(Int), h5model.getLengthUnits(), "m")
  if length(model_x) < 2
    @error "$model_xkey can't have less than 2 points"
    return
  end
  
  # in case model is not sorted
  model_x_ind, model_x = h5geo.sortv(model_x)
  v = v[model_x_ind .+ 1, :]   # C++ returned indexes starts from 0
  
  # Set up model structure
  n = size(v)   # (x,y,z) or (x,z)
  d = (h5model.getSampRate("m") * (-1), model_x[2]-model_x[1])
  o = (model_x[1], h5model.getSRD("m") * (-1))
  
  # Slowness squared [s^2/km^2]
  m = (1f0 ./ v).^2
  
  return JUDI.Model(n, d, o, m)
end