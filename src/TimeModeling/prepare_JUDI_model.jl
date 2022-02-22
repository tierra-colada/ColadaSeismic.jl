

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
  d = (abs(h5model.getSampRate("m")), model_x[2]-model_x[1])
  o = (model_x[1], h5model.getSRD("m") * (-1))
  
  # Slowness squared [s^2/km^2]
  m = (1f0 ./ v).^2
  
  return JUDI.Model(n, d, o, m)
end



function prepare_JUDI_model_3D(
  h5model::PyCall.PyObject;
  model_xkey::String="CDP_X",
  model_ykey::String="CDP_Y")

  h5geo = pyimport("h5geopy._h5geo")

  keylist = ["INLINE", "XLINE"]
  minlist = [-Inf, -Inf]
  maxlist = [Inf, Inf]
  v, il_xl, ind = h5model.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "km/s")
  if isempty(v)
    @error "Unable to read model. Probably data units incorrect/missing"
    return
  end

  keylist = [model_xkey, model_ykey]
  xy = h5model.getXYTraceHeaders(keylist, ind, "m", true)
  if isempty(xy)
    @error "Unable to read $model_xkey and $model_ykey trace headers. Probably length units incorrect/missing"
    return
  end

  uil_ind, uil, uil_from_size = h5geo.sort_unique(il_xl[:,1])
  uxl_ind, uxl, uxl_from_size = h5geo.sort_unique(il_xl[:,2])
  nsamp = h5model.getNSamp()

  nil = length(uil)
  nxl = length(uxl)

  # if the number of XL of each IL is not constant return false
  if !all(uil_from_size[:,2] .== uil_from_size[1,2])
    @error "Selected traces can't be represented as cube"
    return
  end

  # p1/p2 - second point on the first IL/XL
  origin = Vector{Float32}(undef, 2)
  origin[1] = xy[1,1]
  origin[2] = xy[1,2]

  p1 = Vector{Float32}(undef, 2)
  p2 = Vector{Float32}(undef, 2)
  if nxl > 1
    p1[1] = xy[2,1]
    p1[2] = xy[2,2]
  else
    p1 = origin
  end

  if nil > 1
    p2[1] = xy[nxl+1,1]
    p2[2] = xy[nxl+1,2]
  else
    p2 = origin
  end

  orientation = 0
  if nxl > 1
    orientation = atan((p1[2]-origin[2])/(p1[1]-origin[1]))
  elseif nil > 1
    orientation = atan((p2[2]-origin[2])/(p2[1]-origin[1]))
  end

  origin_rot = rotate_2xN_array(origin, -orientation)
  p1_rot = rotate_2xN_array(p1, -orientation)
  p2_rot = rotate_2xN_array(p2, -orientation)
  dx = 0
  dy = 0
  if nxl > 1
    dx = p1_rot[1]-origin_rot[1]
    if dx < 0
      dx = abs(dx)
      origin_rot[1] = origin_rot[1]-dx
      p1_rot[1] = p1_rot[1]+dx
    end
    dy = p2_rot[2]-origin_rot[2]
    if dy < 0
      dy = abs(dy)
      origin_rot[2] = origin_rot[2]-dy
      p2_rot[2] = p2_rot[2]+dy
    end
  elseif nil > 1
    dx = p2_rot[1]-origin_rot[1]
    if dx < 0
      dx = abs(dx)
      origin_rot[1] = origin_rot[1]-dx
      p2_rot[1] = p2_rot[1]+dx
    end
    dy = p1_rot[2]-origin_rot[2]
    if dy < 0
      dy = abs(dy)
      origin_rot[2] = origin_rot[2]-dy
      p1_rot[2] = p1_rot[2]+dy
    end
  end

  v = reshape(v, nsamp, nxl, nil)
  v = permutedims(v, [2,3,1])
  
  # Set up model structure
  n = size(v)   # (z,x,y)
  d = (dx, dy, abs(h5model.getSampRate("m")))
  o = (origin_rot[1], origin_rot[2], h5model.getSRD("m") * (-1))
  
  # Slowness squared [s^2/km^2]
  m = (1f0 ./ v).^2
  
  return JUDI.Model(n, d, o, m), orientation
end