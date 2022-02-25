function H5ReadPhysicalParameter2D(
  h5obj::PyCall.PyObject;
  h5opt::H5PhPOptions)

  h5geo = pyimport("h5geopy._h5geo")

  if h5opt.phptype == VELOCITY 
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int), "km/s"))
    # Slowness squared [s^2/km^2]
    php = (1f0 ./ php).^2
  elseif h5opt.phptype == DENSITY
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int), "g/cm^3"))
  elseif h5opt.phptype == THETA || h5opt.phptype == PHI
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int), "rad"))
  else
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int)))
  end

  if isempty(php)
    @error "Unable to read PhysicalParameter. Probably data units incorrect/missing"
    return
  end

  if ndims(php) != 2
    @error "PhysicalParameter must be two dimensional array"
    return
  end
  
  x = h5obj.getTraceHeader(h5opt.xkey, 0, typemax(Int), h5obj.getLengthUnits(), "m")
  if length(x) < 2
    @error "$(h5opt.xkey) can't have less than 2 points"
    return
  end
  
  # in case PhysicalParameter is not sorted
  x_ind, x = h5geo.sortv(x)
  php = php[x_ind .+ 1, :]   # C++ returned indexes starts from 0
  
  # Set up PhysicalParameter structure
  n = size(php)   # (x,y,z) or (x,z)
  d = (abs(h5obj.getSampRate("m")), x[2]-x[1])
  o = (x[1], h5obj.getSRD("m") * (-1))
  
  return JUDI.PhysicalParameter(php, n, d, o;)
end



function H5ReadPhysicalParameter3D(
  h5obj::PyCall.PyObject;
  h5opt::H5PhPOptions)

  h5geo = pyimport("h5geopy._h5geo")

  keylist = ["INLINE", "XLINE"]
  minlist = [-Inf, -Inf]
  maxlist = [Inf, Inf]
  if h5opt.phptype == VELOCITY 
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "km/s")
    php = (1f0 ./ php).^2
  elseif h5opt.phptype == DENSITY
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "g/cm^3")
  elseif h5opt.phptype == THETA || h5opt.phptype == PHI
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "rad")
  else
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true)
  end

  if isempty(php)
    @error "Unable to read PhysicalParameter. Probably data units incorrect/missing"
    return
  end

  keylist = [h5opt.xkey, h5opt.ykey]
  xy = h5obj.getXYTraceHeaders(keylist, ind, "m", true)
  if isempty(xy)
    @error "Unable to read $(h5opt.xkey) and $(h5opt.ykey) trace headers. Probably length units incorrect/missing"
    return
  end

  _, uil, uil_from_size = h5geo.sort_unique(il_xl[:,1])
  _, uxl, _ = h5geo.sort_unique(il_xl[:,2])
  nsamp = h5obj.getNSamp()

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

  php = reshape(php, nsamp, nxl, nil)
  php = permutedims(php, [2,3,1])
  
  # Set up PhysicalParameter structure
  n = size(php)   # (x,y,z)
  d = (dx, dy, abs(h5obj.getSampRate("m")))
  o = (origin_rot[1], origin_rot[2], h5obj.getSRD("m") * (-1))
  
  return JUDI.PhysicalParameter(php, n, d, o)
end