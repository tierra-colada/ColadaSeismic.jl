function H5ReadPhysicalParameter2D(
  h5obj::PyCall.PyObject;
  phptype::H5PhPType,
  xkey::String)

  h5geo = pyimport("h5geopy._h5geo")

  if phptype == VELOCITY 
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int), "km/s"))
    # Slowness squared [s^2/km^2]
    php = (1f0 ./ php).^2
  elseif phptype == DENSITY
    php = transpose(h5obj.getTrace(0, typemax(Int), 0, typemax(Int), "g/cm^3"))
  elseif phptype == THETA || phptype == PHI
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
  
  x = h5obj.getTraceHeader(xkey, 0, typemax(Int), h5obj.getLengthUnits(), "m")
  if length(x) < 2
    @error "$(xkey) can't have less than 2 points"
    return
  end
  
  # in case PhysicalParameter is not sorted
  x_ind, x = h5geo.sortv(x)
  php = php[x_ind .+ 1, :]   # C++ returned indexes starts from 0
  
  # Set up PhysicalParameter structure
  n = size(php)   # (x,y,z) or (x,z)
  d = (abs(h5obj.getSampRate("m")), x[2]-x[1])
  o = (x[1], h5obj.getSRD("m") * (-1))
  
  return JUDI.PhysicalParameter(php, n, d, o;), 0   # the second parameter is the orientation
end


function H5ReadPhysicalParameter3D(
  h5obj::PyCall.PyObject;
  phptype::H5PhPType,
  xkey::String,
  ykey::String)

  h5geo = pyimport("h5geopy._h5geo")

  keylist = ["INLINE", "XLINE"]
  minlist = [-Inf, -Inf]
  maxlist = [Inf, Inf]
  if !h5obj.hasPKeySort(keylist[1])
    h5obj.addPKeySort(keylist[1])
  end
  
  if phptype == VELOCITY 
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "km/s")
    php = (1f0 ./ php).^2
  elseif phptype == DENSITY
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "g/cm^3")
  elseif phptype == THETA || phptype == PHI
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true, "rad")
  else
    php, il_xl, ind = h5obj.getSortedData(keylist, minlist, maxlist, 0, typemax(Int), true)
  end

  if isempty(php)
    @error "Unable to read PhysicalParameter. Probably data units incorrect/missing"
    return
  end

  keylist = [xkey, ykey]
  xy = h5obj.getXYTraceHeaders(keylist, ind, "m", true)
  if isempty(xy)
    @error "Unable to read $(xkey) and $(ykey) trace headers. Probably length units incorrect/missing"
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

  p1_rot = rotate_2xN_array(p1-origin, -orientation) + origin
  p2_rot = rotate_2xN_array(p2-origin, -orientation) + origin
  dx = 0
  dy = 0
  if nxl > 1
    dx = p1_rot[1]-origin[1]
    if dx < 0
      dx = abs(dx)
      origin[1] = origin[1]-dx
      p1_rot[1] = p1_rot[1]+dx
    end
    dy = p2_rot[2]-origin[2]
    if dy < 0
      dy = abs(dy)
      origin[2] = origin[2]-dy
      p2_rot[2] = p2_rot[2]+dy
    end
  elseif nil > 1
    dx = p2_rot[1]-origin[1]
    if dx < 0
      dx = abs(dx)
      origin[1] = origin[1]-dx
      p2_rot[1] = p2_rot[1]+dx
    end
    dy = p1_rot[2]-origin[2]
    if dy < 0
      dy = abs(dy)
      origin[2] = origin[2]-dy
      p1_rot[2] = p1_rot[2]+dy
    end
  end

  php = reshape(php, nsamp, nxl, nil)
  php = permutedims(php, [2,3,1])
  
  # Set up PhysicalParameter structure
  n = size(php)   # (x,y,z)
  d = (dx, dy, abs(h5obj.getSampRate("m")))
  o = (origin[1], origin[2], h5obj.getSRD("m") * (-1))
  
  return JUDI.PhysicalParameter(php, n, d, o), orientation
end


function H5ReadPhysicalParameter(
  h5obj::PyCall.PyObject;
  phptype::H5PhPType,
  xkey::String,
  ykey::String)

  h5geo = pyimport("h5geopy._h5geo")

  if h5obj.getSurveyType() == h5geo.SurveyType.TWO_D
    return H5ReadPhysicalParameter2D(h5obj, phptype=phptype, xkey=xkey)
  elseif h5obj.getSurveyType() == h5geo.SurveyType.THREE_D
    return H5ReadPhysicalParameter3D(h5obj, phptype=phptype, xkey=xkey, ykey=ykey)
  else 
    @error "H5Object is neither TWO_D nor THREE_D"
    return
  end
end


function H5WritePhysicalParameter(;
  cntName::String,
  objName::String,
  cntCreationType::PyCall.PyObject,
  objCreationType::PyCall.PyObject,
  php::JUDI.PhysicalParameter)

  if isnothing(cntName) || length(cntName) < 1
    @error "Container name is empty"
    return
  end

  if isnothing(objName) || length(objName) < 1
    @error "Object name is empty"
    return
  end

  if isnothing(php) || length(php.data) < 1
    @error "PhysicalParameter is empty"
    return
  end

  colada = pyimport("colada")
  h5geo = pyimport("h5geopy._h5geo")

  cnt = h5geo.createSeisContainerByName(cntName, cntCreationType)
  if isnothing(cnt)
    @error "Unable to create SeisContainer: $cntName"
    return
  end

  authName = colada.Util().CRSAuthName()
  authCode = string(colada.Util().CRSCode())

  p = h5geo.SeisParam()
  p.spatialReference = "$authName:$authCode"
  p.lengthUnits = colada.Util().lengthUnits()
  p.temporalUnits = colada.Util().timeUnits()
  p.domain = h5geo.Domain.TVD
  p.dataType = h5geo.SeisDataType.STACK

  dims = length(size(php.data))
  TRACE = Matrix{Float32}
  x0 = 0.0; y0 = 0.0; z0 = 0.0
  dx = 0.0; dy = 0.0; dz = 0.0
  nx = 1; ny = 1; nz = 1
  if dims == 2
    p.surveyType = h5geo.SurveyType.TWO_D
    x0 = php.o[1]; z0 = -php.o[2]
    dx = php.d[1]; dz = php.d[2]
    nx = php.n[1]; nz = php.n[2]
    nTrc = nx
    nSamp = nz
    TRACE = permutedims(php.data, (2,1))
  elseif dims == 3
    p.surveyType = h5geo.SurveyType.THREE_D
    x0 = php.o[1]; y0 = php.o[2]; z0 = -php.o[3]
    dx = php.d[1]; dy = php.d[2]; dz = php.d[3]
    nx = php.n[1]; ny = php.n[2]; nz = php.n[3]
    nTrc = nx*ny
    nSamp = nz
    TRACE = reshape(permutedims(php.data, (3,2,1)), (nz, nx*ny))
  else
    @error "PhysicalParameter dimensions is neither 2D or 3D"
    return
  end

  p.nTrc = nTrc
  p.nSamp = nSamp
  
  if nTrc < 1000
    p.trcChunk = nTrc
  else 
    p.trcChunk = 1000
  end

  seis = cnt.createSeis(objName, p, objCreationType)
  if isnothing(seis)
    @error "Unable to create Seis: $objName"
    return
  end

  val = seis.generateSTKGeometry(x0, dx, nx, y0, dy, ny, z0)
  if isnothing(val)
    @error "Unable to create generate STK geometry for Seis: $objName"
    return
  end

  val = seis.writeTrace(TRACE)
  if isnothing(val)
    @error "Unable to write traces: $objName"
    return
  end

  seis.getH5File().flush()
end