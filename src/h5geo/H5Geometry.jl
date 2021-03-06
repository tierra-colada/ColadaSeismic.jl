import JUDI.Geometry, JUDI.get_nsrc, JUDI.n_samples, JUDI.subsample, JUDI.compareGeometry, JUDI.isequal

export H5GeometryOOC

mutable struct H5GeometryOOC{T} <: Geometry{T}
	container::H5SeisCon
  # JUDI GeometryOOC members start
  dt::Array{T,1}
  nt::Array{<:Integer,1}
  t::Array{T,1}
  nrec::Array{<:Integer,1}
  # JUDI GeometryOOC members end
	key::String
	xkey::String
	ykey::String
	zkey::String
  do_coord_transform::Bool
end

# mutable struct MyStruct{T}
#   a::Array{<:Integer,1}
#   b::Array{<:Integer,1}
#   c::Array{T,1}
# end

# a = ones(Int64, 2)
# b = ones(Integer, 3)
# c = ones(Float32, 3)

# MyStruct(a, b, c)

function H5GeometryOOC(;
  container::H5SeisCon,
  key::String,
  xkey::String,
  ykey::String,
  zkey::String,
  do_coord_transform::Bool)
  if isnothing(container)
    @error "H5SeisCon is Nothing\n"
    return
  end

  if isnothing(key)
    @error "key is Nothing\n"
    return
  end

  if key != "source" && key != "receiver"
    @error "key must be either 'source' or 'receiver'\n"
    return
  end

  if isnothing(xkey) || isnothing(ykey) || isnothing(zkey)
    @error "xkey/ykey/zkey is Nothing\n"
    return
  end

  nsrc = length(container)
  
  dt = zeros(Float32, nsrc) .+ Float32(abs(container.seis.getSampRate("ms")))
  nt = zeros(Integer, nsrc) .+ Integer(container.seis.getNSamp())
  t = zeros(Float32, nsrc) .+ Float32.((nt.-1).*dt)

  nrec = zeros(Integer, nsrc)
  i = 1
  for val in container.pkeyvals
    nrec[i] = Integer(container.seis.getPKeyTraceSize(container.pkey, val, val))
    i+=1
  end

  return H5GeometryOOC(container, dt, nt, t, nrec,
                      key, xkey, ykey, zkey, do_coord_transform)
end

######################## shapes easy access ################################
get_nsrc(g::H5GeometryOOC) = length(g.container.pkeyvals)

function n_samples(g::H5GeometryOOC, info::Info)
  return n_samples(g, info.nsrc)
end

n_samples(g::H5GeometryOOC, nsrc::Integer) = sum([g.nrec[j]*g.nt[j] for j=1:nsrc])

# # Subsample out-of-core geometry structure
# function subsample(geometry::H5GeometryOOC, srcnum::Number)
#   con = deepcopy(geometry.container)
#   con.pkeyvals = [con.pkeyvals[srcnum]]
#   return H5GeometryOOC(container=con, key=geometry.key,
#                       xkey=geometry.xkey, ykey=geometry.ykey, zkey=geometry.zkey, 
#                       do_coord_transform=geometry.do_coord_transform)
# end

# # srcnum maybe vector
# function subsample(geometry::H5GeometryOOC, srcnum::AbstractArray)
#   con = deepcopy(geometry.container)
#   con.pkeyvals = con.pkeyvals[srcnum]
#   return H5GeometryOOC(container=con, key=geometry.key,
#                       xkey=geometry.xkey, ykey=geometry.ykey, zkey=geometry.zkey, 
#                       do_coord_transform=geometry.do_coord_transform)
# end

# getindex out-of-core geometry structure
function getindex(geometry::H5GeometryOOC{T}, srcnum::JUDI.RangeOrVec) where T
  con = deepcopy(geometry.container)
  con.pkeyvals = con.pkeyvals[srcnum]
  return H5GeometryOOC(container=con, key=geometry.key,
                      xkey=geometry.xkey, ykey=geometry.ykey, zkey=geometry.zkey, 
                      do_coord_transform=geometry.do_coord_transform)
end

getindex(geometry::H5GeometryOOC{T}, srcnum::Integer) where T = getindex(geometry, srcnum:srcnum)

# Compare geometries
function compareGeometry(geometry_A::H5GeometryOOC, geometry_B::H5GeometryOOC)
  if geometry_A.key != geometry_B.key
    return false
  end

  if geometry_A.xkey != geometry_B.xkey || 
    geometry_A.ykey != geometry_B.ykey ||
    geometry_A.zkey != geometry_B.zkey
    return false
  end

  if !geometry_A.container.seis.isEqual(geometry_B.container.seis)
    return false
  end

  if geometry_A.container.pkey != geometry_B.container.pkey
    return false
  end

  if geometry_A.dt != geometry_B.dt
    return false
  end

  if geometry_A.nt != geometry_B.nt
    return false
  end

  if geometry_A.t != geometry_B.t
    return false
  end

  if geometry_A.nrec != geometry_B.nrec
    return false
  end

  return geometry_A.container.pkeyvals == geometry_B.container.pkeyvals
end

function compareGeometry(geometry_A::H5GeometryOOC, geometry_B::GeometryOOC)
  if geometry_A.key != geometry_B.key
    return false
  end

  if geometry_A.dt != geometry_B.dt
    return false
  end

  if geometry_A.nt != geometry_B.nt
    return false
  end

  if geometry_A.t != geometry_B.t
    return false
  end

  return geometry_A.nrec == geometry_B.nrec
end

compareGeometry(geometry_A::GeometryOOC, geometry_B::H5GeometryOOC) = compareGeometry(geometry_B, geometry_A)
compareGeometry(geometry_A::H5GeometryOOC, geometry_B::Geometry) = true
compareGeometry(geometry_A::Geometry, geometry_B::H5GeometryOOC) = true

######################## shapes easy access ################################
# Load geometry from out-of-core Geometry container
function Geometry(geometry::H5GeometryOOC)
  if isnothing(geometry.container.seis) || isnothing(geometry.container.seis) ||
    (geometry.key != "source" && geometry.key != "receiver")
    @error "Geometry key is Nothing or it is neither 'source' nor 'receiver'\n"
    return
  end

  if isnothing(geometry.container.seis)
    @error "Seis is Nothing\n"
    return
  end

  if isnothing(geometry.xkey) || 
    isnothing(geometry.ykey) || 
    isnothing(geometry.zkey)
    @error "xkey/ykey/zkey is Nothing\n"
    return
  end

  if isnothing(geometry.container.pkeyvals)
    @error "pkeyvals is Nothing\n"
    return
  end

  nsrc = length(geometry.container)
  xCell = Vector{Vector{Float32}}(undef, nsrc)
  yCell = Vector{Vector{Float32}}(undef, nsrc)
  zCell = Vector{Vector{Float32}}(undef, nsrc)
  dtCell = Vector{Float32}(undef, nsrc)
  ntCell = Vector{Integer}(undef, nsrc)
  tCell = Vector{Float32}(undef, nsrc)
  keylist = [geometry.container.pkey]
  for block in 1:nsrc
    minlist = [geometry.container.pkeyvals[block]]
    maxlist = [geometry.container.pkeyvals[block]]
    _, _, ind = geometry.container.seis.getSortedData(keylist, minlist, maxlist, 0, 0)
    if length(ind) < 1
      continue
    end
    if geometry.key == "source"
      ind = [ind[1]]
    end
    xy = Float32.(geometry.container.seis.getXYTraceHeaders(
      [geometry.xkey, geometry.ykey], ind, "m", geometry.do_coord_transform))

    xy[:,1] = xy[:,1] .- model_origin_x
    xy[:,2] = xy[:,2] .- model_origin_y
    xy = rotate_Nx2_array(xy, -model_orientation)
    xy[:,1] = xy[:,1] .+ model_origin_x
    xy[:,2] = xy[:,2] .+ model_origin_y

    xCell[block] = xy[:,1]
    yCell[block] = xy[:,2]
    zCell[block] = Float32.(geometry.container.seis.getTraceHeader(
      [geometry.zkey], ind, [geometry.container.seis.getLengthUnits()], ["m"])[:])
    dtCell[block] = geometry.dt[block]
    ntCell[block] = geometry.nt[block]
    tCell[block] = geometry.t[block]
  end

  return GeometryIC{Float32}(xCell, yCell, zCell, dtCell, ntCell, tCell)
end

function Geometry(;
  container::H5SeisCon,
  xkey::String, 
  ykey::String, 
  zkey::String, 
  do_coord_transform::Bool)

  if isnothing(container.seis)
    @error "Seis is Nothing\n"
    return
  end

  if isnothing(container.pkey) || isnothing(xkey) || isnothing(ykey) || isnothing(zkey)
    @error "pkey/xkey/ykey/zkey is Nothing\n"
    return
  end

  if isnothing(container.pkeyvals)
    @error "pkeyvals is Nothing\n"
    return
  end

  dt = Float32(abs(container.seis.getSampRate("ms")))
  nt = Integer(container.seis.getNSamp())
  t = Float32((nt-1)*dt)

  nsrc = length(container.pkeyvals)
  xCell = Vector{Vector{Float32}}(undef, nsrc)
  yCell = Vector{Vector{Float32}}(undef, nsrc)
  zCell = Vector{Vector{Float32}}(undef, nsrc)
  dtCell = Vector{Float32}(undef, nsrc)
  ntCell = Vector{Integer}(undef, nsrc)
  tCell = Vector{Float32}(undef, nsrc)
  keylist = [container.pkey]
  for block in 1:nsrc
    minlist = [container.pkeyvals[block]]
    maxlist = [container.pkeyvals[block]]
    _, _, ind = container.seis.getSortedData(keylist, minlist, maxlist, 0, 0)
    xy = Float32.(container.seis.getXYTraceHeaders([xkey, ykey], ind, "m", do_coord_transform))

    xy[:,1] = xy[:,1] .- model_origin_x
    xy[:,2] = xy[:,2] .- model_origin_y
    xy = rotate_Nx2_array(xy, -model_orientation)
    xy[:,1] = xy[:,1] .+ model_origin_x
    xy[:,2] = xy[:,2] .+ model_origin_y

    xCell[block] = xy[:,1]
    yCell[block] = xy[:,2]
    zCell[block] = Float32.(container.seis.getTraceHeader([zkey], ind, [container.seis.getLengthUnits()], ["m"])[:])
    dtCell[block] = dt
    ntCell[block] = nt
    tCell[block] = t
  end

  return GeometryIC{Float32}(xCell, yCell, zCell, dtCell, ntCell, tCell)
end