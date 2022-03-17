import JUDI.Geometry, JUDI.get_nsrc, JUDI.n_samples, JUDI.subsample, JUDI.compareGeometry, JUDI.isequal

export H5GeometryOOC

mutable struct H5GeometryOOC <: Geometry{Float32}
	h5geo::PyCall.PyObject
	container::H5SeisCon
	key::String
	xkey::String
	ykey::String
	zkey::String
	model_orientation::Number
end

function H5GeometryOOC(;
  h5geo::PyCall.PyObject,
  container::H5SeisCon,
  key::String,
  xkey::String,
  ykey::String,
  zkey::String,
  model_orientation::Number)
  if isnothing(container)
    @error "H5SeisCon is Nothing"
    return
  end

  if isnothing(key)
    @error "key is Nothing"
    return
  end

  if key != "source" && key != "receiver"
    @error "key must be either 'source' or 'receiver'"
    return
  end

  if isnothing(xkey) || isnothing(ykey) || isnothing(zkey)
    @error "xkey/ykey/zkey is Nothing"
    return
  end

  return H5GeometryOOC(h5geo, container, key, xkey, ykey, zkey, model_orientation)
end

######################## shapes easy access ################################
get_nsrc(g::H5GeometryOOC) = length(g.container.pkeyvals)

function n_samples(g::H5GeometryOOC, info::Info)
  return n_samples(g, info.nsrc)
end

# if work slow then I need to add 'nsamp' as a member var of H5GeometryOOC
function n_samples(g::H5GeometryOOC, nsrc::Integer)
  nt = Integer(g.container.seis.getNSamp())
  nSamp = 0
  for j=1:nsrc
    nSamp += g.container.seis.getPKeyTraceSize(g.container.pkey, g.container.pkeyvals[j], g.container.pkeyvals[j])
  end
  return nSamp*nt
end

# Subsample out-of-core geometry structure
function subsample(geometry::H5GeometryOOC, srcnum)
  con = deepcopy(geometry.container)
  con.pkeyvals = [con.pkeyvals[srcnum]]
  return H5GeometryOOC(h5geo=geometry.h5geo, container=con, key=geometry.key,
                      xkey=geometry.xkey, ykey=geometry.ykey, zkey=geometry.zkey, 
                      model_orientation=geometry.model_orientation)
end

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

  return geometry_A.container.pkeyvals == geometry_B.container.pkeyvals
end

isequal(geometry_A::H5GeometryOOC, geometry_B::H5GeometryOOC) = compareGeometry(geometry_A, geometry_B)

compareGeometry(geometry_A::H5GeometryOOC, geometry_B::Geometry) = true
compareGeometry(geometry_A::Geometry, geometry_B::H5GeometryOOC) = true

######################## shapes easy access ################################
# Load geometry from out-of-core Geometry container
function Geometry(geometry::H5GeometryOOC)
  if isnothing(geometry.container.seis) || isnothing(geometry.container.seis) ||
    (geometry.key != "source" && geometry.key != "receiver")
    @error "Geometry key is Nothing or it is neither 'source' nor 'receiver'"
    return
  end

  if isnothing(geometry.container.seis)
    @error "Seis is Nothing"
    return
  end

  if isnothing(geometry.xkey) || 
    isnothing(geometry.ykey) || 
    isnothing(geometry.zkey)
    @error "xkey/ykey/zkey is Nothing"
    return
  end

  if isnothing(geometry.container.pkeyvals)
    @error "pkeyvals is Nothing"
    return
  end

  dt = Float32(abs(geometry.container.seis.getSampRate("ms")))
  nt = Integer(geometry.container.seis.getNSamp())
  t = Float32((nt-1)*dt)

  nsrc = length(geometry.container.pkeyvals)
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
      [geometry.xkey, geometry.ykey], ind, "m", false))  # MUST BE TRUE

    xCell[block] = xy[:,1]
    yCell[block] = xy[:,2]
    zCell[block] = Float32.(geometry.container.seis.getTraceHeader(
      [geometry.zkey], ind, [geometry.container.seis.getLengthUnits()], ["m"])[:])
    dtCell[block] = dt
    ntCell[block] = nt
    tCell[block] = t
  end

  return GeometryIC{Float32}(xCell, yCell, zCell, dtCell, ntCell, tCell)
end

# function Geometry(con::H5SeisCon)
#   if isnothing(con.seis)
#     @error "Seis is Nothing"
#     return
#   end

#   if isnothing(con.pkey) || isnothing(con.xkey) || isnothing(con.ykey) || isnothing(con.zkey)
#     @error "pkey/xkey/ykey/zkey is Nothing"
#     return
#   end

#   if isnothing(con.pkeyvals)
#     @error "pkeyvals is Nothing"
#     return
#   end

#   dt = Float32(abs(con.seis.getSampRate("ms")))
#   nt = Integer(con.seis.getNSamp())
#   t = Float32((nt-1)*dt)

#   nsrc = length(con.pkeyvals)
#   xCell = Vector{Vector{Float32}}(undef, nsrc)
#   yCell = Vector{Vector{Float32}}(undef, nsrc)
#   zCell = Vector{Vector{Float32}}(undef, nsrc)
#   dtCell = Vector{Float32}(undef, nsrc)
#   ntCell = Vector{Integer}(undef, nsrc)
#   tCell = Vector{Float32}(undef, nsrc)
#   keylist = [con.pkey]
#   for block in 1:nsrc
#     minlist = [con.pkeyvals[block]]
#     maxlist = [con.pkeyvals[block]]
#     _, _, ind = con.seis.getSortedData(keylist, minlist, maxlist, 0, 0)
#     xy = Float32.(con.seis.getXYTraceHeaders([con.xkey, con.ykey], ind, "m", false))  # MUST BE TRUE

#     xCell[block] = xy[:,1]
#     yCell[block] = xy[:,2]
#     zCell[block] = Float32.(con.seis.getTraceHeader([con.zkey], ind, [con.seis.getLengthUnits()], ["m"])[:])
#     dtCell[block] = dt
#     ntCell[block] = nt
#     tCell[block] = t
#   end

#   return GeometryIC{Float32}(xCell, yCell, zCell, dtCell, ntCell, tCell)
# end