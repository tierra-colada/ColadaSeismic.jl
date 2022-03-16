import JUDI.get_data, JUDI.judiVector, JUDI.convert


function judiVector(con::H5SeisCon)
  nSamp = con.seis.getNSamp()
  nTrc = 0
  for val in con.pkeyvals
    nTrc += con.seis.getPKeyTraceSize(con.pkey, val, val)
  end
  m = nSamp*nTrc
  n = 1
  nsrc = length(con.pkeyvals)

  geometry = Geometry(con)
  dataCell = Array{H5SeisCon,1}(undef, nsrc)
  for i = 1:nsrc
    pkeyvals = Array{Float64,1}(undef,1)
    pkeyvals[1] = con.pkeyvals[i]
    dataCell[i] = H5SeisCon(seis=con.seis, 
                            xkey=con.xkey, ykey=con.ykey, zkey=con.zkey, 
                            pkey=con.pkey, pkeyvals=pkeyvals)
  end

  return judiVector{Float32, H5SeisCon}("Julia seismic data container",m,n,nsrc,geometry,dataCell)
end


function judiVector(geometry::JUDI.Geometry, con::H5SeisCon)
  nSamp = con.seis.getNSamp()
  nTrc = 0
  for val in con.pkeyvals
    nTrc += con.seis.getPKeyTraceSize(con.pkey, val, val)
  end
  m = nSamp*nTrc
  n = 1
  nsrc = length(con.pkeyvals)

  # fill data vector with pointers to data location
  dataCell = Array{H5SeisCon,1}(undef, nsrc)
  for j=1:nsrc
      dataCell[j] = split(con,j)
  end

  return judiVector{Float32, H5SeisCon}("Julia seismic data container",m,n,nsrc,geometry,dataCell)
end


function judiVector(geometry::JUDI.Geometry, conCell::Array{H5SeisCon})
  nSamp = conCell[1].seis.getNSamp()
  nTrc = 0
  nsrc = 0
  for con in conCell
    for val in con.pkeyvals
      nsrc += length(con.pkeyvals)
      nTrc += con.seis.getPKeyTraceSize(con.pkey, val, val)
    end
  end
  m = nSamp*nTrc
  n = 1

  # fill data vector with pointers to data location
  dataCell = Array{H5SeisCon,1}(undef, nsrc)
  for j=1:nsrc
      dataCell[j] = conCell[j]
  end

  return judiVector{Float32, H5SeisCon}("Julia seismic data container",m,n,nsrc,geometry,dataCell)
end

##########################################################

vec(x::H5SeisCon) = vec(x[1].data)

##########################################################

for opo=[:+, :-, :*, :/]
  @eval begin
    $opo(a::judiVector{avDT, H5SeisCon}, ::T) where {avDT, T<:Real} = throw(DomainError(a, "Addition for OOC judiVectors not supported."))
    $opo(::T, b::judiVector{bvDT, H5SeisCon}) where {bvDT, T<:Real} = throw(DomainError(b, "Addition for OOC judiVectors not supported."))
    $opo(::judiVector{avDT, H5SeisCon}, b::judiVector{bvDT, H5SeisCon}) where {avDT, bvDT} = throw(DomainError(b, "Addition for OOC judiVectors not supported."))
  end
end

for ipop=[:lmul!, :rmul!, :rdiv!, :ldiv!]
  @eval begin
    $ipop(a::judiVector{avDT, H5SeisCon}, ::T) where {avDT, T<:Real} = throw(DomainError(a, "Addition for OOC judiVectors not supported."))
    $ipop(::T, b::judiVector{bvDT, H5SeisCon}) where {bvDT, T<:Real} = throw(DomainError(b, "Addition for OOC judiVectors not supported."))
    $ipop(::judiVector{avDT, H5SeisCon}, b::judiVector{bvDT, H5SeisCon}) where {avDT, bvDT} = throw(DomainError(b, "Addition for OOC judiVectors not supported."))
  end
end

function get_data(x::judiVector{T, H5SeisCon}) where T
  shots = Array{Array{Float32, 2}, 1}(undef, x.nsrc)
  rec_geometry = Geometry(x.geometry)
  for j=1:x.nsrc
    shots[j] = convert(Array{Float32, 2}, x.data[j][1].data)
  end
  return judiVector(rec_geometry, shots)
end

##########################################################
convert(::Type{Matrix{T}}, x::H5SeisCon) where T = convert(Matrix{T}, x[1].data)