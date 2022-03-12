import JUDI.Geometry


function Geometry(con::H5SeisCon)
  if isnothing(con.seis)
    @error "Seis is Nothing"
    return
  end

  if isnothing(con.pkey) || isnothing(con.xkey) || isnothing(con.ykey) || isnothing(con.zkey)
    @error "pkey/xkey/ykey/zkey is Nothing"
    return
  end

  if isnothing(con.pkeyvals)
    @error "pkeyvals is Nothing"
    return
  end

  dt = Float32(abs(con.seis.getSampRate("ms")))
  nt = Integer(con.seis.getNSamp())
  t = Float32((nt-1)*dt)

  nsrc = length(con.pkeyvals)
  xCell = Vector{Vector{Float32}}(undef, nsrc)
  yCell = Vector{Vector{Float32}}(undef, nsrc)
  zCell = Vector{Vector{Float32}}(undef, nsrc)
  dtCell = Vector{Float32}(undef, nsrc)
  ntCell = Vector{Integer}(undef, nsrc)
  tCell = Vector{Float32}(undef, nsrc)
  keylist = [con.pkey]
  for block in 1:nsrc
    minlist = [con.pkeyvals[block]]
    maxlist = [con.pkeyvals[block]]
    _, _, ind = con.seis.getSortedData(keylist, minlist, maxlist, 0, 0)
    xy = Float32.(con.seis.getXYTraceHeaders([con.xkey, con.ykey], ind, "m", false))  # MUST BE TRUE

    xCell[block] = xy[:,1]
    yCell[block] = xy[:,2]
    zCell[block] = Float32.(con.seis.getTraceHeader([con.zkey], ind, [con.seis.getLengthUnits()], ["m"])[:])
    dtCell[block] = dt
    ntCell[block] = nt
    tCell[block] = t
  end

  return GeometryIC{Float32}(xCell, yCell, zCell, dtCell, ntCell, tCell)
end