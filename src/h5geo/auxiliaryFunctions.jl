import JUDI.get_computational_nt

function get_computational_nt(srcGeometry::H5GeometryOOC, recGeometry::H5GeometryOOC, model::Model; dt=nothing)
  # Determine number of computational time steps
  nsrc = length(srcGeometry.container)
  nt = Array{Integer}(undef, nsrc)
  dtComp = calculate_dt(model; dt=dt)
  dt_src = abs(srcGeometry.container.seis.getSampRate("ms"))
  nt_src = srcGeometry.container.seis.getNSamp()
  dt_rec = abs(recGeometry.container.seis.getSampRate("ms"))
  nt_rec = recGeometry.container.seis.getNSamp()
  for j=1:nsrc
    ntRec = Int(dt_rec*(nt_rec-1) ÷ dtComp) + 1
    ntSrc = Int(dt_src*(nt_src-1) ÷ dtComp) + 1
    nt[j] = max(ntRec, ntSrc)
  end
  return nt
end
