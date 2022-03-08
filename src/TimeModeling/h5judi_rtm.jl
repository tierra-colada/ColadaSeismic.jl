function H5RTM(;
  model::JUDI.Model,
  opt::JUDI.Options,
  src_frq::Float32,
  srcGeometry, 
  dobs::JUDI.judiVector)

  # setup wavelet
  wavelet = ricker_wavelet(srcGeometry.t[1], srcGeometry.dt[1], src_frq)
  q = judiVector(srcGeometry, wavelet)

  nsrc = length(srcGeometry.xloc)

  # Set up info structure for linear operators
  ntComp = get_computational_nt(srcGeometry, dobs.geometry, model)
  info = Info(prod(model.n), nsrc, ntComp)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(info, dobs.geometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, q.geometry)
  J = judiJacobian(Pr*F0*adjoint(Ps), q)

  # RTM
  rtm = adjoint(J)*dobs
end