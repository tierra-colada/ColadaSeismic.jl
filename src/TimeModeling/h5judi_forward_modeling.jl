function H5ForwardModeling(;
  model::JUDI.Model,
  opt::JUDI.Options,
  src_frq::Float32,
  srcGeometry, recGeometry)

  # setup wavelet
  wavelet = ricker_wavelet(srcGeometry.t[1], srcGeometry.dt[1], src_frq)
  q = judiVector(srcGeometry, wavelet)

  nsrc = length(srcGeometry.xloc)

  # Set up info structure for linear operators
  ntComp = get_computational_nt(srcGeometry, recGeometry, model)
  info = Info(prod(model.n), nsrc, ntComp)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(info, recGeometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, srcGeometry)

  # Nonlinear modeling
  dobs = Pr*F*adjoint(Ps)*q
end