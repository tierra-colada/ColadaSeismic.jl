function H5ForwardModeling(;
  model::JUDI.Model,
  opt::JUDI.Options,
  q::JUDI.judiVector,
  recGeometry)

  nsrc = length(q.geometry.container)

  # Set up info structure for linear operators
  ntComp = get_computational_nt(q.geometry, recGeometry, model)
  info = Info(prod(model.n), nsrc, ntComp)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(info, recGeometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, q.geometry)

  # Nonlinear modeling
  dobs = Pr*F*adjoint(Ps)*q
end