function H5RTM(;
  model::JUDI.Model,
  opt::JUDI.Options,
  q::JUDI.judiVector, 
  dobs::JUDI.judiVector)

  nsrc = length(q.geometry.container)

  # Set up info structure for linear operators
  ntComp = get_computational_nt(q.geometry, dobs.geometry, model)
  info = Info(prod(model.n), nsrc, ntComp)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(info, dobs.geometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, q.geometry)
  J = judiJacobian(Pr*F*adjoint(Ps), q)

  # RTM
  rtm = adjoint(J)*dobs
end