function H5ForwardModeling(;
  model::JUDI.Model,
  opt::JUDI.Options,
  q::JUDI.judiVector,
  info::JUDI.Info,
  recGeometry)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(info, recGeometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, q.geometry)
  
  dobs = Pr*F*adjoint(Ps)*q
end