function H5ForwardModeling(;
  model::JUDI.Model,
  opt::JUDI.JUDIOptions,
  q::JUDI.judiVector,
  recGeometry)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(recGeometry)
  F = judiModeling(model; options=opt)
  Ps = judiProjection(q.geometry)

  dobs = Pr*F*adjoint(Ps)*q
  @info "Modeling finished"
end