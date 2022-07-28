function H5RTM(;
  model::JUDI.Model,
  opt::JUDI.JUDIOptions,
  q::JUDI.judiVector, 
  dobs::JUDI.judiVector)

  ###################################################################################################

  # Setup operators
  Pr = judiProjection(dobs.geometry)
  F = judiModeling(model; options=opt)
  Ps = judiProjection(q.geometry)
  J = judiJacobian(Pr*F*adjoint(Ps), q)

  # RTM
  rtm = adjoint(J)*dobs

  H5WritePhysicalParameter(
    cntName="$(opt.file_path)/rtm.h5",
    objName="rtm",
    cntCreationType=h5geo.CreationType.OPEN_OR_CREATE,
    objCreationType=h5geo.CreationType.CREATE_UNDER_NEW_NAME,
    php=rtm
  )

  @info "RTM finished\n"
end