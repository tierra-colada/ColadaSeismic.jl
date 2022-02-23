
# src_frq in KHz
function h5wavemodeling2d(
  h5model::PyCall.PyObject, h5geom::PyCall.PyObject;
  h5model_opt::H5ModelOptions=H5ModelOptions(),
  h5geom_opt::H5GeomOptions=H5GeomOptions(),
  opt::JUDI.Options=nothing,
  src_frq::Float32=0.01f0,
  src_dt::Number=2,
  src_nt::Integer=50,
  rec_dt::Number=2,
  rec_nt::Integer=50)

  # pylogging = pyimport("logging")
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  qt = nothing
  slicer = nothing
  try
    qt = pyimport("qt")
    slicer = pyimport("slicer")
  catch
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling"
    return
  end

  if src_frq < 0
    @error "Source frequency is negative"
    return
  end

  if src_dt <= 0 || src_nt < 1 || rec_dt <= 0 || rec_nt < 1
    @error "Source/Receiver timings negative"
    return
  end

  model = prepare_JUDI_model_2D(h5model, h5opt=h5model_opt)
  if isnothing(model)
    @error "Unable to prepare JUDI Model"
    return
  end

  usrc_x = h5geom.getPKeyValues(h5geom_opt.src_xkey, h5geom.getLengthUnits(), "m")
  if isempty(usrc_x)
    @error "No $h5geom_opt.src_xkey sorting"
    return
  end

  # receiver sampling and recording time (JUDI works with positive values)
  src_t1 = (src_nt-1)*src_dt

  # setup wavelet
  wavelet = ricker_wavelet(src_t1, src_dt, src_frq)

  nsrc = length(usrc_x)
  progressDialog = nothing
  if !isnothing(qt)
    progressDialog = qt.QProgressDialog("Forward modeling", "Abort", 1, nsrc+1)
  end

  for i in 1:nsrc
    if !isnothing(qt)
      progressDialog.setValue(i)
      slicer.app.processEvents()
    end

    h5geom_opt.src_xkey_min = usrc_x[i]
    h5geom_opt.src_xkey_max = usrc_x[i]
    srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_2D(
      h5geom, h5opt=h5geom_opt,
      src_dt=src_dt,
      src_nt=src_nt,
      rec_dt=rec_dt,
      rec_nt=rec_nt)

    if isnothing(srcGeometry) || isnothing(recGeometry) || isnothing(indCellVec)
      @error "Unable to prepare Source and Receiver JUDI Geometry"
      return
    end

    q = judiVector(srcGeometry, wavelet)

    # Set up info structure for linear operators
    ntComp = get_computational_nt(srcGeometry, recGeometry, model)
    info = Info(prod(model.n), 1, ntComp)

    ###################################################################################################

    # Setup operators
    Pr = judiProjection(info, recGeometry)
    F = judiModeling(info, model; options=opt)
    Ps = judiProjection(info, srcGeometry)

    # Nonlinear modeling
    dobs = Pr*F*adjoint(Ps)*q
    # qad = Ps*adjoint(F)*adjoint(Pr)*dobs

    for ii in 1:length(indCellVec[1])
      h5geom.writeTrace(dobs.data[1][:,ii], indCellVec[1][ii])
    end
    h5geom.getH5File().flush(true)

    if !isnothing(qt) && progressDialog.wasCanceled
      break
    end
  end

  if !isnothing(qt)
    progressDialog.setValue(nsrc+1)
    slicer.app.processEvents()
  end
end



# src_frq in KHz
function h5wavemodeling3d(
  h5model::PyCall.PyObject, h5geom::PyCall.PyObject;
  h5model_opt::H5ModelOptions=H5ModelOptions(),
  h5geom_opt::H5GeomOptions=H5GeomOptions(),
  opt::JUDI.Options=nothing,
  src_frq::Float32=0.01f0,
  src_dt::Number=2,
  src_nt::Integer=50,
  rec_dt::Number=2,
  rec_nt::Integer=50)

  # pylogging = pyimport("logging")
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  qt = nothing
  slicer = nothing
  try
    qt = pyimport("qt")
    slicer = pyimport("slicer")
  catch
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling"
    return
  end

  if src_frq < 0
    @error "Source frequency is negative"
    return
  end

  if src_dt <= 0 || src_nt < 1 || rec_dt <= 0 || rec_nt < 1
    @error "Source/Receiver timings negative"
    return
  end

  model, h5geom_opt.model_orientation = prepare_JUDI_model_3D(h5model, h5opt=h5model_opt)
  if isnothing(model)
    @error "Unable to prepare JUDI Model"
    return
  end

  keylist = [h5geom_opt.src_xkey, h5geom_opt.src_ykey]
  minlist = [-Inf, -Inf]
  maxlist = [Inf, Inf]
  ~, src_xy, ~ = h5geom.getSortedData(keylist, minlist, maxlist, 0, 0, true, "", "m", true)
  if isempty(src_xy)
    @error "Unable to get sorted src_x src_y headers. Probably $h5geom_opt.src_xkey is missing"
    return
  end

  ~, usrc_xy, ~ = h5geo.sort_rows_unique(src_xy)

  # receiver sampling and recording time (JUDI works with positive values)
  src_t1 = (src_nt-1)*src_dt

  # setup wavelet
  wavelet = ricker_wavelet(src_t1, src_dt, src_frq)

  nsrc = size(usrc_xy)[1]
  progressDialog = nothing
  if !isnothing(qt)
    progressDialog = qt.QProgressDialog("Forward modeling", "Abort", 1, nsrc+1)
  end

  for i in 1:nsrc
    if !isnothing(qt)
      progressDialog.setValue(i)
      slicer.app.processEvents()
    end

    h5geom_opt.src_xkey_min = usrc_x[i]
    h5geom_opt.src_xkey_max = usrc_x[i]
    h5geom_opt.src_ykey_min = usrc_y[i]
    h5geom_opt.src_ykey_max = usrc_y[i]
    srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_3D(
      h5geom, h5opt=h5geom_opt,
      src_dt=src_dt,
      src_nt=src_nt,
      rec_dt=rec_dt,
      rec_nt=rec_nt)
    if isnothing(srcGeometry) || isnothing(recGeometry) || isnothing(indCellVec)
      @error "Unable to prepare Source and Receiver JUDI Geometry"
      return
    end

    q = judiVector(srcGeometry, wavelet)

    # Set up info structure for linear operators
    ntComp = get_computational_nt(srcGeometry, recGeometry, model)
    info = Info(prod(model.n), 1, ntComp)

    ###################################################################################################

    # Setup operators
    Pr = judiProjection(info, recGeometry)
    F = judiModeling(info, model; options=opt)
    Ps = judiProjection(info, srcGeometry)

    # Nonlinear modeling
    dobs = Pr*F*adjoint(Ps)*q
    # qad = Ps*adjoint(F)*adjoint(Pr)*dobs

    for ii in 1:length(indCellVec[1])
      h5geom.writeTrace(dobs.data[1][:,ii], indCellVec[1][ii])
    end
    h5geom.getH5File().flush(true)

    if !isnothing(qt) && progressDialog.wasCanceled
      break
    end
  end

  if !isnothing(qt)
    progressDialog.setValue(nsrc+1)
    slicer.app.processEvents()
  end
end



# src_frq in KHz
function h5wavemodeling3d_segy(
  h5model::PyCall.PyObject, h5geom::PyCall.PyObject;
  h5model_opt::H5ModelOptions=H5ModelOptions(),
  h5geom_opt::H5GeomOptions=H5GeomOptions(),
  opt::JUDI.Options=nothing,
  src_frq::Float32=0.01f0,
  src_dt::Number=2,
  src_nt::Integer=50,
  rec_dt::Number=2,
  rec_nt::Integer=50)

  # pylogging = pyimport("logging")
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  qt = nothing
  slicer = nothing
  try
    qt = pyimport("qt")
    slicer = pyimport("slicer")
  catch
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling"
    return
  end

  if src_frq < 0
    @error "Source frequency is negative"
    return
  end

  if src_dt <= 0 || src_nt < 1 || rec_dt <= 0 || rec_nt < 1
    @error "Source/Receiver timings negative"
    return
  end

  model, h5geom_opt.model_orientation = prepare_JUDI_model_3D(h5model, h5opt=h5model_opt)
  if isnothing(model)
    @error "Unable to prepare JUDI Model"
    return
  end

  # receiver sampling and recording time (JUDI works with positive values)
  src_t1 = (src_nt-1)*src_dt

  # setup wavelet
  wavelet = ricker_wavelet(src_t1, src_dt, src_frq)

  srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_3D(
    h5geom, h5opt=h5geom_opt,
    src_dt=src_dt,
    src_nt=src_nt,
    rec_dt=rec_dt,
    rec_nt=rec_nt)
  if isnothing(srcGeometry) || isnothing(recGeometry) || isnothing(indCellVec)
    @error "Unable to prepare Source and Receiver JUDI Geometry"
    return
  end

  q = judiVector(srcGeometry, wavelet)

  # Set up info structure for linear operators
  ntComp = get_computational_nt(srcGeometry, recGeometry, model)
  info = Info(prod(model.n), length(srcGeometry.xloc), ntComp)

  ###################################################################################################

  opt.save_data_to_disk = true
  opt.file_path = "/home/kerim/Documents/Colada_prj/my JUDI test"
  opt.file_name = "test1"

  # Setup operators
  Pr = judiProjection(info, recGeometry)
  F = judiModeling(info, model; options=opt)
  Ps = judiProjection(info, srcGeometry)

  # Nonlinear modeling
  dobs = Pr*F*adjoint(Ps)*q
  # qad = Ps*adjoint(F)*adjoint(Pr)*dobs
    
end