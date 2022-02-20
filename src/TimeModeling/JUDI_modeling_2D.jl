include("prepare_JUDI_model.jl")
include("prepare_JUDI_geometry.jl")

# frq in KHz
function h5wavemodeling2d(
  h5model::PyCall.PyObject, h5geom::PyCall.PyObject;
  opt::JUDI.Options=nothing,
  model_xkey::String="CDP_X",
  src_xkey::String="SRCX",
  src_zkey::String="SES",
  rec_xkey::String="GRPX",
  rec_zkey::String="RGE",
  frq::Float32=0.01f0)

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

  model = prepare_JUDI_model_2D(h5model, model_xkey=model_xkey)
  if isnothing(model)
    @error "Unable to prepare JUDI Model"
    return
  end

  usrc_x = h5geom.getPKeyValues(src_xkey, h5geom.getLengthUnits(), "m")
  if isempty(usrc_x)
    @error "No $src_xkey sorting"
    return
  end

  # receiver sampling and recording time
  timeR = h5geom.getLastSample(0, "ms") * (-1)   # receiver recording time [ms]
  dtR = h5geom.getSampRate("ms") * (-1)    # receiver sampling interval [ms]

  # source sampling and number of time steps
  timeS = timeR  # ms
  dtS = dtR   # ms

  # setup wavelet
  wavelet = ricker_wavelet(timeS, dtS, frq)

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

    srcGeometry, recGeometry, indCellVec = prepare_JUDI_geometry_2D(
      h5geom,
      src_xkey_min=usrc_x[i],
      src_xkey_max=usrc_x[i],
      src_xkey=src_xkey,
      src_zkey=src_zkey,
      rec_xkey=rec_xkey,
      rec_zkey=rec_zkey)

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