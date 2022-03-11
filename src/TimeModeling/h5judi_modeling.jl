# operation_type: Forward Modeling, RTM, FWI, TWRI
function H5Modeling(;
  # models
  h5vel::PyCall.PyObject,
  h5density::PyCall.PyObject=nothing,
  h5epsilon::PyCall.PyObject=nothing,
  h5delta::PyCall.PyObject=nothing,
  h5tetha::PyCall.PyObject=nothing,
  h5phi::PyCall.PyObject=nothing,
  model_xkey::String="CDP_X",
  model_ykey::String="CDP_Y",

  # geometry
  h5geom::PyCall.PyObject=nothing,
  geom_segy_files::Array{String,1}=nothing,
  geom_length_units::String,
  geom_temporal_units::String,
  geom_crs::String,
  geom_src_xkey::String,
  geom_src_ykey::String,
  geom_src_zkey::String,
  geom_rec_xkey::String,
  geom_rec_ykey::String,
  geom_rec_zkey::String,

  # judi
  opt::JUDI.Options,
  src_frq::Float32=0.01f0,

  # FWI constraints
  vmin::Number,
  vmax::Number,

  computation_type::String,
  fwi_niter::Int=nothing,
  fwi_batchsize::Int=nothing,
  twri_opt::TWRIOptions=nothing
  )

  # pylogging = pyimport("logging")
  h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")

  if isnothing(h5vel)
    @error "No velocity model"
    return
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling"
    return
  end

  phpvel, model_orientation = H5ReadPhysicalParameter3D(
    h5vel, phptype=VELOCITY, xkey=model_xkey, ykey=model_ykey)
  if isnothing(phpvel)
    @error "Unable to read PhysicalParameter: VELOCITY"
    return
  end
  model = JUDI.Model(phpvel.n, phpvel.d, phpvel.o, phpvel.data)

  if !isnothing(h5density)
    phpdensity, model_orientation = H5ReadPhysicalParameter3D(
      h5density, phptype=DENSITY, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdensity)
      @error "Unable to read PhysicalParameter: DENSITY"
      return
    end
  end

  if !isnothing(h5epsilon)
    phpepsilon, model_orientation = H5ReadPhysicalParameter3D(
      h5epsilon, phptype=EPSILON, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpepsilon)
      @error "Unable to read PhysicalParameter: EPSILON"
      return
    end
  end

  if !isnothing(h5delta)
    phpdelta, model_orientation = H5ReadPhysicalParameter3D(
      h5delta, phptype=DELTA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdelta)
      @error "Unable to read PhysicalParameter: DELTA"
      return
    end
  end

  if !isnothing(h5tetha)
    phptheta, model_orientation = H5ReadPhysicalParameter3D(
      h5tetha, phptype=THETA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phptheta)
      @error "Unable to read PhysicalParameter: THETA"
      return
    end
  end

  if !isnothing(h5phi)
    phpphi, model_orientation = H5ReadPhysicalParameter3D(
      h5phi, phptype=PHI, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpphi)
      @error "Unable to read PhysicalParameter: PHI"
      return
    end
  end

  model = JUDI.Model(phpvel.n, phpvel.d, phpvel.o, phpvel.data, 
          rho=phpdensity.data, epsilon=phpepsilon.data, delta=phpdelta.data, 
          theta=phptheta.data, phi=phpphi.data)

  # Geometry
  if isnothing(h5geom) && isnothing(geom_segy_files)
    @error "Neither geometry H5Seis nor SEGY files were provided"
    return
  end

  if !isnothing(geom_segy_files)
    geom_rec_xkey_judi = h5geo2judiTraceHeaderName(geom_rec_xkey)
    if isnothing(geom_rec_xkey_judi)
      @error "Unable to map h5geo trace header: $geom_rec_xkey to judi"
      return
    end
    geom_rec_ykey_judi = h5geo2judiTraceHeaderName(geom_rec_ykey)
    if isnothing(geom_rec_ykey_judi)
      @error "Unable to map h5geo trace header: $geom_rec_ykey to judi"
      return
    end
    geom_rec_zkey_judi = h5geo2judiTraceHeaderName(geom_rec_zkey)
    if isnothing(geom_rec_zkey_judi)
      @error "Unable to map h5geo trace header: $geom_rec_zkey to judi"
      return
    end
    geom_src_zkey_judi = h5geo2judiTraceHeaderName(geom_src_zkey)
    if isnothing(geom_src_zkey_judi)
      @error "Unable to map h5geo trace header: $geom_src_zkey to judi"
      return
    end
    # before using there mapped keys I need to modify JUDI Geometry.jl to allow to use custom headers
    container = segy_scan(geom_segy_files, ["GroupX", "GroupY", "RecGroupElevation", "SourceSurfaceElevation", "dt"])
    d_obs = judiVector(container)
    d_obs.geometry.h5geo = h5geo
    d_obs.geometry.crs = geom_crs
    d_obs.geometry.lengthUnits = geom_length_units
    d_obs.geometry.temporalUnits = geom_temporal_units
    d_obs.geometry.model_origin_x = model.o[1]
    d_obs.geometry.model_origin_y = model.o[2]
    d_obs.geometry.model_orientation = model_orientation
    recGeometry = d_obs.geometry

    # set up source
    srcGeometry = Geometry(container; key = "source")
    srcGeometry.h5geo = h5geo
    srcGeometry.crs = geom_crs
    srcGeometry.lengthUnits = geom_length_units
    srcGeometry.temporalUnits = geom_temporal_units
    srcGeometry.model_origin_x = model.o[1]
    srcGeometry.model_origin_y = model.o[2]
    srcGeometry.model_orientation = model_orientation
  else
    dt = h5geom.getSampRate("ms")
    nt = h5geom.getNSamp()
    srcGeometry, recGeometry, indCellVec = H5ReadGeometry(
      h5geom,
      src_xkey=geom_src_xkey,
      src_ykey=geom_src_ykey,
      src_zkey=geom_src_zkey,
      rec_xkey=geom_rec_xkey,
      rec_ykey=geom_rec_ykey,
      rec_zkey=geom_rec_zkey,
      src_dt=dt,
      src_nt=nt,
      rec_dt=dt,
      rec_nt=nt)
    if isnothing(srcGeometry) || isnothing(srcGeometry) || isnothing(indCellVec)
      @error "Unable to prepare Source and Receiver JUDI Geometry objects from H5Seis"
      return
    end
  end

  # setup wavelet
  wavelet = ricker_wavelet(srcGeometry.t[1], srcGeometry.dt[1], src_frq)
  q = judiVector(srcGeometry, wavelet)

  if computation_type == "Forward Modeling"
    H5ForwardModeling(model=model, opt=opt, q=q, recGeometry=recGeometry)
  elseif computation_type == "RTM"
    H5RTM(model=model, opt=opt, q=q, dobs=dobs)
  elseif computation_type == "FWI"
    H5FWI(model=model, opt=opt, q=q, dobs=dobs, 
          niterations=fwi_niter, batchsize=fwi_batchsize, 
          vmin=vmin, vmax=vmax)
  elseif computation_type == "TWRI"
    @error "TWRI not ready yet"
  end

end