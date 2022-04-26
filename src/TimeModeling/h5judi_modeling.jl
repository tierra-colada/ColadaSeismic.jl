# operation_type: Forward Modeling, RTM, FWI, TWRI
function H5Modeling(;
  # models
  h5vel,
  h5density=nothing,
  h5epsilon=nothing,
  h5delta=nothing,
  h5tetha=nothing,
  h5phi=nothing,
  model_xkey::String="CDP_X",
  model_ykey::String="CDP_Y",

  # geometry
  h5geom=nothing,
  geom_src_pkey::String="SP",
  geom_src_xkey::String="SRCX",
  geom_src_ykey::String="SRCY",
  geom_src_zkey::String="SES",
  geom_rec_xkey::String="GRPX",
  geom_rec_ykey::String="GRPY",
  geom_rec_zkey::String="RGE",

  # judi
  opt::JUDI.Options,
  src_frq::Number,

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

  @info "H5Modeling started\n"

  if isnothing(h5vel)
    @error "No velocity model\n"
    return
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling\n"
    return
  end

  phpvel, model_orientation = H5ReadPhysicalParameter(
    h5vel, phptype=VELOCITY, xkey=model_xkey, ykey=model_ykey)
  if isnothing(phpvel)
    @error "Unable to read PhysicalParameter: VELOCITY\n"
    return
  end
  model = JUDI.Model(phpvel.n, phpvel.d, phpvel.o, phpvel.data)

  # declare PHP vars to be able use thm as arguments (or unknown var error will appear)
  phpdensity = nothing
  phpepsilon = nothing
  phpdelta = nothing
  phptheta = nothing
  phpphi = nothing

  if !isnothing(h5density)
    phpdensity, model_orientation = H5ReadPhysicalParameter(
      h5density, phptype=DENSITY, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdensity)
      @error "Unable to read PhysicalParameter: DENSITY\n"
      return
    end
  end

  if !isnothing(h5epsilon)
    phpepsilon, model_orientation = H5ReadPhysicalParameter(
      h5epsilon, phptype=EPSILON, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpepsilon)
      @error "Unable to read PhysicalParameter: EPSILON\n"
      return
    end
  end

  if !isnothing(h5delta)
    phpdelta, model_orientation = H5ReadPhysicalParameter(
      h5delta, phptype=DELTA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdelta)
      @error "Unable to read PhysicalParameter: DELTA\n"
      return
    end
  end

  if !isnothing(h5tetha)
    phptheta, model_orientation = H5ReadPhysicalParameter(
      h5tetha, phptype=THETA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phptheta)
      @error "Unable to read PhysicalParameter: THETA\n"
      return
    end
  end

  if !isnothing(h5phi)
    phpphi, model_orientation = H5ReadPhysicalParameter(
      h5phi, phptype=PHI, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpphi)
      @error "Unable to read PhysicalParameter: PHI\n"
      return
    end
  end

  model = JUDI.Model(phpvel.n, phpvel.d, phpvel.o, phpvel.data, 
          rho=isnothing(phpdensity) ? nothing : phpdensity.data, 
          epsilon=isnothing(phpepsilon) ? nothing : phpepsilon.data, 
          delta=isnothing(phpdelta) ? nothing : phpdelta.data, 
          theta=isnothing(phptheta) ? nothing : phptheta.data, 
          phi=isnothing(phpphi) ? nothing : phpphi.data)

  # Geometry
  if isnothing(h5geom)
    @error "Geometry H5Seis is NULL\n"
    return
  end

  model_origin_x = model.o[1]
  model_origin_y = 0
  if length(model.o) == 3
    model_origin_y = model.o[2]
  end

  con = H5SeisCon(seis=h5geom, pkey=geom_src_pkey)
  recGeometry = H5GeometryOOC(
    h5geo=h5geo, 
    container=con, 
    key="receiver", 
    xkey=geom_rec_xkey, 
    ykey=geom_rec_ykey, 
    zkey=geom_rec_zkey, 
    do_coord_transform=true, 
    model_origin_x=model_origin_x, 
    model_origin_y=model_origin_y, 
    model_orientation=model_orientation)
  srcGeometry = H5GeometryOOC(
    h5geo=h5geo, 
    container=con, 
    key="source", 
    xkey=geom_src_xkey, 
    ykey=geom_src_ykey, 
    zkey=geom_src_zkey, 
    do_coord_transform=false, 
    model_origin_x=model_origin_x, 
    model_origin_y=model_origin_y, 
    model_orientation=model_orientation)

  if isnothing(srcGeometry) || isnothing(srcGeometry)
    @error "Unable to prepare Source and Receiver JUDI Geometry objects from H5Seis\n"
    return
  end

  # in case computation type is RTM/FWI dobs is observed data
  dobs = judiVector(recGeometry, con)

  dt = h5geom.getSampRate("ms")
  nt = h5geom.getNSamp()
  t = dt*(nt-1)

  # print modeling info
  @info "computation_type: $computation_type\n"
  @info "model settings:\n"
  @info "model_xkey: $model_xkey\n"
  @info "model_ykey: $model_ykey\n"
  @info "model_origin_x: $model_origin_x\n"
  @info "model_origin_y: $model_origin_y\n"
  @info "model_orientation: $model_orientation\n"
  @info "geometry settings:\n"
  @info "geom_src_xkey: $geom_src_xkey\n"
  @info "geom_src_ykey: $geom_src_ykey\n"
  @info "geom_src_zkey: $geom_src_zkey\n"
  @info "t0: 0\n"
  @info "t1: $t\n"
  @info "dt: $dt\n"
  @info "src_frq (kHz): $src_frq"
  @info "nsrc: $(length(con))\n"
  @info "pkey: $(con.pkey)\n"
  @info "pkeyvals: $(con.pkeyvals)\n"

  # setup wavelet
  wavelet = ricker_wavelet(t, dt, src_frq)
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
    @error "TWRI not ready yet\n"
  end

end