export H5Modeling

# operation_type: Forward Modeling, RTM, FWI, TWRI
function H5Modeling(;
  # models
  h5vel,
  h5den=nothing,
  h5qf=nothing,
  model_nb::Real=40,  # number of absorbing boundaries points on each side
  h5eps=nothing,
  h5delta=nothing,
  h5tetha=nothing,
  h5phi=nothing,
  model_xkey::String="CDP_X",
  model_ykey::String="CDP_Y",

  # geometry
  geom_con=nothing,
  geom_src_pkey::String="SP",
  geom_src_xkey::String="SRCX",
  geom_src_ykey::String="SRCY",
  geom_src_zkey::String="SES",
  geom_rec_xkey::String="GRPX",
  geom_rec_ykey::String="GRPY",
  geom_rec_zkey::String="RGE",

  # judi
  opt::JUDI.JUDIOptions,

  # common settings
  computation_type::String,
  nt::Real,
  dt::Number,
  save_as::SaveAs,
  do_coord_transform::Bool,
  spatial_reference::String,
  
  # LSRTM
  lsrtm_niter::Int=nothing,
  lsrtm_batchsize::Int=nothing,
  # FWI
  fwi_niter::Int=nothing,
  fwi_batchsize::Int=nothing,
  fwi_vmin::Number,
  fwi_vmax::Number,
  # TWRI
  twri_opt::TWRIOptions=nothing,
  )

  global h5geo = pyimport("h5geopy._h5geo")
  h5geo.sr.setSpatialReferenceFromUserInput(spatial_reference)
  h5geo.sr.setLengthUnits("m")
  h5geo.sr.setTemporalUnits("ms")

  global spatial_reference = spatial_reference
  global save_as = save_as

  @info "H5Modeling started\n"

  if isnothing(h5vel)
    @error "No velocity model\n"
    return
  end

  if isnothing(opt)
    @error "No JUDI Options for modeling\n"
    return
  end

  if !isnothing(opt.dt_comp) && opt.dt_comp <= 0
    @error "dt_comp can't be less or equal to zero\n"
    return
  end

  # prepare physical parameters
  phpvel, global model_orientation = H5ReadPhysicalParameter(
    h5vel, phptype=VELOCITY, xkey=model_xkey, ykey=model_ykey)
  if isnothing(phpvel)
    @error "Unable to read PhysicalParameter: VELOCITY\n"
    return
  end

  # declare PHP vars to be able use thm as arguments (or unknown var error will appear)
  phpdensity = nothing
  phpqualityFactor = nothing
  phpepsilon = nothing
  phpdelta = nothing
  phptheta = nothing
  phpphi = nothing

  if !isnothing(h5den)
    phpdensity, _ = H5ReadPhysicalParameter(
      h5den, phptype=DENSITY, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdensity)
      @error "Unable to read PhysicalParameter: DENSITY\n"
      return
    end
  end

  if !isnothing(h5qf)
    phpqualityFactor, _ = H5ReadPhysicalParameter(
      h5qf, phptype=QUALITYFACTOR, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpqualityFactor)
      @error "Unable to read PhysicalParameter: QUALITYFACTOR\n"
      return
    end
  end

  if !isnothing(h5eps)
    phpepsilon, _ = H5ReadPhysicalParameter(
      h5eps, phptype=EPSILON, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpepsilon)
      @error "Unable to read PhysicalParameter: EPSILON\n"
      return
    end
  end

  if !isnothing(h5delta)
    phpdelta, _ = H5ReadPhysicalParameter(
      h5delta, phptype=DELTA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpdelta)
      @error "Unable to read PhysicalParameter: DELTA\n"
      return
    end
  end

  if !isnothing(h5tetha)
    phptheta, _ = H5ReadPhysicalParameter(
      h5tetha, phptype=THETA, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phptheta)
      @error "Unable to read PhysicalParameter: THETA\n"
      return
    end
  end

  if !isnothing(h5phi)
    phpphi, _ = H5ReadPhysicalParameter(
      h5phi, phptype=PHI, xkey=model_xkey, ykey=model_ykey)
    if isnothing(phpphi)
      @error "Unable to read PhysicalParameter: PHI\n"
      return
    end
  end

  model = JUDI.Model(phpvel.n, phpvel.d, phpvel.o, phpvel.data, 
          rho=isnothing(phpdensity) ? nothing : phpdensity.data, 
          qp=isnothing(phpqualityFactor) ? nothing : phpqualityFactor.data, 
          epsilon=isnothing(phpepsilon) ? nothing : phpepsilon.data, 
          delta=isnothing(phpdelta) ? nothing : phpdelta.data, 
          theta=isnothing(phptheta) ? nothing : phptheta.data, 
          phi=isnothing(phpphi) ? nothing : phpphi.data,
          nb=model_nb)

  # Geometry
  if isnothing(geom_con)
    @error "Geometry container is NULL\n"
    return
  end

  global model_origin_x = model.o[1]
  global model_origin_y = length(model.o == 2) ? 0.0 : model.o[2]
  global survey_type = length(model.n) == 2 ? h5geo.SurveyType.TWO_D : h5geo.SurveyType.THREE_D

  if geom_con isa SeisCon
    segy_depth_key_rec = h5geo2judiTraceHeaderName(geom_rec_zkey)
    if isnothing(segy_depth_key_rec)
      @error "Invalid geometry receiver ZKEY"
      return
    end
    segy_depth_key_src = h5geo2judiTraceHeaderName(geom_src_zkey)
    if isnothing(segy_depth_key_src)
      @error "Invalid geometry source ZKEY"
      return
    end
    recGeometry = Geometry(geom_con; key="receiver", segy_depth_key=segy_depth_key_rec)
    srcGeometry = Geometry(geom_con; key="source", segy_depth_key=segy_depth_key_src)
  elseif geom_con isa H5SeisCon
    recGeometry = H5GeometryOOC(
      container=geom_con, 
      key="receiver", 
      xkey=geom_rec_xkey, 
      ykey=geom_rec_ykey, 
      zkey=geom_rec_zkey, 
      do_coord_transform=do_coord_transform)
    srcGeometry = H5GeometryOOC(
      container=geom_con, 
      key="source", 
      xkey=geom_src_xkey, 
      ykey=geom_src_ykey, 
      zkey=geom_src_zkey, 
      do_coord_transform=do_coord_transform)
  else
    @error "Geometry container is neither SeisCon nor H5SeisCon\n"
    return
  end

  # prepare seis_out for storing traces
  if save_as == SaveAs::H5SEIS
    global seis_out_cnt = h5geo.createSeisContainerByName(
      joinpath(opt.file_path, opt.file_name), 
      h5geo.CreationType.OPEN_OR_CREATE)
    if isnothing(seis_out_cnt)
      @error "Unable to open or create seis container for storing output data\n"
      return
    end
  end

  if isnothing(srcGeometry) || isnothing(srcGeometry)
    @error "Unable to prepare Source and Receiver JUDI Geometry objects from H5Seis\n"
    return
  end

  # setup wavelet
  t = Float32((nt-1)*dt)
  wavelet = ricker_wavelet(t, dt, opt.f0)
  q = judiVector(srcGeometry, wavelet)

  # JUDI (Devito) doesn't create a dir to save data
  if opt.save_data_to_disk
    mkpath(opt.file_path)
  end

  # simply to print the value
  ntComp = get_computational_nt(srcGeometry, recGeometry, model, dt=opt.dt_comp)
  dtComp = calculate_dt(model)

  # print modeling info
  @info "computation_type: $computation_type\n"
  @info "spatial_reference: $spatial_reference\n"
  @info "save_as: $save_as\n"
  @info "nt: $nt"
  @info "ntComp: $ntComp\n"
  @info "dt: $dt"
  @info "dtComp: $dtComp\n"
  @info "model settings:\n"
  @info "model.o: $(model.o)\n"
  @info "model.n: $(model.n)\n"
  @info "model.d: $(model.d)\n"
  @info "model_xkey: $model_xkey\n"
  @info "model_ykey: $model_ykey\n"
  @info "model_origin_x: $model_origin_x\n"
  @info "model_origin_y: $model_origin_y\n"
  @info "survey_type: $survey_type\n"
  @info "model_orientation: $model_orientation\n"
  @info "geometry settings:\n"
  @info "geom_src_xkey: $geom_src_xkey\n"
  @info "geom_src_ykey: $geom_src_ykey\n"
  @info "geom_src_zkey: $geom_src_zkey\n"
  @info "geom_rec_xkey: $geom_rec_xkey\n"
  @info "geom_rec_ykey: $geom_rec_ykey\n"
  @info "geom_rec_zkey: $geom_rec_zkey\n"
  @info "dt: $dt\n"
  @info "t: $t\n"
  @info "src_frq (kHz): $(opt.f0)"
  @info "nsrc: $(length(con))\n"
  @info "pkey: $(con.pkey)\n"
  @info "pkeyvals: $(con.pkeyvals)\n"
  @info "options: $opt\n"

  if computation_type == "Forward Modeling"
    H5ForwardModeling(model=model, opt=opt, q=q, recGeometry=recGeometry)
  elseif computation_type == "RTM"
    dobs = judiVector(recGeometry, con)
    mkpath(opt.file_path)
    H5RTM(model=model, opt=opt, q=q, dobs=dobs)
  elseif computation_type == "LSRTM"
    @error "LSRTM not ready yet\n"
  elseif computation_type == "FWI"
    dobs = judiVector(recGeometry, con)
    mkpath(opt.file_path)
    H5FWI(model=model, opt=opt, q=q, dobs=dobs, 
          niterations=fwi_niter, batchsize=fwi_batchsize, 
          vmin=fwi_vmin, vmax=fwi_vmax)
  elseif computation_type == "TWRI"
    @error "TWRI not ready yet\n"
  end

end