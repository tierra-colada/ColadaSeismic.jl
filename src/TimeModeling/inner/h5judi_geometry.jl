function H5ReadGeometry2D(
  h5obj::PyCall.PyObject;
  h5opt::H5GeomOptions=H5GeomOptions(),
  src_dt::Number=2,
  src_nt::Integer=50,
  rec_dt::Number=2,
  rec_nt::Integer=50)

  if src_dt <= 0 || src_nt < 1 || rec_dt <= 0 || rec_nt < 1
    @error "Source/Receiver timings negative"
    return
  end
  
  usrc_x = h5obj.getPKeyValues(h5opt.src_xkey, h5obj.getLengthUnits(), "m")
  if isempty(usrc_x)
    @error "Unable to get src_xkey PKey values. Probably $(h5opt.src_xkey) is sorting"
    return
  end

  usrc_x = usrc_x[(usrc_x .>= h5opt.src_xkey_min) .& (usrc_x .<= h5opt.src_xkey_max)] # brackets are necessary
  if isempty(usrc_x)
    @error "Search $(h5opt.src_xkey) >= $(h5opt.src_xkey_min) and $(h5opt.src_xkey) >= $(h5opt.src_xkey_max) resulted in an empty array"
    return
  end
  
  # receiver sampling and recording time
  src_t1 = (src_nt-1)*src_dt
  rec_t1 = (rec_nt-1)*rec_dt

  nsrc = length(usrc_x)
  keylist = [h5opt.src_xkey]
  srcXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srcYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srcZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  indCellVec = Vector{Vector{Int64}}(undef, nsrc)
  for i in 1:nsrc
    minlist = [usrc_x[i]]
    maxlist = [usrc_x[i]]
    _, _, indCellVec[i] = h5obj.getSortedData(keylist, minlist, maxlist, 0, 0)
  
    # Set up source geometry (cell array with source locations for each shot)
    srcXCellVec[i] = [usrc_x[i]]
    srcYCellVec[i] = [0]
    srcZCellVec[i] = -h5obj.getTraceHeader(h5opt.src_zkey, indCellVec[i][1], 1, h5obj.getLengthUnits(), "m")[:]

    # Set up receiver structure
    recXCellVec[i] = h5obj.getTraceHeader([h5opt.rec_xkey], indCellVec[i], [h5obj.getLengthUnits()], ["m"])[:]  # [:] converts matrix to vector
    recYCellVec[i] = zeros(length(recXCellVec[i]))
    recZCellVec[i] = -h5obj.getTraceHeader([h5opt.rec_zkey], indCellVec[i], [h5obj.getLengthUnits()], ["m"])[:]
  end

  srcGeometry = Geometry(srcXCellVec, srcYCellVec, srcZCellVec; dt=src_dt, t=src_t1)
  recGeometry = Geometry(recXCellVec, recYCellVec, recZCellVec; dt=rec_dt, t=rec_t1)
  return srcGeometry, recGeometry, indCellVec
end



function H5ReadGeometry3D(
  h5obj::PyCall.PyObject;
  h5opt::H5GeomOptions=H5GeomOptions(),
  src_dt::Number=2,
  src_nt::Integer=50,
  rec_dt::Number=2,
  rec_nt::Integer=50)

  h5geo = pyimport("h5geopy._h5geo")

  if src_dt <= 0 || src_nt < 1 || rec_dt <= 0 || rec_nt < 1
    @error "Source/Receiver timings negative"
    return
  end

  keylist = [h5opt.src_xkey, h5opt.src_ykey]
  minlist = [h5opt.src_xkey_min, h5opt.src_ykey_min]
  maxlist = [h5opt.src_xkey_max, h5opt.src_ykey_max]
  _, src_xy, _ = h5obj.getSortedData(keylist, minlist, maxlist, 0, 0, true, "", "m", true)
  if isempty(src_xy)
    @error "Unable to get sorted src_x src_y headers. Probably $(h5opt.src_xkey) is missing"
    return
  end

  _, usrc_xy, _ = h5geo.sort_rows_unique(src_xy)

  # receiver sampling and recording time
  src_t1 = (src_nt-1)*src_dt
  rec_t1 = (rec_nt-1)*rec_dt
  
  nsrc = size(usrc_xy)[1]
  keylist = [h5opt.src_xkey, h5opt.src_ykey]
  srcXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srcYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srcZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  indCellVec = Vector{Vector{Int64}}(undef, nsrc)
  for i in 1:nsrc
    minlist = [usrc_xy[i,1], usrc_xy[i,2]]
    maxlist = [usrc_xy[i,1], usrc_xy[i,2]]
    _, _, indCellVec[i] = h5obj.getSortedData(keylist, minlist, maxlist, 0, 0)
  
    # Set up source geometry (cell array with source locations for each shot)
    src_xy_rot = rotate_2xN_array(usrc_xy[i,:], -h5opt.model_orientation)
    srcXCellVec[i] = [src_xy_rot[1]]
    srcYCellVec[i] = [src_xy_rot[2]]
    srcZCellVec[i] = -h5obj.getTraceHeader(h5opt.src_zkey, indCellVec[i][1], 1, h5obj.getLengthUnits(), "m")[:]

    # Set up receiver structure
    rec_xy = h5obj.getTraceHeader([h5opt.rec_xkey, h5opt.rec_ykey], indCellVec[i], [h5obj.getLengthUnits(), h5obj.getLengthUnits()], ["m", "m"])
    rec_xy_rot = rotate_Nx2_array(rec_xy, -h5opt.model_orientation)
    recXCellVec[i] = rec_xy_rot[:,1]
    recYCellVec[i] = rec_xy_rot[:,2]
    recZCellVec[i] = -h5obj.getTraceHeader([h5opt.rec_zkey], indCellVec[i], [h5obj.getLengthUnits()], ["m"])[:]
  end

  srcGeometry = Geometry(srcXCellVec, srcYCellVec, srcZCellVec; dt=src_dt, t=src_t1)
  recGeometry = Geometry(recXCellVec, recYCellVec, recZCellVec; dt=rec_dt, t=rec_t1)
  return srcGeometry, recGeometry, indCellVec
end