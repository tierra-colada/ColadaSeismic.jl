

function prepare_JUDI_geometry_2D(
  h5geom::PyCall.PyObject;
  src_xkey_min::Float64=-Inf,
  src_xkey_max::Float64=Inf,
  src_xkey::String="SRCX",
  src_zkey::String="SES",
  rec_xkey::String="GRPX",
  rec_zkey::String="RGE")

  # pylogging = pyimport("logging")
  # h5gt = pyimport("h5gtpy._h5gt")
  # h5geo = pyimport("h5geopy._h5geo")
  
  usrc_x = h5geom.getPKeyValues(src_xkey, h5geom.getLengthUnits(), "m")
  if isempty(usrc_x)
    @error "Unable to get src_xkey PKey values. Probably $src_xkey is sorting"
    return
  end

  usrc_x = usrc_x[(usrc_x .>= src_xkey_min) .& (usrc_x .<= src_xkey_max)] # brackets are necessary
  if isempty(usrc_x)
    @error "Search $src_xkey >= $src_xkey_min and $src_xkey >= $src_xkey_max resulted in an empty array"
    return
  end
  
  # receiver sampling and recording time
  timeR = abs(h5geom.getLastSample(0, "ms"))   # receiver recording time [ms]
  dtR = abs(h5geom.getSampRate("ms"))          # receiver sampling interval [ms]
  
  # source sampling and number of time steps
  timeS = timeR  # ms
  dtS = dtR   # ms
  
  nsrc = length(usrc_x)
  keylist = [src_xkey]
  srcXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  # srcYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srcZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  # recYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  indCellVec = Vector{Vector{Int64}}(undef, nsrc)
  for i in 1:nsrc
    minlist = [usrc_x[i]]
    maxlist = [usrc_x[i]]
    ~, ~, indCellVec[i] = h5geom.getSortedData(keylist, minlist, maxlist, 0, 0)
  
    # Set up source geometry (cell array with source locations for each shot)
    srcXCellVec[i] = [usrc_x[i]]
    # srcYCellVec[i] = [0]
    srcZCellVec[i] = -h5geom.getTraceHeader(src_zkey, indCellVec[i][1], 1, h5geom.getLengthUnits(), "m")[:]

    # Set up receiver structure
    recXCellVec[i] = h5geom.getTraceHeader([rec_xkey], indCellVec[i], [h5geom.getLengthUnits()], ["m"])[:]  # [:] converts matrix to vector
    # recYCellVec[i] = zeros(length(recXCellVec[i]))
    recZCellVec[i] = -h5geom.getTraceHeader([rec_zkey], indCellVec[i], [h5geom.getLengthUnits()], ["m"])[:]
  end

  srcGeometry = Geometry(srcXCellVec, 0.0, srcZCellVec; dt=dtS, t=timeS)
  recGeometry = Geometry(recXCellVec, 0.0, recZCellVec; dt=dtR, t=timeR)
  return srcGeometry, recGeometry, indCellVec
end



function prepare_JUDI_geometry_3D(
  h5geom::PyCall.PyObject;
  src_xkey_min::Float64=-Inf,
  src_xkey_max::Float64=Inf,
  src_ykey_min::Float64=-Inf,
  src_ykey_max::Float64=Inf,
  src_xkey::String="SRCX",
  src_ykey::String="SRCY",
  src_zkey::String="SES",
  rec_xkey::String="GRPX",
  rec_ykey::String="GRPY",
  rec_zkey::String="RGE",
  model_orientation::Number = 0.0)

  # pylogging = pyimport("logging")
  # h5gt = pyimport("h5gtpy._h5gt")
  h5geo = pyimport("h5geopy._h5geo")
  

  keylist = [src_xkey, src_ykey]
  minlist = [src_xkey_min, src_ykey_min]
  maxlist = [src_xkey_max, src_ykey_max]
  ~, src_xy, ~ = h5geom.getSortedData(keylist, minlist, maxlist, 0, 0, true, "", "m", true)
  if isempty(src_xy)
    @error "Unable to get sorted src_x src_y headers. Probably $src_xkey is missing"
    return
  end

  ~, usrc_xy, ~ = h5geo.sort_rows_unique(src_xy)

  # receiver sampling and recording time
  timeR = abs(h5geom.getLastSample(0, "ms"))   # receiver recording time [ms]
  dtR = abs(h5geom.getSampRate("ms"))          # receiver sampling interval [ms]
  
  # source sampling and number of time steps
  timeS = timeR  # ms
  dtS = dtR   # ms
  
  nsrc = size(usrc_xy)[1]
  keylist = [src_xkey, src_ykey]
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
    ~, ~, indCellVec[i] = h5geom.getSortedData(keylist, minlist, maxlist, 0, 0)
  
    # Set up source geometry (cell array with source locations for each shot)
    src_xy_rot = rotate_2xN_array(usrc_xy[i,:], -model_orientation)
    srcXCellVec[i] = [src_xy_rot[1]]
    srcYCellVec[i] = [src_xy_rot[2]]
    srcZCellVec[i] = -h5geom.getTraceHeader(src_zkey, indCellVec[i][1], 1, h5geom.getLengthUnits(), "m")[:]

    # Set up receiver structure
    rec_xy = h5geom.getTraceHeader([rec_xkey, rec_ykey], indCellVec[i], [h5geom.getLengthUnits(), h5geom.getLengthUnits()], ["m", "m"])
    rec_xy_rot = rotate_Nx2_array(rec_xy, -model_orientation)
    recXCellVec[i] = rec_xy_rot[:,1]
    recYCellVec[i] = rec_xy_rot[:,2]
    recZCellVec[i] = -h5geom.getTraceHeader([rec_zkey], indCellVec[i], [h5geom.getLengthUnits()], ["m"])[:]
  end

  srcGeometry = Geometry(srcXCellVec, srcYCellVec, srcZCellVec; dt=dtS, t=timeS)
  recGeometry = Geometry(recXCellVec, recYCellVec, recZCellVec; dt=dtR, t=timeR)
  return srcGeometry, recGeometry, indCellVec
end