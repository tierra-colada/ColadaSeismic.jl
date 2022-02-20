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
    @error "No $src_xkey sorting"
    return
  end

  usrc_x = usrc_x[(usrc_x .>= src_xkey_min) .& (usrc_x .<= src_xkey_max)] # brackets are necessary
  if isempty(usrc_x)
    @error "Search $src_xkey >= $src_xkey_min and $src_xkey >= $src_xkey_max resulted in an empty array"
    return
  end
  
  # receiver sampling and recording time
  timeR = h5geom.getLastSample(0, "ms") * (-1)   # receiver recording time [ms]
  dtR = h5geom.getSampRate("ms") * (-1)    # receiver sampling interval [ms]
  
  # source sampling and number of time steps
  timeS = timeR  # ms
  dtS = dtR   # ms
  
  nsrc = length(usrc_x)
  keylist = [src_xkey, rec_xkey]
  srcXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  # srcYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  srczCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recXCellVec = Vector{Vector{Float32}}(undef, nsrc)
  # recYCellVec = Vector{Vector{Float32}}(undef, nsrc)
  recZCellVec = Vector{Vector{Float32}}(undef, nsrc)
  indCellVec = Vector{Vector{Int64}}(undef, nsrc)
  for i in 1:nsrc
    minlist = [usrc_x[i], -Inf]
    maxlist = [usrc_x[i], Inf]
    ~, ~, indCellVec[i] = h5geom.getSortedData(keylist, minlist, maxlist, 0, 0)
  
    # Set up source geometry (cell array with source locations for each shot)
    srcXCellVec[i] = [usrc_x[i]]
    # srcYCellVec[i] = [0]
    srczCellVec[i] = -h5geom.getTraceHeader(src_zkey, indCellVec[i][1], 1, h5geom.getLengthUnits(), "m")[:]

    # Set up receiver structure
    recXCellVec[i] = h5geom.getTraceHeader([rec_xkey], indCellVec[i], [h5geom.getLengthUnits()], ["m"])[:]  # [:] converts matrix to vector
    # recYCellVec[i] = zeros(length(recXCellVec[i]))
    recZCellVec[i] = -h5geom.getTraceHeader([rec_zkey], indCellVec[i], [h5geom.getLengthUnits()], ["m"])[:]
  end

  srcGeometry = Geometry(srcXCellVec, 0.0, srczCellVec; dt=dtS, t=timeS)
  recGeometry = Geometry(recXCellVec, 0.0, recZCellVec; dt=dtR, t=timeR)
  return srcGeometry, recGeometry, indCellVec
end