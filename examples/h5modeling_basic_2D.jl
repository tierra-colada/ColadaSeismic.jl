# This example demonstrates how to use custom h5geo data containers within JUDI.
# The idea of the example is to generate SEGY files using JUDI API
# and after that use generated SEGY files to prepare HDF5 containers with h5geo library.
# After the h5geo containers are created they are planned to be used to calculate FWI, RSTM, TWRI.
# The example creates 'ColadaSeismic/examles/tmp' folder and stores all the calculated data there.
# For now I tested that with JUDI 'time-tests' branch wich is not finished yet and 
# it throws exception when JUDI.Options::save_data_to_disk is 'true'.
# Thus I distribute previously modeled SEGY files with 'master' branch


using ColadaSeismic, JUDI, SegyIO, LinearAlgebra, PyCall, PyPlot

# 'h5map_segy' maps SEGY to HDF5 h5geo::H5Seis object so that the same API is provided
# but no any data is read from SEGY
function h5map_segy(segy_files::Array{String})
  if length(segy_files) < 1
    @error "No SEGY files were provided"
    return
  end

  h5geo = pyimport("h5geopy._h5geo")
  if isnothing(h5geo)
    @error "Unable to import h5geo"
    return
  end

  container_name = splitext(segy_files[1])[1]
  container_name *= "_mapped.h5"
  container = h5geo.createSeisContainerByName(container_name, h5geo.CreationType.CREATE_UNDER_NEW_NAME)
  if isnothing(container)
    @error "Unable to create HDF5 container: $container_name"
    return
  end

  # parameters to create seis object
  p = h5geo.SeisParam()
  # Cant understand why in julia inheritance from H5BaseObject doesn't work
  # p.lengthUnits = "cm"
  # p.temporalUnits = "microsecond"
  p.domain = h5geo.Domain.TWT
  p.dataType = h5geo.SeisDataType.PRESTACK
  p.surveyType = h5geo.SurveyType.TWO_D
  p.trcChunk = 5
  p.segyFiles = segy_files
  p.mapSEGY = true

  @info "segy_files[1]: $(segy_files[1])"
  @info "segy_files[2]: $(segy_files[2])"

  seis_name = "mapped_segy"
  seis = container.createSeis(seis_name, p, h5geo.CreationType.CREATE_OR_OVERWRITE)
  if isnothing(seis)
    @error "Unable to create seis"
    return
  end

  seis.setLengthUnits("cm")
  seis.setTemporalUnits("microsecond")
  return seis
end

# 'h5read_segy' read SEGY files to HDF5 h5geo::H5Seis object
function h5read_segy(segy_files::Array{String})
  if length(segy_files) < 1
    @error "No SEGY files were provided"
    return
  end

  h5geo = pyimport("h5geopy._h5geo")
  if isnothing(h5geo)
    @error "Unable to import h5geo"
    return
  end

  container_name = splitext(segy_files[1])[1]
  container_name *= "_read.h5"
  container = h5geo.createSeisContainerByName(container_name, h5geo.CreationType.CREATE_UNDER_NEW_NAME)
  if isnothing(container)
    @error "Unable to create HDF5 container: $container_name"
    return
  end

  # parameters to create seis object
  p = h5geo.SeisParam()
  # Cant understand why in julia inheritance from H5BaseObject doesn't work
  # p.lengthUnits = "cm"
  # p.temporalUnits = "microsecond"
  p.domain = h5geo.Domain.TWT
  p.dataType = h5geo.SeisDataType.PRESTACK
  p.surveyType = h5geo.SurveyType.TWO_D
  # explicitely set 'nSamp' so that chunking along 1-st dim is equal to 'nSamp'
  p.nSamp = h5geo.getSEGYNSamp(segy_files[1], h5geo.getSEGYEndian(segy_files[1]))
  # 'nTrc' must be bigger that 'trcChunk'
  p.nTrc = 10
  p.trcChunk = 5

  @info "segy_files[1]: $(segy_files[1])"
  @info "segy_files[2]: $(segy_files[2])"
  @info "p.nSamp: $(p.nSamp)"

  seis_name = "read_segy"
  seis = container.createSeis(seis_name, p, h5geo.CreationType.CREATE_OR_OVERWRITE)
  if isnothing(seis)
    @error "Unable to create seis"
    return
  end

  seis.readSEGYTraces(segy_files)
  seis.setLengthUnits("cm")
  seis.setTemporalUnits("microsecond")
  return seis
end

# Set up model structure
n = (120, 100)   # (x,y,z) or (x,z)
d = (10., 10.)
o = (0., 0.)

# Velocity [km/s]
v = ones(Float32,n) .+ 0.5f0
v0 = ones(Float32,n) .+ 0.5f0
v[:,Int(round(end/2)):end] .= 3.5f0
rho = (v0 .+ .5f0) ./ 2

# Slowness squared [s^2/km^2]
m = (1f0 ./ v).^2
m0 = (1f0 ./ v0).^2
dm = vec(m - m0)

# Setup info and model structure
nsrc = 2	# number of sources
model = Model(n, d, o, m)
model0 = Model(n, d, o, m0)

# Set up receiver geometry
nxrec = 120
xrec = range(50f0, stop=1150f0, length=nxrec)
yrec = 0f0
zrec = range(50f0, stop=50f0, length=nxrec)

# receiver sampling and recording time
timeR = 1000f0   # receiver recording time [ms]
dtR = 2f0    # receiver sampling interval [ms]

# Set up receiver structure
recGeometry = Geometry(xrec, yrec, zrec; dt=dtR, t=timeR, nsrc=nsrc)

# Set up source geometry (cell array with source locations for each shot)
xsrc = convertToCell(range(400f0, stop=800f0, length=nsrc))
ysrc = convertToCell(range(0f0, stop=0f0, length=nsrc))
zsrc = convertToCell(range(200f0, stop=200f0, length=nsrc))

# source sampling and number of time steps
timeS = 1000f0  # ms
dtS = 2f0   # ms

# Set up source structure
srcGeometry = Geometry(xsrc, ysrc, zsrc; dt=dtS, t=timeS)

# setup wavelet
f0 = 0.01f0     # kHz
wavelet = ricker_wavelet(timeS, dtS, f0)
q = judiVector(srcGeometry, wavelet)

# Set up info structure for linear operators
ntComp = get_computational_nt(srcGeometry, recGeometry, model)
info = Info(prod(n), nsrc, ntComp)

###################################################################################################

segy_path = @__DIR__
segy_path *= "/tmp/"
segy_name = "shot"

# rm(segy_path, force=true, recursive=true)
# mkpath(segy_path)

# Write shots as segy files to disk (gives error in 'time-tests' branch)
# opt = Options(optimal_checkpointing=false, isic=false, subsampling_factor=2, dt_comp=1.0,
#               save_data_to_disk=true, file_path=segy_path, file_name=segy_name)
opt = Options(optimal_checkpointing=false, isic=false, subsampling_factor=2, dt_comp=1.0)

# Setup operators
Pr = judiProjection(info, recGeometry)
F = judiModeling(info, model; options=opt)
F0 = judiModeling(info, model0; options=opt)
Ps = judiProjection(info, srcGeometry)
J = judiJacobian(Pr*F0*adjoint(Ps), q)

# Nonlinear modeling
dobs = Pr*F*adjoint(Ps)*q

# modeled SEGY files
files = readdir(segy_path, join=true)
segy_files = Array{String}(undef, 0)
for file in files
  if splitext(file)[2] == ".segy"
    push!(segy_files, file)
  end
end

# Ether map SEGY (only ieee32 supported) or read (any 4-byte format) it

# seis = h5map_segy(segy_files)
seis = h5read_segy(segy_files)
if isnothing(seis)
  @error "Unable to create seis"
  return
end

# To retrieve sorted data we need to prepare sorting from primary keys (PKey)
# 'h5geo' uses trace header names that differs from the ones that uses SegyIO
if !seis.hasPKeySort("SRCX")
  seis.addPKeySort("SRCX")
end

# create H5SeisCon and judiVector from it
con = H5SeisCon(seis=seis, pkey="SRCX")
h5dobs = judiVector(con)

# as long as SEGY stores coordinates in Int32 and h5geo stores it in Float64
# we need to set previously calculated 'recGeometry' to prevent round off mismatch error (needed only for this example)
h5dobs.geometry = recGeometry

##################################################################################################

# Continue colaculations  (FWI, RSTM etc) with 'judiVector{H5SeisCon}'

# Adjoint with h5dobs
pdeFull = Ps*adjoint(F)*adjoint(Pr)
qad = pdeFull*h5dobs

# Linearized modeling
dD = J*dm
# Adjoint jacobian
rtm = adjoint(J)*dD

# evaluate FWI objective function with h5dobs
f, g = fwi_objective(model0, q, h5dobs; options=opt)

imshow(g'); title("FWI (g')")
savefig(segy_path * "/g.png")

# TWRI with h5dobs
f, gm, gy = twri_objective(model0, q, h5dobs, nothing; options=opt, optionswri=TWRIOptions(params=:all))
f, gm = twri_objective(model0, q, h5dobs, nothing; options=Options(frequencies=[[.009, .011], [.008, .012]]),
                                                 optionswri=TWRIOptions(params=:m))

imshow(gm'); title("TWRI (gm')")
savefig(segy_path * "/gm.png")

# evaluate LSRTM objective function with h5dobs
fj, gj = lsrtm_objective(model0, q, dD, dm; options=opt)
fjn, gjn = lsrtm_objective(model0, q, h5dobs, dm; nlind=true, options=opt)

imshow(gj'); title("LSRTM (gj')")
savefig(segy_path * "/gj.png")
imshow(gjn'); title("LSRTM (gjn')")
savefig(segy_path * "/gjn.png")

# By extension, lsrtm_objective is the same as fwi_objecive when `dm` is zero
# And with computing of the residual. Small noise can be seen in the difference
# due to floating point roundoff errors with openMP, but running with 
# OMP_NUM_THREAS=1 (no parllelism) produces the exact (difference == 0) same result
# gjn2 == g
fjn2, gjn2 = lsrtm_objective(model0, q, h5dobs, 0f0.*dm; nlind=true, options=opt)

imshow(gjn2'); title("LSRTM (gjn2')")
savefig(segy_path * "/gjn2.png")
