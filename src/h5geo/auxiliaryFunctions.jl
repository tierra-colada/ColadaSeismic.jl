import JUDI.get_computational_nt, JUDI.time_modeling

function get_computational_nt(srcGeometry::H5GeometryOOC, recGeometry::H5GeometryOOC, model::Model; dt=nothing)
  # Determine number of computational time steps
  nsrc = length(srcGeometry.container)
  nt = Array{Integer}(undef, nsrc)
  dtComp = calculate_dt(model; dt=dt)
  for j=1:nsrc
    dt_rec = recGeometry.dt[j]
    nt_rec = recGeometry.nt[j]
    dt_src = srcGeometry.dt[j]
    nt_src = srcGeometry.nt[j]
    ntRec = Int(dt_rec*(nt_rec-1) ÷ dtComp) + 1
    ntSrc = Int(dt_src*(nt_src-1) ÷ dtComp) + 1
    nt[j] = max(ntRec, ntSrc)
  end
  return nt
end

# Setup time-domain linear or nonlinear foward and adjoint modeling and interface to devito
function time_modeling(model_full::Model, srcGeometry::H5GeometryOOC, srcData::JUDI.ArrayOrNot,
                      recGeometry::H5GeometryOOC, recData::JUDI.ArrayOrNot, 
                      dm::JUDI.PhysOrNot, op::Symbol, options::JUDIOptions)
  # Load full geometry for out-of-core geometry containers
  recGeometryIC = Geometry(recGeometry)
  srcGeometryIC = Geometry(srcGeometry)

  # Reutrn directly for J*0
  if op==:born
    if norm(dm) == 0 && options.return_array == false
      return judiVector(recGeometryIC, zeros(Float32, recGeometryIC.nt[1], length(recGeometryIC.xloc[1])))
    elseif norm(dm) == 0 && options.return_array == true
      return vec(zeros(Float32, recGeometryIC.nt[1], length(recGeometryIC.xloc[1])))
    end
  end

  # limit model to area with sources/receivers
  if options.limit_m == true
    model = deepcopy(model_full)
    model, dm = limit_model_to_receiver_area(srcGeometryIC, recGeometryIC, model, options.buffer_size; pert=dm)
  else
    model = model_full
  end

  # Set up Python model structure
  modelPy = devito_model(model, options, dm)

  # Remove receivers outside the modeling domain (otherwise leads to segmentation faults)
  recGeometryIC, recData = remove_out_of_bounds_receivers(recGeometryIC, recData, model)

  # Devito interface
  argout = devito_interface(modelPy, srcGeometryIC, srcData, recGeometryIC, recData, dm, options)
  # Extend gradient back to original model size
  if op==:adjoint_born && options.limit_m==true
    argout = extend_gradient(model_full, model, argout)
  end

  if options.save_data_to_disk
    nSamp = size(argout.data[1])[1]
    nTrc = size(argout.data[1])[2]

    (length(srcGeometry.container[1].ind) != nTrc) && (@error "Number of indexes doesn't match the number of traces\n")
    if @time !srcGeometry.container.seis.writeTrace(argout.data[1],srcGeometry.container[1].ind,0)
      @warn "Unable to write traces for pkey: $(srcGeometry.container[1].pkey)\n
            pkeyval: $(srcGeometry.container[1].pkeyvals[1])\n"
    end
    srcGeometry.container.seis.getH5File().flush()
  end

  return argout
end

# Saving to disk utilities
# save_to_disk(shot, srcGeometry, srcData, options, ::Val) = shot
# save_to_disk(shot::judiVector, srcGeometry, srcData, options, ::Val{false}) = shot

# function save_to_disk(data::) 
#   shot.geometry[1]

#   container = write_shot_record(srcGeometry, srcData, shot.geometry[1], shot.data[1], options)
#   return judiVector(container)
# end

# function save_to_disk(shot::judiVector) 
#   shot.geometry[1]

#   container = write_shot_record(srcGeometry, srcData, shot.geometry[1], shot.data[1], options)
#   return judiVector(container)
# end

# function write_shot_record(srcGeometry::GeometryIC, srcData, recGeometry::GeometryIC, recData, options)
#   q = judiVector(srcGeometry, srcData)
#   d = judiVector(recGeometry, recData)
#   pos = [srcGeometry.xloc[1][1], srcGeometry.yloc[1][1],  srcGeometry.zloc[1][1]]
#   pos = join(["_"*string(trunc(p; digits=2)) for p in pos])
#   file = join([string(options.file_name), pos,".segy"])
#   block_out = judiVector_to_SeisBlock(d, q)
#   segy_write(join([options.file_path,"/",file]), block_out)
#   container = scan_file(join([options.file_path,"/",file]),
#                         ["GroupX", "GroupY", "dt", "SourceSurfaceElevation", "RecGroupElevation"];
#                         chunksize=256)
#   return container
# end
