# vmin/vmax in km/s
function H5FWI(;
  model::JUDI.Model,
  opt::JUDI.Options,
  q::JUDI.judiVector,
  dobs::JUDI.judiVector,
  niterations::Int,
  batchsize::Int,
  vmin::Number,
  vmax::Number)

  h5geo = pyimport("h5geopy._h5geo")
  
  # Bound constraints: Slowness squared [s^2/km^2]
  mmin = vec(ones(Float32, model.n) .* (1f0 / vmax)^2)
  mmax = vec(ones(Float32, model.n) .* (1f0 / vmin)^2)

  ############################### FWI ###########################################
  F = judiModeling(deepcopy(model), q.geometry, dobs.geometry; options=opt)

  # Optimization parameters
  if niterations < 1
    @error "Number of iterations can't be less that 1: $niterations\n"
    return
  end

  if batchsize < 1
    @error "Batch size can't be less that 1: $batchsize\n"
    return
  end

  # Prevent to access out of boundary elements
  if batchsize > dobs.nsrc
    batchsize = dobs.nsrc
  end

  # NLopt objective function
  println("No.  ", "fval         ", "norm(gradient)")
  function f!(x,grad)

      # Update model
      model.m .= convert(Array{Float32, 2}, reshape(x, model.n))

      # Seclect batch and calculate gradient
      i = randperm(dobs.nsrc)[1:batchsize]
      fval, gradient = fwi_objective(model, q[i], dobs[i])

      # Reset gradient in water column to zero
      gradient = reshape(gradient, model.n)
      gradient[:, 1:21] .= 0f0
      grad[1:end] = vec(gradient)

      global count; count += 1
      println(count, "    ", fval, "    ", norm(grad))
      return convert(Float64, fval)
  end

  # Optimization parameters
  nlopt = Opt(:LD_LBFGS, prod(model.n))
  lower_bounds!(nlopt, mmin); upper_bounds!(nlopt, mmax)
  min_objective!(nlopt, f!)
  maxeval!(nlopt, parse(Int, get(ENV, "NITER", "10")))
  (minf, minx, ret) = optimize(nlopt, vec(model.m.data))

  H5WritePhysicalParameter(
    cntName="$(opt.file_path)/fwi.h5",
    objName="fwi",
    cntCreationType=h5geo.CreationType.OPEN_OR_CREATE,
    objCreationType=h5geo.CreationType.CREATE_UNDER_NEW_NAME,
    php=sqrt.(1f0./model.m)
  )
end