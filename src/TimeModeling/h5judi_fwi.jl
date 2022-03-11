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
  
  # Bound constraints: Slowness squared [s^2/km^2]
  mmin = vec(ones(Float32, model.n) .* (1f0 / vmax)^2)
  mmax = vec(ones(Float32, model.n) .* (1f0 / vmin)^2)

  ############################### FWI ###########################################
  F = judiModeling(deepcopy(model), q.geometry, dobs.geometry; options=opt)

  # Optimization parameters
  if niterations < 1
    @error "Number of iterations can't be less that 1: $niterations"
    return
  end

  if batchsize < 1
    @error "Batch size can't be less that 1: $batchsize"
    return
  end

  fhistory_SGD = zeros(Float32, niterations)

  # Projection operator for bound constraints
  proj(x) = reshape(median([vec(mmin) vec(x) vec(mmax)]; dims=length(model.n)), model.n)
  ls = BackTracking(order=3, iterations=10, )

  # Main loop
  for j=1:niterations
    # get fwi objective function value and gradient
    i = randperm(dobs.nsrc)[1:batchsize]
    fval, gradient = fwi_objective(model, q[i], dobs[i])
    p = -gradient/norm(gradient, Inf)
    
    println("FWI iteration no: ",j,"; function value: ",fval)
    fhistory_SGD[j] = fval

    # linesearch
    function ϕ(α)
      F.model.m .= proj(model.m .+ α * p)
      misfit = .5*norm(F[i]*q[i] - dobs[i])^2
      @show α, misfit
      return misfit
    end
    step, fval = ls(ϕ, 1f-1, fval, dot(gradient, p))

    # Update model and bound projection
    model.m .= proj(model.m .+ step .* p)
  end

  figure(); imshow(sqrt.(1f0./adjoint(model.m))); title("FWI with SGD")
end