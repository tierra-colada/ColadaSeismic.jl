import JUDI.split, JUDI.subsample

export H5SeisCon

mutable struct H5SeisCon
  seis::PyCall.PyObject
  pkey::String
  # 'pkeyvals' must be Float64 as roundoff may cause problems 
  # when retrieving data from 'seis' (getPKeyTraceSize for example)
  pkeyvals::Array{Float64,1}
end

function H5SeisCon(;
  seis::PyCall.PyObject=nothing,
  pkey::String="SP",
  pkeyvals=nothing)
  if isnothing(seis)
    @error "Seis is Nothing\n"
    return
  end

  if isnothing(pkey)
    @error "pkey is Nothing\n"
    return
  end

  if !seis.hasPKeySort(pkey)
    @info "Updating trace header limits and adding pkey sort: $pkey..."
    seis.updateTraceHeaderLimits(Int(1e7))
    seis.addPKeySort(pkey)
  end

  # must go after 'addPKeySort'
  if isnothing(pkeyvals)
    pkeyvals = seis.getPKeyValues(pkey)
  end

  return H5SeisCon(seis, pkey, pkeyvals)
end

size(con::H5SeisCon) = size(con.pkeyvals)
length(con::H5SeisCon) = length(con.pkeyvals)

function getindex(con::H5SeisCon, a::TA) where {TA<:Union{Array{<:Integer,1}, AbstractRange, Integer}} 
  read_H5SeisCon(con, a)
end

function getindex(con::H5SeisCon, a::Colon) 
  read_H5SeisCon(con, 1:length(con))
end

function split(s::H5SeisCon, inds::Union{Vector{Ti}, AbstractRange{Ti}}) where {Ti<:Integer}
  return H5SeisCon(s.seis, s.pkey, view(s.pkeyvals, inds))
end

split(s::H5SeisCon, inds::Integer) = split(s, [inds])

subsample(a::Array{H5SeisCon,1}, i::Int64) = a[i]

iterate(S::H5SeisCon, state::Integer=1) = state > length(S) ? nothing : (S[state], state+1)