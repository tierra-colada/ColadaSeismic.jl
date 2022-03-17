import JUDI.split, JUDI.subsample

export H5SeisCon

mutable struct H5SeisCon
  seis::PyCall.PyObject
  pkey::String
  pkeyvals::Array{Float64,1}
end

function H5SeisCon(;
  seis::PyCall.PyObject=nothing,
  pkey::String="SP",
  pkeyvals=nothing)
  if isnothing(seis)
    @error "Seis is Nothing"
    return
  end

  if isnothing(pkey)
    @error "pkey is Nothing"
    return
  end

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