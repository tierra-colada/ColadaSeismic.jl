export read_H5SeisCon

function read_H5SeisCon(con::H5SeisCon, blocks::Array{Int,1})
  if isnothing(con.seis)
    @error "Seis is Nothing\n"
    return
  end

  if isnothing(con.pkey)
    @error "pkey is Nothing\n"
    return
  end

  if isnothing(con.pkeyvals)
    @error "pkeyvals is Nothing\n"
    return
  end

  ind = Array{Integer,1}(undef, 0)
  # data = Array{Array{Float32,2},1}(undef, length(blocks))
  data = Array{Float32,2}(undef, con.seis.getNSamp(), 0)  # without preallocating rows I'm unable to `hcat`
  keylist = [con.pkey]
  for block in blocks
    minlist = [con.pkeyvals[block]]
    maxlist = [con.pkeyvals[block]]
    TRACE, _, block_ind = con.seis.getSortedData(keylist, minlist, maxlist)
    # data[block], _, block_ind = con.seis.getSortedData(keylist, minlist, maxlist)
    ind = vcat(ind, block_ind)
    data = hcat(data, TRACE)
  end

  return H5Seis(con.seis, ind, data)
end


# function read_H5SeisCon(con::H5SeisCon, blocks::Array{Int,1})
#   if isnothing(con.seis)
#     @error "Seis is Nothing\n"
#     return
#   end

#   if isnothing(con.pkey)
#     @error "pkey is Nothing\n"
#     return
#   end

#   if isnothing(con.pkeyvals)
#     @error "pkeyvals is Nothing\n"
#     return
#   end

#   data = Array{Array{Float32,2},1}(undef, length(blocks))
#   keylist = [con.pkey]
#   for block in blocks
#     minlist = [con.pkeyvals[block]]
#     maxlist = [con.pkeyvals[block]]
#     data[block], _, _ = con.seis.getSortedData(keylist, minlist, maxlist)
#   end

#   return data
# end

# RANGES & INT
function read_H5SeisCon(con::H5SeisCon, blocks::TR) where {TR<:AbstractRange}
  read_H5SeisCon(con, Array(blocks))
end
function read_H5SeisCon(con::H5SeisCon, blocks::Integer)
  read_H5SeisCon(con, [blocks])
end
