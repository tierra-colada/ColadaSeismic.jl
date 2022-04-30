mutable struct H5Seis
  seis::PyCall.PyObject
  # types must be compatible with h5geo seis IO data
  ind::Array{UInt64,1}
  data::Array{Float32,2}
end
