mutable struct H5Seis
  seis::PyCall.PyObject
  ind::Array{Integer,1}
  data::Array{Float32,2}
end
