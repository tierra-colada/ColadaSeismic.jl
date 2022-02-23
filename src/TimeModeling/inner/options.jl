mutable struct H5ModelOptions
  model_xkey::String
  model_ykey::String
end

H5ModelOptions(;
  model_xkey="CDP_X",
  model_ykey="CDP_Y"
) = H5ModelOptions(
  model_xkey,
  model_ykey)



mutable struct H5GeomOptions
  src_xkey_min::Float64
  src_xkey_max::Float64
  src_ykey_min::Float64
  src_ykey_max::Float64
  src_xkey::String
  src_ykey::String
  src_zkey::String
  rec_xkey::String
  rec_ykey::String
  rec_zkey::String
  model_orientation::Number
end

H5GeomOptions(;
  src_xkey_min=-Inf,
  src_xkey_max=Inf,
  src_ykey_min=-Inf,
  src_ykey_max=Inf,
  src_xkey="SRCX",
  src_ykey="SRCY",
  src_zkey="SES",
  rec_xkey="GRPX",
  rec_ykey="GRPY",
  rec_zkey="RGE",
  model_orientation=0.0
) = H5GeomOptions(
  src_xkey_min,
  src_xkey_max,
  src_ykey_min,
  src_ykey_max,
  src_xkey,
  src_ykey,
  src_zkey,
  rec_xkey,
  rec_ykey,
  rec_zkey,
  model_orientation)