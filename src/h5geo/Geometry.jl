import JUDI.Geometry

# Ovewritten to be able to rotate coordinates
# Load geometry from out-of-core Geometry container and apply rotation
function Geometry(geometry::GeometryOOC)
  nsrc = length(geometry.container)

  # read either source or receiver geometry
  if geometry.key=="source"
      params = ["SourceX","SourceY",geometry.segy_depth_key,"dt","ns"]
      gt = Float32
  elseif geometry.key=="receiver"
      params = ["GroupX","GroupY",geometry.segy_depth_key,"dt","ns"]
      gt = Array{Float32, 1}
  else
      throw("Specified keyword not supported")
  end
  xloc = Array{gt, 1}(undef, nsrc)
  yloc = Array{gt, 1}(undef, nsrc)
  zloc = Array{gt, 1}(undef, nsrc)
  dt = Array{Float32}(undef, nsrc); nt = Array{Integer}(undef, nsrc); t = Array{Float32}(undef, nsrc)

  for j=1:nsrc

      header = read_con_headers(geometry.container[j], params, 1)
      if geometry.key=="source"
          xloc[j] = convert(gt, get_header(header, params[1])[1])
          yloc[j] = convert(gt, get_header(header, params[2])[1])
          zloc[j] = abs.(convert(gt,get_header(header, params[3])[1]))
      else
          xloc[j] = convert(gt, get_header(header, params[1]))
          yloc[j] = convert(gt, get_header(header, params[2]))
          zloc[j] = abs.(convert(gt, get_header(header, params[3])))
      end

      # Rotate coordinates
      xloc[j] = xloc[j] .- model_origin_x
      yloc[j] = yloc[j] .- model_origin_y
      xloc[j], yloc[j] = rotate_arrays(xloc[j], yloc[j], -model_orientation)
      xloc[j] = xloc[j] .+ model_origin_x
      yloc[j] = yloc[j] .+ model_origin_y

      dt[j] = get_header(header, params[4])[1]/1f3
      nt[j] = convert(Integer, get_header(header, params[5])[1])
      t[j] =  (nt[j]-1)*dt[j]
  end
  if geometry.key == "source"
      xloc = convertToCell(xloc)
      yloc = convertToCell(yloc)
      zloc = convertToCell(zloc)
  end
  return GeometryIC(xloc,yloc,zloc,dt,nt,t)
end