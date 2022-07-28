import JUDI.get_computational_nt, JUDI.time_modeling, JUDI.write_shot_record

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

function write_shot_record(srcGeometry::GeometryIC, srcData, recGeometry::GeometryIC, recData, options)
  nTrc = length(recGeometry.xloc[1])
  nTrc < 1 && return
  
  if save_as == SEGY
    q = judiVector(srcGeometry, srcData)
    d = judiVector(recGeometry, recData)
    pos = [srcGeometry.xloc[1][1], srcGeometry.yloc[1][1],  srcGeometry.zloc[1][1]]
    pos = join(["_"*string(trunc(p; digits=2)) for p in pos])
    file = join([string(options.file_name), pos,".segy"])
    block_out = judiVector_to_SeisBlock(d, q)
    segy_write(join([options.file_path,"/",file]), block_out)
    container = scan_file(join([options.file_path,"/",file]),
                          ["GroupX", "GroupY", "dt", "SourceSurfaceElevation", "RecGroupElevation"];
                          chunksize=256)
    return container
  elseif save_as == SaveAs::H5SEIS
    from_trace = 0
    src_x = repeat(srcGeometry.xloc[1], nTrc)
    src_y = repeat(srcGeometry.yloc[1], nTrc)
    src_z = repeat(srcGeometry.zloc[1], nTrc)
    rec_x = recGeometry.xloc[1]
    rec_y = recGeometry.yloc[1] == nTrc ? recGeometry.yloc[1] : repeat(recGeometry.yloc[1], nTrc)
    rec_z = recGeometry.zloc[1]
    if isnothing(seis_out)
      p = h5geo.SeisParam()
      p.spatialReference = spatial_reference
      p.lengthUnits = "m"
      p.temporalUnits = "ms"
      p.domain = h5geo.Domain.TWT
      p.dataType = h5geo.SeisDataType.PRESTACK
      p.surveyType = survey_type
      p.nSamp = recGeometry.nt
      p.nTrc = nTrc
      p.trcChunk = nTrc

      # splitext removes extension
      global seis_out = seis_out_cnt.createSeis(splitext(options.file_name)[1], p, h5geo.CreationType.CREATE_OR_OVERWRITE)
      if isnothing(seis_out)
        @error "Unable to create or overwrite H5Seis to store calculated data"
        return
      end
      seis_out.setSampRate(recGeometry.dt)
    else
      from_trace = seis_out.getNTrc()
      seis_out.setNTrc(from_trace+nTrc)
    end

    status = seis_out.writeTrace(recData, from_trace, 0)
    if !status
      @error "Unable to write traces"
      return
    end

    status &= seis_out.writeTraceHeader("SRCX", src_x, from_trace)
    status &= seis_out.writeTraceHeader("SRCY", src_y, from_trace)
    status &= seis_out.writeTraceHeader("SES", src_z, from_trace)
    status &= seis_out.writeTraceHeader("GRPX", rec_x, from_trace)
    status &= seis_out.writeTraceHeader("GRPY", rec_y, from_trace)
    status &= seis_out.writeTraceHeader("RGE", rec_z, from_trace)
    if !status
      @error "Unable to write trace headers"
      return
    end

    return H5SeisCon(seis=seis_out)
  end
end
