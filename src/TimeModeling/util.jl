
export rotate_Nx2_array, rotate_2xN_array, rotate_arrays, h5geo2judiTraceHeaderName, judi2h5geoTraceHeaderName
@enum H5PhPType begin
  VELOCITY = 1
  DENSITY = 2
  QUALITYFACTOR = 3
  EPSILON = 4
  DELTA = 5
  THETA = 6
  PHI = 7
end

@enum SaveAs begin
  SEGY = 1
  H5SEIS = 2
end

function rotate_Nx2_array(xy, ϕ)
  return xy * [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
end

function rotate_2xN_array(xy, ϕ)
  return [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)] * xy
end

function rotate_arrays(x, y, ϕ)
  x_calc = x*cos(ϕ)+y*sin(ϕ)
  y_calc = -x*sin(ϕ)+y*cos(ϕ)
  return x_calc, y_calc
end

function h5geo2judiTraceHeaderName(h5geoHdrName::String)
  h5geo = pyimport("h5geopy._h5geo")
  _, h5geoShortHeaderNames = h5geo.getTraceHeaderNames()

  ind = findall(x->x==h5geoHdrName, h5geoShortHeaderNames)
  if length(ind) < 1
    return
  end

  judiHeaderNames = judiTraceHeaders()
  if ind[1] > length(judiHeaderNames)
    return
  end

  return judiHeaderNames[ind[1]]
end

function judi2h5geoTraceHeaderName(judiHdrName::String)
  h5geo = pyimport("h5geopy._h5geo")
  judiHeaderNames = judiTraceHeaders()

  ind = findall(x->x == judiHdrName, judiHeaderNames)
  if length(ind) < 1
    return
  end

  _, h5geoShortHeaderNames = h5geo.getTraceHeaderNames()
  if ind[1] > length(h5geoShortHeaderNames)
    return
  end

  return h5geoShortHeaderNames[ind[1]]
end

function judiTraceHeaders()
  return [
    "TraceNumWithinLine"             
    "TraceNumWithinFile"             
    "FieldRecord"                    
    "TraceNumber"                    
    "EnergySourcePoint"              
    "CDP"                            
    "CDPTrace"                       
    "TraceIDCode"                    
    "NSummedTraces"                  
    "NStackedTraces"                 
    "DataUse"                        
    "Offset"                         
    "RecGroupElevation"              
    "SourceSurfaceElevation"         
    "SourceDepth"                    
    "RecDatumElevation"              
    "SourceDatumElevation"           
    "SourceWaterDepth"               
    "GroupWaterDepth"                
    "ElevationScalar"                
    "RecSourceScalar"                
    "SourceX"                        
    "SourceY"                        
    "GroupX"                         
    "GroupY"                         
    "CoordUnits"                     
    "WeatheringVelocity"             
    "SubWeatheringVelocity"          
    "UpholeTimeSource"               
    "UpholeTimeGroup"                
    "StaticCorrectionSource"         
    "StaticCorrectionGroup"          
    "TotalStaticApplied"             
    "LagTimeA"                       
    "LagTimeB"                       
    "DelayRecordingTime"             
    "MuteTimeStart"                  
    "MuteTimeEnd"                    
    "ns"                             
    "dt"                             
    "GainType"                       
    "InstrumentGainConstant"         
    "InstrumntInitialGain"           
    "Correlated"                     
    "SweepFrequencyStart"            
    "SweepFrequencyEnd"              
    "SweepLength"                    
    "SweepType"                      
    "SweepTraceTaperLengthStart"     
    "SweepTraceTaperLengthEnd"       
    "TaperType"                      
    "AliasFilterFrequency"           
    "AliasFilterSlope"               
    "NotchFilterFrequency"           
    "NotchFilterSlope"               
    "LowCutFrequency"                
    "HighCutFrequency"               
    "LowCutSlope"                    
    "HighCutSlope"                   
    "Year"                           
    "DayOfYear"                      
    "HourOfDay"                      
    "MinuteOfHour"                   
    "SecondOfMinute"                 
    "TimeCode"                       
    "TraceWeightingFactor"           
    "GeophoneGroupNumberRoll"        
    "GeophoneGroupNumberTraceStart"  
    "GeophoneGroupNumberTraceEnd"    
    "GapSize"                        
    "OverTravel"                     
    "CDPX"                           
    "CDPY"                           
    "Inline3D"                       
    "Crossline3D"                    
    "ShotPoint"                      
    "ShotPointScalar"                
    "TraceValueMeasurmentUnit"   
    "TransductionConstnatMantissa"  
    "TransductionConstantPower"     
    "TransductionUnit"              
    "TraceIdentifier"               
    "ScalarTraceHeader"             
    "SourceType"                    
    "SourceEnergyDirectionMantissa" 
    "SourceEnergyDirectionExponent" 
    "SourceMeasurmentMantissa"      
    "SourceMeasurementExponent"     
    "SourceMeasurmentUnit"          
    "Unassigned1"                   
    "Unassigned2"                               
    ]
end


# function h5geoTraceHeaders()
#   return ["SEQWL", "SEQWR", "FFID", "TRCFLD", "SP", "CDP", "TRCNUM", "TRCID", "NVST", "NHST", "DU", "DSREG", "RGE", "SES", "SDBS", "DERG", "DES", "WDS", "WGD", "SAED", "SAC", "SRCX", "SRCY", "GRPX", "GRPY", "UNITS", "WVEL", "SVEL", "UTSRC", "UTGRP", "SECSCOR", "GRPSCOR", "TSA", "LAGTA", "LAGTB", "DELRECT", "MTSTART", "MTEND", "NSMP", "SI", "GTFI", "IG", "IGC", "CORREL", "SFSTART", "SFEND", "SLEN", "STYP", "SSTRLS", "SSTLE", "TTYP", "AFF", "AFS", "NFF", "NFS", "LOCF", "HOCF", "LOCS", "HICS", "YEAR", "DAY", "HOUR", "MINUTE", "SCE", "TMBS", "TWF", "GGNSW", "GGN1ST", "GGNLST", "GAPSZ", "OAWT", "CDP_X", "CDP_Y", "INLINE", "XLINE", "SPN", "SPS", "TVMU"]
# end
