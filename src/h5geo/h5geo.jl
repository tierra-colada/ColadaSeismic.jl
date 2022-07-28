# Import Base functions (copied from JUDI)
import Base.*, Base./, Base.+, Base.-
import Base.copy!, Base.copy
import Base.sum, Base.ndims, Base.reshape, Base.fill!, Base.axes, Base.dotview
import Base.eltype, Base.length, Base.size, Base.iterate, Base.show, Base.display, Base.showarg
import Base.maximum, Base.minimum, Base.push!
import Base.Broadcast.broadcasted, Base.BroadcastStyle, Base.Broadcast.DefaultArrayStyle, Base.Broadcast, Base.broadcast!
import Base.getindex, Base.setindex!, Base.firstindex, Base.lastindex
import Base.similar, Base.isapprox, Base.isequal
import Base.materialize!, Base.materialize
import Base.promote_shape, Base.diff, Base.cumsum, Base.cumsum!

# h5geo includes
include("H5Seis.jl")
include("H5SeisCon.jl")
include("read_H5SeisCon.jl")
include("H5Geometry.jl")
include("H5judiVector.jl")
include("Geometry.jl")
include("auxiliaryFunctions.jl")