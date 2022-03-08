module B

using A

include("TimeModeling/B.jl")

end # module B


using A
using B

a()