module HurToolbox

# greet() = print("Hello World!")

# @vars n n1 n2 n3 t

using Pkg
l = ["LinearAlgebra","SymPy"]
for item in l
    Pkg.add(item)
end

using LinearAlgebra, SymPy

# global n=Sym("n")
# global n1=Sym("n1")
# global n2=Sym("n2")
# global n3=Sym("n3")
# global t=Sym("t")

export	@HurDefineRF,
		@HurDefineGeneralizedCoordinates,
		@HurDefineTime,
		HurGetIndexGlobalRF,
		HurFindIndex,
		HurTranspose,
		HurDefineDCM,
		HurDefineDCMRelative,
		HurConstructTriadsConversion,
		HurRotationMatrix,
		HurRow2Matrix,
		HurMatrix2Row,
		HurMakeSymmetricMatrix,
		HurUnifyTriadPool,
		HurUnifyTriads,
		HurGetRelativeDCM

include("API.jl")

end # module
