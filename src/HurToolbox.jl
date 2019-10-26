module HurToolbox

# greet() = print("Hello World!")

# @vars n n1 n2 n3 t

using Pkg
l = ["LinearAlgebra","SymPy"]
for item in l
    Pkg.add(item)
end

using LinearAlgebra, SymPy

@vars n n1 n2 n3

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
