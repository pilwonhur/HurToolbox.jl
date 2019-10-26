module HurToolbox

# greet() = print("Hello World!")

@vars n n1 n2 n3 t

export	n,n1,n2,n3,t,
		@HurDefineRF,
		@HurDefineGeneralizedCoordinates,
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

using Pkg
l = ["LinearAlgebra","SymPy"]
for item in l
    Pkg.add(item)
end

using LinearAlgebra, SymPy

include("API.jl")

end # module
