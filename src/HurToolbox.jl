module HurToolbox

# greet() = print("Hello World!")

export	HurDefineRF,
		HurDefineGeneralizedCoordinates,
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
