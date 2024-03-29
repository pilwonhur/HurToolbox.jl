module HurToolbox

# greet() = print("Hello World!")

# @vars n n1 n2 n3 t

using Pkg
l = ["LinearAlgebra","SymPy"]
for item in l
    Pkg.add(item)
end

using LinearAlgebra, SymPy

@vars n n1 n2 n3 t

export	@HurInitialize,
		@HurDefineRF,
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
		HurGetRelativeDCM,
		HurCoordTriads,
		HurVectorDiff,
		HurAppendRF2Coord,
		HurCrossCoord,
		HurSetAngularVel,
		HurGetAngularVel,
		HurSimplify,
		HurDiff,
		HurUnifyTriadsCoord,
		HurMatrixVectorProduct,
		HurDot,
		HurSetSimplify,
		HurMatrixVectorProductTriads,
		HurGetMMatrix,
		HurGetCMatrix,
		HurGetGVector,
		HurSetLagrangian,
		HurGetELEquationFromLagrangian,
		HurDefineNonConservativeForces,
		HurDefineCOMPos,
		HurGetNumGlobalRF,
		HurSetAngularAcc,
		HurGetAngularAcc,
		HurGetLinearCOMVel,
		HurSetCOMAcc,
		HurGetLinearCOMAcc,
		HurSetCOMVel,
		HurKinematics,
		HurTurnOnSimplify,
		HurTurnOffSimplify,
		HurDefineMass,
		HurDefineInertia,
		HurDefineVertical,
		HurGetAngularMomentum,
		HurGetLinearMomentum,
		HurGetInertiaTensor,
		HurELEquation,
		HurGetLagrangian,
		HurGetKineticE,
		HurGetPotentialE,
		HurDefineOtherPotentialE,
		HurDefineVariableList,
		HurELInverse



include("API.jl")

end # module
