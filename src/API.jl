# Due to programming difference, procedure used in Mathematica cannot be used in
# Julia. So, there are a slight differences.
# First, Newtonian reference frame "n" should be defined manually.
# Also, time "t" should be defined manually.
# Finally, I have to run HurConstructTriadsConversion() to explicitely generate 
# the triads conversion table.

# Other miscellaneous things: [] -> () in functions.
# No direct access to the global variables. Ex) HurToolbox.HurGlobalRF 

#=
Version History
v0.2 (10/26/2019)
The following were added.
HurVectorDiff,HurAppendRF2Coord,HurCrossCoord,HurSetAngularVel,HurGetAngularVel,
HurSimplify,HurDiff,HurUnifyTriadsCoord

v0.1 (10/26/2019)
Basic implementation of definitions of RF, GC, time, DCM, Conversion between triads.
These are enough for manual usage of toolbox for derivation of dynamics EOM 

@HurDefineRF,@HurDefineGeneralizedCoordinates,@HurDefineTime,HurGetIndexGlobalRF,
HurFindIndex,HurTranspose,HurDefineDCM,HurDefineDCMRelative,HurConstructTriadsConversion,
HurRotationMatrix,HurRow2Matrix,HurMatrix2Row,HurMakeSymmetricMatrix,HurUnifyTriadPool,
HurUnifyTriads,HurGetRelativeDCM,HurDiff

=#

HurGlobalG = Sym(undef)
HurGlobalTime = Sym(undef)
HurGlobalRF = Array{Sym}(undef,1)
HurGlobalGeneralizedCoordinates = Array{SymFunction}(undef,0)
HurGlobalListTriads = Array{Sym}(undef,1,3)
HurGlobalTriadsConversion = Array{Pair{Sym,Sym}}(undef,0,0)
HurGlobalDCM = Array{Sym}(undef,1,9)
HurGlobalAngularVel=Array{Sym}(undef,1)
HurGlobalAngularAcc=Array{Sym}(undef,1)
HurGlobalSimplify = true
HurGlobalELEquation=Array{Sym}(undef,1)
HurGlobalLagrangian=Array{Sym}(undef,1)
HurGlobalPotentialE=Array{Sym}(undef,1)
HurGlobalOtherPotentialE=Array{Sym}(undef,1)
HurGlobalKineticE=Array{Sym}(undef,1)
HurGlobalMMatrix = Array{Sym}(undef,0,0)
HurGlobalCMatrix = Array{Sym}(undef,0,0)
HurGlobalGVector = Array{Sym}(undef,0)
HurGlobalRayleighDissipationE = Array{Sym}(undef,0)
HurGlobalNonConservativeForces=Array{Sym}(undef,1)
HurGlobalCOMPos = Array{Sym}(undef,0)
HurGlobalCOMVel	= Array{Sym}(undef,0)
HurGlobalCOMAcc = Array{Sym}(undef,0)
HurGlobalMass = Array{Sym}(undef,0)
HurGlobalInertia = Array{Sym}(undef,1,6)
HurGlobalVertical = Sym(undef)
HurGlobalLinearMomentum = Array{Sym}(undef,0)
HurGlobalAngularMomentum = Array{Sym}(undef,0)
HurGlobalVariableList = Array{Sym}(undef,0)

# push!(HurGlobalRF,n)
HurGlobalListTriads[1,:]=[n1 n2 n3];
HurGlobalTime=t;
# HurGlobalTriadsConversion=[n1=>n1 n2=>n2 n3=>n3]
HurGlobalDCM[1,:]=[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0];
HurGlobalInertia[1,:]=[1.0 0.0 0.0 1.0 0.0 1.0];
HurGlobalAngularVel[1]=0;

# HurGlobalListTriads=vcat(HurGlobalListTriads,[n1 n2 n3])

# push!: one at a time
# append! can be altogether

# Symbol: julia nativge type
# Sym:    SymPy type
# If these two are compared, they are not the same


# p=n1+n2+n3
# p(HurGlobalTriadsConversion[2,1],HurGlobalTriadsConversion[2,2],HurGlobalTriadsConversion[2,3])
macro HurDefineTime(x) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block);
	@vars t
	# push!(tmp.args,:($(esc(x))=$(esc(symbols(x))) ));
	push!(tmp.args,:($(esc(Symbol("t")))=$(esc(t) )));
	# push!(tmp.args,:($(esc(HurGlobalTime))=$(esc(Sym(x))) ))
	# push!(tmp.args, :(push!($(esc(HurGlobalTime)),$(esc(x))))  )
	# global HurGlobalTime=Sym(string(x));
	global HurGlobalTime=t;
	return tmp;
end

macro HurInitialize()
	tmp=Expr(:block);
	@vars g t n n1 n2 n3
	push!(tmp.args,:($(esc(Symbol("t")))=$(esc(t) )));
	push!(tmp.args,:($(esc(Symbol("g")))=$(esc(g) )));
	push!(tmp.args,:($(esc(Symbol("n")))=$(esc(n) )));
	push!(tmp.args,:($(esc(Symbol("n1")))=$(esc(n1) )));
	push!(tmp.args,:($(esc(Symbol("n2")))=$(esc(n2) )));
	push!(tmp.args,:($(esc(Symbol("n3")))=$(esc(n3) )));

	global HurGlobalG = g
	global HurGlobalTime = t
	# global HurGlobalTime = t

	global HurGlobalRF[1] = n
	global HurGlobalListTriads[1,:]=[n1 n2 n3]
	global HurGlobalDCM[1,:]=[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]
	return tmp;
end

# function HurDefineRF(rfs...)
# 	n=length(rfs)

# 	global HurGlobalRF = Array{Sym}(undef,n)

# 	for rf in rfs
# 		if HurGetIndexGlobalRF(rf)==0
# 			push!(HurGlobalRF,symbols(rf))
# 		end

# 	end
# end

macro HurDefineRF(x...) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block)
	for xx in x
		# avoid redefinition of RFs
		ind=HurGetIndexGlobalRF(xx)	# HurGlobalRF contains Sym, whereas xx is Symbol
		if ind==0
			push!(tmp.args,:($(esc(xx))=$(esc(symbols(xx))) ));
			push!(tmp.args, :(push!($(esc(HurGlobalRF)),$(esc(xx))))  );
			
			global HurGlobalDCM=vcat(HurGlobalDCM,[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]);
			
			global HurGlobalListTriads=vcat(HurGlobalListTriads,[n1 n2 n3]);
			n,=size(HurGlobalListTriads);
			for i=1:3
				tempp=string(xx)*string(i);
				temp=Symbol(tempp);
				push!(tmp.args,:($(esc(temp))=$(esc(symbols(temp)))));

				HurGlobalListTriads[n,i]=Sym(tempp);
			end
		end
	end

	m=length(x)+1

	global HurGlobalRayleighDissipationE = Array{Sym}(undef,m)
	global HurGlobalLagrangian = Array{Sym}(undef,m)
	global HurGlobalPotentialE = Array{Sym}(undef,m)
	global HurGlobalOtherPotentialE = Array{Sym}(undef,m)
	global HurGlobalKineticE = Array{Sym}(undef,m)
	global HurGlobalCOMPos = Array{Sym}(undef,m)
	global HurGlobalCOMVel = Array{Sym}(undef,m)
	global HurGlobalCOMAcc = Array{Sym}(undef,m)
	global HurGlobalAngularVel = Array{Sym}(undef,m)
	global HurGlobalAngularAcc = Array{Sym}(undef,m)
	global HurGlobalMass = Array{Sym}(undef,m)
	global HurGlobalInertia = Array{Sym}(undef,m,6)
	global HurGlobalLinearMomentum = Array{Sym}(undef,m)
	global HurGlobalAngularMomentum = Array{Sym}(undef,m)

	return tmp;
end

# macro HurDefineRF1(x...) # ... is the way to handle tuples in the argument.
# 	tmp=Expr(:block)
# 	for xx in x
# 		# avoid redefinition of RFs
# 		ind=HurGetIndexGlobalRF(xx)	# HurGlobalRF contains Sym, whereas xx is Symbol
# 		if ind==0
# 			push!(tmp.args,:($(esc(xx))=$(esc(symbols(xx))) ));
# 			push!(tmp.args, :(push!($(esc(HurGlobalRF)),$(esc(xx))))  );
			
# 			global HurGlobalDCM=vcat(HurGlobalDCM,[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]);
# 			HurGlobalDCM[1,:]=[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0];

# 			# push!(HurGlobalAngularVel,0);

# 			global HurGlobalListTriads=vcat(HurGlobalListTriads,[n1 n2 n3]);
# 			n,=size(HurGlobalListTriads);
# 			for i=1:3
# 				tempp=string(xx)*string(i);
# 				temp=Symbol(tempp);
# 				push!(tmp.args,:($(esc(temp))=$(esc(symbols(temp)))));
# 				HurGlobalListTriads[n-1,i]=Sym(tempp);
# 			end
# 		end
# 	end
# 	# push!(tmp.args, :(print(length($(esc(HurGlobalRF))))  ))
# 	return tmp;
# end

macro HurDefineGeneralizedCoordinates(gcs...) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block)
	n=length(gcs)
	for gc in gcs
		# avoid redefinition of RFs
		push!(tmp.args,:($(esc(gc))=$(esc(SymFunction(string(gc)))) ));
		push!(tmp.args, :(push!($(esc(HurGlobalGeneralizedCoordinates)),$(esc(gc)))) );
	end

	global HurGlobalELEquation = Array{Sym}(undef,n)
	global HurGlobalMMatrix = Array{Sym}(undef,n,n)
	global HurGlobalCMatrix = Array{Sym}(undef,n,n)
	global HurGlobalGVector = Array{Sym}(undef,n)
	global HurGlobalNonConservativeForces = Array{Sym}(undef,n)
	
	return tmp;
end

function HurGetIndexGlobalRF(rf)
	return HurFindIndex(HurGlobalRF,rf)
end

function HurFindIndex(Arr,expr)
	index=0
	count=0
	for val in Arr
		count+=1
		if val==Sym(string(expr))
			index=count
			break
		end
	end
	return index
	# ind=findall(x->x==expr,Arr)
	# return ind[1]
end

function HurTranspose(mat)
	n,m=size(mat)
	matp=typeof(mat)(undef,m,n)
	for i=1:n
		for j=1:n
			matp[i,j]=mat[j,i]
		end
	end
	return matp
end

function HurDefineDCM(rf, dcm)
	HurDefineDCMRelative(rf, HurGlobalRF[1], dcm);	
end

function HurDefineDCMRelative(rf1, rf2, dcm)
	rot=HurMatrix2Row(HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:])*dcm);


	if HurGlobalSimplify
		HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]=[simplify(rot[1]) simplify(rot[2]) simplify(rot[3]) simplify(rot[4]) simplify(rot[5]) simplify(rot[6]) simplify(rot[7]) simplify(rot[8]) simplify(rot[9])];
	else
		HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]=[rot[1] rot[2] rot[3] rot[4] rot[5] rot[6] rot[7] rot[8] rot[9]];
	end
end

function HurConstructTriadsConversion()
	n=length(HurGlobalRF);
	global HurGlobalTriadsConversion=Array{Pair{Sym,Sym}}(undef,n,n*3);
	global HurGlobalListTriads
	global HurGlobalDCM

	global HurGlobalSimplify
	tempsimplify=HurGlobalSimplify
	HurGlobalSimplify=true
	for i=1:n
		for j=1:n
			Rot=HurSimplify(HurUnifyTriadPool(HurGlobalRF[i], HurGlobalRF[j]));
			# tuple((temp[i] for i=1:3)...)
			tmpTriads=Rot*HurGlobalListTriads[i,:];
			HurGlobalTriadsConversion[i,3*(j-1)+1]=HurGlobalListTriads[j,1]=>tmpTriads[1];
			HurGlobalTriadsConversion[i,3*(j-1)+2]=HurGlobalListTriads[j,2]=>tmpTriads[2];
			HurGlobalTriadsConversion[i,3*(j-1)+3]=HurGlobalListTriads[j,3]=>tmpTriads[3];
		end
	end
	HurGlobalSimplify=tempsimplify

	m,=size(HurGlobalListTriads)
	if m>n
		HurGlobalListTriads=HurGlobalListTriads[setdiff(1:end,m),:]
	end

	m,=size(HurGlobalDCM)
	k,=size(HurGlobalGeneralizedCoordinates)
	if m>n
		HurGlobalDCM=HurGlobalDCM[setdiff(1:end,m),:]
	end
	m,=size(HurGlobalGeneralizedCoordinates)
	global HurGlobalLagrangian
	[HurGlobalLagrangian[i]=0 for i=1:n]
	global HurGlobalRayleighDissipationE
	[HurGlobalRayleighDissipationE[i]=0 for i=1:n]
	global HurGlobalPotentialE
	[HurGlobalPotentialE[i]=0 for i=1:n]
	global HurGlobalOtherPotentialE
	[HurGlobalOtherPotentialE[i]=0 for i=1:n]
	global HurGlobalKineticE
	[HurGlobalKineticE[i]=0 for i=1:n]
	global HurGlobalCOMPos
	[HurGlobalCOMPos[i]=0 for i=1:n]
	global HurGlobalCOMVel
	[HurGlobalCOMVel[i]=0 for i=1:n]
	global HurGlobalCOMAcc
	[HurGlobalCOMAcc[i]=0 for i=1:n]
	global HurGlobalAngularVel
	[HurGlobalAngularVel[i]=0 for i=1:n]
	global HurGlobalAngularAcc
	[HurGlobalAngularAcc[i]=0 for i=1:n]
	global HurGlobalMass
	[HurGlobalMass[i]=0 for i=1:n]
	global HurGlobalInertia
	[HurGlobalInertia[j,i]=0 for j in 1:n, i in 1:6]
	global HurGlobalVertical
	HurGlobalVertical=HurGlobalListTriads[1,3]
	global HurGlobalLinearMomentum
	[HurGlobalLinearMomentum[i]=0 for i=1:n]
	global HurGlobalAngularMomentum
	[HurGlobalAngularMomentum[i]=0 for i=1:n]

	global HurGlobalNonConservativeForces
	[HurGlobalNonConservativeForces[i]=0 for i=1:k]
	global HurGlobalELEquation
	[HurGlobalELEquation[i]=0 for i=1:k]	

	# m,=size(HurGlobalAngularVel)
	# if m>n
	# 	HurGlobalAngularVel=HurGlobalAngularVel[setdiff(1:end,m)]
	# end
end

function HurRotationMatrix(ang,dir)
	dirv=findall(x->x!=0,dir)
	if length(dirv)==1 # simple rotation
		if dirv[1]==1 # rotation about x axis
			rot=[1.0 0.0 0.0;0.0 cos(ang) -sin(ang);0.0 sin(ang) cos(ang)];
		elseif dirv[1]==2 # rotation about y axis
			rot=[cos(ang) 0.0 sin(ang);0.0 1.0 0.0;-sin(ang) 0.0 cos(ang)];
		else # rotation about z axis
			rot=[cos(ang) -sin(ang) 0.0;sin(ang) cos(ang) 0.0;0.0 0.0 1.0];
		end
	else
		rot=[1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0];
	end
	return rot;
end
	
function HurRow2Matrix(v)
	return [v[1] v[2] v[3];v[4] v[5] v[6];v[7] v[8] v[9]];
end

function HurMatrix2Row(Mat)
	return [Mat[1,1] Mat[1,2] Mat[1,3] Mat[2,1] Mat[2,2] Mat[2,3] Mat[3,1] Mat[3,2] Mat[3,3]]; #reshape(Mat',1,:)
end

function HurMakeSymmetricMatrix(v)
	return [v[1] v[2] v[3];v[2] v[4] v[5];v[3] v[5] v[6]];
end

function HurUnifyTriadPool(rf1, rf2)
	rot1=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]);
	rot2=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:]);
	return HurTranspose(rot2)*rot1;
end

function HurUnifyTriads(v,rf)
	n=length(HurGlobalRF);
	temp= v((HurGlobalTriadsConversion[HurGetIndexGlobalRF(rf),i] for i=1:3*n)...);	
	coef1=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),1]));
	coef2=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),2]));
	coef3=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),3]));
	temp=coef1*HurGlobalListTriads[HurGetIndexGlobalRF(rf),1]+coef2*HurGlobalListTriads[HurGetIndexGlobalRF(rf),2]+coef3*HurGlobalListTriads[HurGetIndexGlobalRF(rf),3]
	return temp
end

function HurUnifyTriadsCoord(v,rf)
	n=length(HurGlobalRF);
	temp= v((HurGlobalTriadsConversion[HurGetIndexGlobalRF(rf),i] for i=1:3*n)...);	
	coef1=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),1]));
	coef2=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),2]));
	coef3=HurSimplify(diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),3]));
	return [coef1,coef2,coef3,rf];
end

function HurGetRelativeDCM(rf1,rf2)
	rot1=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]);
	rot2=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:]);
	r=HurTranspose(rot2)*rot1;
	if HurGlobalSimplify
		return [simplify(r[1,1]) simplify(r[1,2]) simplify(r[1,3]);simplify(r[2,1]) simplify(r[2,2]) simplify(r[2,3]);simplify(r[3,1]) simplify(r[3,2]) simplify(r[3,3])];
	else
		return [r[1,1] r[1,2] r[1,3];r[2,1] r[2,2] r[2,3];r[3,1] r[3,2] r[3,3]];
	end
end

function HurDiff(mat,var)	# for testing only
	matp=copy(mat)
	n=length(mat)
	for i=1:n
		if n==1
			matp=diff(mat,var)
		else
			matp[i]=diff(mat[i],var)
		end
		
	end
	return matp
end

function HurSimplify(mat)	# for testing only
	matp=copy(mat)
	n=length(mat)
	if HurGlobalSimplify
		for i=1:n
			if typeof(mat[i]) == Sym
				if n==1
					matp=simplify(mat)
				else
					matp[i]=simplify(mat[i])
				end
			end
		end
	end
	return matp
end

function HurGetAngularVel(rf1,rf2)
	rot5=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]);
	Sw=HurDiff(rot5, HurGlobalTime)*HurTranspose(rot5);
	ww=(-Sw[2,3])*HurGlobalListTriads[1,1]+(Sw[1,3])*HurGlobalListTriads[1,2]+(-Sw[1,2])*HurGlobalListTriads[1,3];
	ww=HurSimplify(ww);
	www=HurUnifyTriads(ww,rf2);
	if HurGlobalSimplify
		www=HurSimplify(www)
	end
	HurSetAngularVel(rf1,www);
	return www;
end

function HurGetAngularAcc(rf1,rf2)
	alpha=HurVectorDiff( HurGlobalAngularVel[HurGetIndexGlobalRF(rf1)], HurGlobalRF[1],rf2)
	# alpha=HurUnifyTriads(HurCoordTriads(HurAppendRF2Coord[D[HurUnifyTriadsCoord[HurGetAngularVel[rf1,HurGlobalRF[[1]] ],HurGlobalRF[[1]] ][[1;;3]],Global`t],HurGlobalRF[[1]] ] ),rf2)
	# alpha1=If[HurGlobalSimplify, Simplify[alpha], alpha];
	if HurGlobalSimplify
		alpha=HurSimplify(alpha)
	end
  	HurSetAngularAcc(rf1,alpha)
  	alpha
end
  
function HurSetAngularVel(rf,w)
	global HurGlobalAngularVel
	HurGlobalAngularVel[HurGetIndexGlobalRF(rf)]=w;
end

function HurSetAngularAcc(rf,alpha)
	global HurGlobalAngularAcc
	HurGlobalAngularAcc[HurGetIndexGlobalRF(rf)] = alpha;
end

function HurCrossCoord(v1, v2, rf)
	temp=cross(HurUnifyTriadsCoord(v1,rf)[1:3],HurUnifyTriadsCoord(v2,rf)[1:3]);
	temp=HurSimplify(temp)
	return HurAppendRF2Coord(temp,rf)
end

function HurAppendRF2Coord(coord, rf)
	push!(coord,rf);
	return coord;
end

function HurVectorDiff(v,rf1,rf2)
	df2dvdt=HurDiff(HurUnifyTriadsCoord(v,rf2)[1:3],HurGlobalTime)
	www=HurUnifyTriads(HurGetAngularVel(rf2,rf2)-HurGetAngularVel(rf1,rf2),rf2);
	wcrossv=HurCrossCoord(www,v,rf2);
	df1dvdt=df2dvdt+wcrossv[1:3];

	return HurCoordTriads(HurAppendRF2Coord(df1dvdt,rf2));
	# return df1dvdt
end

function HurCoordTriads(v...)
	n=length(v)
	if n==1
		vec=sum(v[1][1:3].*HurGlobalListTriads[HurGetIndexGlobalRF(v[1][4]),:]);
	elseif n==2
		vec=sum(v[1][:].*HurGlobalListTriads[HurGetIndexGlobalRF(v[2]),:]);
	end
	if HurGlobalSimplify
		vec=HurSimplify(vec)
	end
	return vec;
end

function HurMatrixVectorProduct(mat,vec,rf)
	temp=mat*HurUnifyTriadsCoord(vec,rf)[1:3]
	if HurGlobalSimplify
		temp=HurSimplify(temp)
	end 
    return temp 
end

function HurMatrixVectorProductTriads(mat,vec,rf) 
    return HurCoordTriads(HurMatrixVectorProduct(mat,vec,rf),rf)
end

function HurDot(v1, v2)
	tempD=sum(HurUnifyTriadsCoord(v1, HurGlobalRF[1])[1:3].*HurUnifyTriadsCoord(v2, HurGlobalRF[1])[1:3]);
	if HurGlobalSimplify
		return simplify(tempD)
	else
		return tempD
	end
end

function HurSetSimplify(flag)
	global HurGlobalSimplify=flag
end

function HurTurnOnSimplify()
	global HurGlobalSimplify=true
end

function HurTurnOffSimplify()
	global HurGlobalSimplify=false
end

function HurGetMMatrix()
	n=length(HurGlobalGeneralizedCoordinates);
	global HurGlobalMMatrix
	for i=1:n
		for j=1:n
			qdd=diff(HurGlobalGeneralizedCoordinates[j](HurGlobalTime),HurGlobalTime,HurGlobalTime);
			HurGlobalMMatrix[i,j]=diff(HurGlobalELEquation[i],qdd);
		end
	end
	if HurGlobalSimplify
		HurGlobalMMatrix=HurSimplify(HurGlobalMMatrix)
	end
	return HurGlobalMMatrix
end

function HurGetCMatrix()
	gcs=HurGlobalGeneralizedCoordinates
	ti=HurGlobalTime
	n=length(gcs);
	global HurGlobalCMatrix
	for k=1:n
		for j=1:n
			HurGlobalCMatrix[k,j]=sum([1/2*(diff(HurGlobalMMatrix[k,j], gcs[i](ti))+diff(HurGlobalMMatrix[k,i], gcs[j](ti))-diff(HurGlobalMMatrix[i,j], gcs[k](ti)))*diff(gcs[i](ti),ti) for i=1:n])
		end
	end
	if HurGlobalSimplify
		HurGlobalCMatrix=HurSimplify(HurGlobalCMatrix)
	end
	return HurGlobalCMatrix
end

function HurGetGVector()
	gcs=HurGlobalGeneralizedCoordinates
	ti=HurGlobalTime
	n=length(gcs);
	global HurGlobalGVector=[diff(sum(HurGlobalPotentialE+HurGlobalOtherPotentialE),gcs[i](ti)) for i=1:n];
	if HurGlobalSimplify
		HurGlobalGVector=HurSimplify(HurGlobalGVector)
	end
	return HurGlobalGVector
end

function HurDefineMass(rf, m)
	global HurGlobalMass[HurGetIndexGlobalRF(rf)] = m
end

function HurDefineInertia(rf, II)
	global HurGlobalInertia
	[HurGlobalInertia[HurGetIndexGlobalRF(rf),i] = II[i] for i in 1:6]
end

function HurDefineVertical(v)
	global HurGlobalVertical=v
end

function HurGetELEquationFromLagrangian()
	gcs=HurGlobalGeneralizedCoordinates;
	ngcs=length(gcs);
	ti=HurGlobalTime
	L=sum(HurGlobalLagrangian);
	DE=sum(HurGlobalRayleighDissipationE);
	global HurGlobalELEquation
	for i=1:ngcs
		HurGlobalELEquation[i]=diff(diff(L,diff(gcs[i](ti),ti)),ti)-diff(L,gcs[i](ti))+diff(DE,diff(gcs[i](ti),ti))-HurGlobalNonConservativeForces[i];
	end
	if HurGlobalSimplify
		HurGlobalELEquation=HurSimplify(HurGlobalELEquation)
	end
	return HurGlobalELEquation;
end

function HurELEquation()
	nrfs=HurGetNumGlobalRF()
	if nrfs==1
		error("There is only one reference frame!")
	else
		for i in 2:nrfs
			if HurGlobalMass[i]!=0
				HurGetLinearMomentum(HurGlobalRF[i],HurGlobalRF[i]);
        		HurGetAngularMomentum(HurGlobalRF[i],HurGlobalRF[i]);        
			end
		end
		HurGetLagrangian(HurGlobalRF[2:nrfs]);
    	# HurGetELEquation(HurGlobalGeneralizedCoordinates);
    	HurGetELEquationFromLagrangian()
    	HurGetMMatrix();
    	HurGetCMatrix();
    	HurGetGVector();
	end
	return HurGlobalELEquation
end
# Ib = HurMakeSymmetricMatrix([mb*kb^2, 0, 0, 0, 0, mb*kb^2])

function HurDefineVariableList(var...)
	global HurGlobalVariableList=[var[i] for i=1:length(var)]
end

function HurELInverse()
	n=length(HurGlobalELEquation)
	# tempequations=Flatten[ Table[HurGlobalELEquation[[i]]==0,{i,n} ] ];	
	if length(HurGlobalVariableList)==0
		variables=[diff(HurGlobalGeneralizedCoordinates[i](HurGlobalTime),HurGlobalTime,HurGlobalTime) for i in 1:length(HurGlobalGeneralizedCoordinates)]
	else
		variables=HurGlobalVariableList
	end
	temp=solve(HurGlobalELEquation,variables)
	if HurGlobalSimplify
		temp1=[simplify(temp[variables[i]]) for i in 1:length(variables)]
	else
		temp1=[temp[variables[i]] for i in 1:length(variables)]
	end
	return temp1
end

function HurDefineOtherPotentialE(rf, pe)
	HurGlobalOtherPotentialE[HurGetIndexGlobalRF(rf)] = pe
	return sum(HurGlobalOtherPotentialE)
end

function HurDefineRayleighDissipationE(rf, de)
	HurGlobalRayleighDissipationE[HurGetIndexGlobalRF(rf)] = de
	return sum(HurGlobalRayleighDissipationE)
end

function HurGetLagrangian(rf...)
	HurGetKineticE(rf);
	HurGetPotentialE(rf);
	global HurGlobalLagrangian=HurGlobalKineticE-HurGlobalPotentialE-HurGlobalOtherPotentialE
	return sum(HurGlobalLagrangian)
end

function HurGetKineticE(rf...)
	# narg=length(rf)
	global HurGlobalKineticE
	rf=HurGlobalRF
	narg=length(rf)
	for i in 2:narg
		KE=1/2*HurGlobalMass[HurGetIndexGlobalRF(rf[i])]*HurDot(HurGlobalCOMVel[HurGetIndexGlobalRF(rf[i])],HurGlobalCOMVel[HurGetIndexGlobalRF(rf[i])])+1/2*HurDot(HurGlobalAngularVel[HurGetIndexGlobalRF(rf[i])],HurGlobalAngularMomentum[HurGetIndexGlobalRF(rf[i])])
		if HurGlobalSimplify
			KE=simplify(KE)
		end
		HurGlobalKineticE[HurGetIndexGlobalRF(rf[i])] = KE
	end
	return sum(HurGlobalKineticE)
end

function HurGetPotentialE(rf...)
	# narg=length(rf)
	global HurGlobalPotentialE
	rf=HurGlobalRF
	narg=length(rf)
	for i in 2:narg
		PE=HurGlobalG*HurGlobalMass[HurGetIndexGlobalRF(rf[i])]*HurDot(HurGlobalCOMPos[HurGetIndexGlobalRF(rf[i])],HurGlobalVertical[1]);
    	if HurGlobalSimplify
			PE=simplify(PE)
		end
    	HurGlobalPotentialE[HurGetIndexGlobalRF(rf[i])] = PE;
	end
	return sum(HurGlobalPotentialE)
end

function HurGetInertiaTensor(rf)
	IItemp=HurGlobalInertia[HurGetIndexGlobalRF(rf),:]
	II=	HurMakeSymmetricMatrix(IItemp)
	return II
end

function HurGetAngularMomentum(rf1, rf2)
	II=HurGetInertiaTensor(rf1)
	H = HurMatrixVectorProductTriads(II, HurGetAngularVel(rf1,rf1), rf1)
	H1=HurUnifyTriads(H,rf2)
	if HurGlobalSimplify
		H1=simplify(H1)
	end
	global HurGlobalAngularMomentum
	HurGlobalAngularMomentum[HurGetIndexGlobalRF(rf1)] = H1
  	return H1
end

function HurGetLinearMomentum(rf1, rf2)
	v=HurGlobalCOMVel[HurGetIndexGlobalRF(rf1)]
	m=HurGlobalMass[HurGetIndexGlobalRF(rf1)]
	ll=HurUnifyTriads(m*v,rf2)
	
	if HurGlobalSimplify
		ll=simplify(ll)
	end
	global HurGlobalLinearMomentum
	HurGlobalLinearMomentum[HurGetIndexGlobalRF(rf1)] = ll
  	return ll
end

function HurSetLagrangian(L,rf)
	global HurGlobalLagrangian
	HurGlobalLagrangian[HurGetIndexGlobalRF(rf)]=L
end

function HurDefineNonConservativeForces(f...)
	n=length(f)
	m,=size(HurGlobalGeneralizedCoordinates)
	if n!=m
		error("The number of nonconservative forces you entered is not the same as the number of generalized coordinates!")
	end
	global HurGlobalNonConservativeForces
	[HurGlobalNonConservativeForces[i]=f[i] for i=1:n]
end

function HurDefineCOMPos(rf,v)
	global HurGlobalCOMPos
	HurGlobalCOMPos[HurGetIndexGlobalRF(rf)] = v
end

function HurGetNumGlobalRF()
	return length(HurGlobalRF)
end

function HurKinematics()
	nrfs=HurGetNumGlobalRF()
	if nrfs==1
		error("There is only one reference frame!")
	else
		for i in 2:nrfs
			HurGetAngularVel(HurGlobalRF[i] , HurGlobalRF[i])
			HurGetAngularAcc(HurGlobalRF[i] , HurGlobalRF[i])
			HurGetLinearCOMVel(HurGlobalRF[i] , HurGlobalRF[i])
			HurGetLinearCOMAcc(HurGlobalRF[i] , HurGlobalRF[i])
		end
	end
end

function HurGetLinearCOMVel(rf1,rf2)
	tempx=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[1]
	tempy=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[2]
	tempz=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[3]
	tempvx=diff(tempx,HurGlobalTime)
	tempvy=diff(tempy,HurGlobalTime)
	tempvz=diff(tempz,HurGlobalTime)

	v=HurUnifyTriads(HurCoordTriads(HurAppendRF2Coord([tempvx,tempvy,tempvz],HurGlobalRF[1])),rf2);
	if HurGlobalSimplify
		v=HurSimplify(v)
	end
	HurSetCOMVel(rf1,v);
	return v
end

function HurSetCOMVel(rf,v)
	global HurGlobalCOMVel
	HurGlobalCOMVel[HurGetIndexGlobalRF(rf)] = v
end

function HurGetLinearCOMAcc(rf1,rf2)
	tempx=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[1]
	tempy=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[2]
	tempz=HurUnifyTriadsCoord(HurGlobalCOMPos[HurGetIndexGlobalRF(rf1)],HurGlobalRF[1])[3]
	tempax=diff(tempx,HurGlobalTime,HurGlobalTime)
	tempay=diff(tempy,HurGlobalTime,HurGlobalTime)
	tempaz=diff(tempz,HurGlobalTime,HurGlobalTime)

	acc=HurUnifyTriads(HurCoordTriads(HurAppendRF2Coord([tempax,tempay,tempaz],HurGlobalRF[1])),rf2);
	# acc1=If[HurGlobalSimplify, Simplify[acc], acc];
	if HurGlobalSimplify
		acc=HurSimplify(acc)
	end
	HurSetCOMAcc(rf1,acc);
	return acc
end

function HurSetCOMAcc(rf,acc)
	global HurGlobalCOMAcc
	HurGlobalCOMAcc[HurGetIndexGlobalRF(rf)] = acc
end