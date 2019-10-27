# Due to programming difference, procedure used in Mathematica cannot be used in
# Julia. So, there are a slight differences.
# First, newtonian reference frame "n" should be defined manually.
# Also, time "t" should be defined manually.
# Finally, I have to run HurConstructTriadsConversion() to explicitely generate 
# the triads conversion table.

# Other miscellaneous things: [] -> () in functions.
# No direct access to the global variables. Ex) HurToolbox.HurGlobalRF 

#=
Version History
v0.2


v0.1 (10/26/2019)
Basic implementation of definitions of RF, GC, time, DCM, Conversion between triads.
These are enough for manual usage of toolbox for derivation of dynamics EOM 

@HurDefineRF,@HurDefineGeneralizedCoordinates,@HurDefineTime,HurGetIndexGlobalRF,
HurFindIndex,HurTranspose,HurDefineDCM,HurDefineDCMRelative,HurConstructTriadsConversion,
HurRotationMatrix,HurRow2Matrix,HurMatrix2Row,HurMakeSymmetricMatrix,HurUnifyTriadPool,
HurUnifyTriads,HurGetRelativeDCM,HurDiff

=#



HurGlobalTime = Sym(undef)
HurGlobalRF = Array{Sym}(undef,0)
HurGlobalGeneralizedCoordinates = Array{SymFunction}(undef,0)
HurGlobalListTriads = Array{Sym}(undef,1,3)
HurGlobalTriadsConversion = Array{Pair{Sym,Sym}}(undef,0,0)
HurGlobalDCM = Array{Sym}(undef,1,9)

# push!(HurGlobalRF,n)
HurGlobalListTriads[1,:]=[n1 n2 n3]
HurGlobalTime=t
# HurGlobalTriadsConversion=[n1=>n1 n2=>n2 n3=>n3]
HurGlobalDCM[1,:]=[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]


# HurGlobalListTriads=vcat(HurGlobalListTriads,[n1 n2 n3])

# push!: one at a time
# append! can be altogether

# Symbol: julia nativge type
# Sym:    SymPy type
# If these two are compared, they are not the same


# p=n1+n2+n3
# p(HurGlobalTriadsConversion[2,1],HurGlobalTriadsConversion[2,2],HurGlobalTriadsConversion[2,3])

macro HurDefineTime(x) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block)

	push!(tmp.args,:($(esc(x))=$(esc(symbols(x))) ))
	# push!(tmp.args,:($(esc(HurGlobalTime))=$(esc(Sym(x))) ))
	# push!(tmp.args, :(push!($(esc(HurGlobalTime)),$(esc(x))))  )
	global HurGlobalTime=Sym(string(x))
	return tmp;
end

macro HurDefineRF(x...) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block)
	for xx in x
		# avoid redefinition of RFs
		ind=HurGetIndexGlobalRF(xx)	# HurGlobalRF contains Sym, whereas xx is Symbol
		if ind==0
			push!(tmp.args,:($(esc(xx))=$(esc(symbols(xx))) ))
			push!(tmp.args, :(push!($(esc(HurGlobalRF)),$(esc(xx))))  )
			
			global HurGlobalDCM=vcat(HurGlobalDCM,[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0])
			HurGlobalDCM[1,:]=[1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0]

			global HurGlobalListTriads=vcat(HurGlobalListTriads,[n1 n2 n3])
			n,=size(HurGlobalListTriads)
			for i=1:3
				tempp=string(xx)*string(i)
				temp=Symbol(tempp)
				push!(tmp.args,:($(esc(temp))=$(esc(symbols(temp)))))
				HurGlobalListTriads[n-1,i]=Sym(tempp)
			end
		end
	end
	# push!(tmp.args, :(print(length($(esc(HurGlobalRF))))  ))
	return tmp;
end

macro HurDefineGeneralizedCoordinates(x...) # ... is the way to handle tuples in the argument.
	tmp=Expr(:block)
	for xx in x
		# avoid redefinition of RFs
		push!(tmp.args,:($(esc(xx))=$(esc(SymFunction(string(xx)))) ))
		push!(tmp.args, :(push!($(esc(HurGlobalGeneralizedCoordinates)),$(esc(xx))))  )
		
		# HurGlobalGeneralizedCoordinates=gcs;
		# HurGlobalELEquation=Table[0,{i,ngcs}];
 		# HurGlobalNonConservativeForces=Table[0,{i,ngcs}];
  		# HurGlobalConstrainedELEquation=Table[0,{i,ngcs}];			
	end
	# push!(tmp.args, :(print(length($(esc(HurGlobalRF))))  ))
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
	HurDefineDCMRelative(rf, HurGlobalRF[1], dcm)	
end

function HurDefineDCMRelative(rf1, rf2, dcm)
	rot=HurMatrix2Row(HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:])*dcm)
	HurGlobalDCM[HurGetIndexGlobalRF(rf1),:]=[simplify(rot[1]) simplify(rot[2]) simplify(rot[3]) simplify(rot[4]) simplify(rot[5]) simplify(rot[6]) simplify(rot[7]) simplify(rot[8]) simplify(rot[9])]
end

function HurConstructTriadsConversion()
	n=length(HurGlobalRF)
	global HurGlobalTriadsConversion=Array{Pair{Sym,Sym}}(undef,n,n*3)
	
	for i=1:n
		for j=1:n
			Rot=HurUnifyTriadPool(HurGlobalRF[i], HurGlobalRF[j])
			# tuple((temp[i] for i=1:3)...)
			tmpTriads=Rot*HurGlobalListTriads[i,:]
			HurGlobalTriadsConversion[i,3*(j-1)+1]=HurGlobalListTriads[j,1]=>tmpTriads[1]
			HurGlobalTriadsConversion[i,3*(j-1)+2]=HurGlobalListTriads[j,2]=>tmpTriads[2]
			HurGlobalTriadsConversion[i,3*(j-1)+3]=HurGlobalListTriads[j,3]=>tmpTriads[3]
		end
	end
end

function HurRotationMatrix(ang,dir)
	dirv=findall(x->x!=0,dir)
	if length(dirv)==1 # simple rotation
		if dirv[1]==1 # rotation about x axis
			rot=[1.0 0.0 0.0;0.0 cos(ang) -sin(ang);0.0 sin(ang) cos(ang)]
		elseif dirv[1]==2 # rotation about y axis
			rot=[cos(ang) 0.0 sin(ang);0.0 1.0 0.0;-sin(ang) 0.0 cos(ang)]
		else # rotation about z axis
			rot=[cos(ang) -sin(ang) 0.0;sin(ang) cos(ang) 0.0;0.0 0.0 1.0]
		end
	else
		rot=[1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0]
	end
	return rot
end
	


function HurRow2Matrix(v)
	return [v[1] v[2] v[3];v[4] v[5] v[6];v[7] v[8] v[9]]
end

function HurMatrix2Row(Mat)
	return [Mat[1,1] Mat[1,2] Mat[1,3] Mat[2,1] Mat[2,2] Mat[2,3] Mat[3,1] Mat[3,2] Mat[3,3]] #reshape(Mat',1,:)
end

function HurMakeSymmetricMatrix(v)
	return [v[1] v[2] v[3];v[2] v[4] v[5];v[3] v[5] v[6]]
end

function HurUnifyTriadPool(rf1, rf2)
	rot1=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf1),:])
	rot2=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:])
	return HurTranspose(rot2)*rot1
end

function HurUnifyTriads(v,rf)
	n=length(HurGlobalRF)
	temp= v((HurGlobalTriadsConversion[HurGetIndexGlobalRF(rf),i] for i=1:3*n)...)	
	coef1=diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),1])
	coef2=diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),2])
	coef3=diff(temp,HurGlobalListTriads[HurGetIndexGlobalRF(rf),3])
	return coef1*HurGlobalListTriads[HurGetIndexGlobalRF(rf),1]+coef2*HurGlobalListTriads[HurGetIndexGlobalRF(rf),2]+coef3*HurGlobalListTriads[HurGetIndexGlobalRF(rf),3]
end

function HurGetRelativeDCM(rf1,rf2)
	rot1=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf1),:])
	rot2=HurRow2Matrix(HurGlobalDCM[HurGetIndexGlobalRF(rf2),:])
	r=HurTranspose(rot2)*rot1;
	return [simplify(r[1,1]) simplify(r[1,2]) simplify(r[1,3]);simplify(r[2,1]) simplify(r[2,2]) simplify(r[2,3]);simplify(r[3,1]) simplify(r[3,2]) simplify(r[3,3]);]
end

function HurDiff()	# for testing only
	diff(HurGlobalDCM[3,1],HurGlobalTime)
end
