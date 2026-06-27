(* ::Package:: *)

current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;


(* ::Subsection:: *)
(*Functions*)


FamilyInfo=Association[];


JSector[vs_]:=If[#>0,1,0]&/@vs;
JProp[vs_]:=Count[vs,_?(Positive)];
JDots[vs_]:=If[#>0,#-1,0]&/@vs//Total;
SortFIs[ints_]:=SortBy[ints,{-JProp[#[[-1]]],-JDots[#[[-1]]],Min[#[[-1]]]}&];


BFFamilyDefine[family_,loopdenominators_,loopmomenta_,externalmomenta_,scalarproducts_]:=Module[{branch},
	FamilyInfo[family]=Association[];
	BFFamily = family;
	BFPropagators = loopdenominators;
	BFLoop = loopmomenta; 
	BFIndepLeg = externalmomenta;
	BFSPToSTU=scalarproducts;
	BFBranch=Values[KeySortBy[PositionIndex[Intersection[#,BFLoop]&/@Variables/@(BFPropagators/.BFIndepLeg->0)],Length]];
	FeynParameters=SortBy[Flatten@Table[BFBranch[[i,j]]->z[i][j],{i,Length@BFBranch},{j,Length@BFBranch[[i]]}],First]//Values;
	FamilyInfo[family,"ABC"]=EvaluateABC[BFPropagators,FeynParameters];
	FamilyInfo[family,"UF"]=EvaluateUF[BFPropagators,FeynParameters];
	FamilyInfo[family,"SR"]=EvaluateSR[family];
	
	genIBPSolution[family];
	genCornerRelation[family];
	KMasterInfo[family]=Complement[KI[family,0,#]&/@Keys[Select[BFSectorType[family],#===1&]]//SortFIs//Reverse,Keys[BFCornerRelation[family]]];
];


EvaluateABC[gpds_,para_]:=Module[{denominator,\[Lambda],A,B,C},
	denominator=gpds . para/.Thread[BFLoop->\[Lambda]*BFLoop];
	
	If[Exponent[denominator,\[Lambda]]=!=2,
		Print["Incorrect input for propagators: not quadratic propagator."];
		Abort[];
	];
	
	{C,B,A}=CoefficientList[denominator,\[Lambda]];
	B=1/2 Coefficient[B,#]&/@BFLoop//Expand;
	A=Table[
		If[i==j, Coefficient[A,BFLoop[[i]]^2],
			1/2 Coefficient[A,BFLoop[[i]]BFLoop[[j]]]
		]
	,{i,Length@BFLoop},{j,Length@BFLoop}];
	{A,B,C}
];


EvaluateUF[gpds_,para_]:=Module[{A,B,C,U,F,G},
	{A,B,C}=EvaluateABC[gpds,para]/.z[i_][j_]:>X[i]z[i][j];
	A=A/.Table[z[i][1]->1-Total[z[i]/@Range[2,Length@BFBranch[[i]]]],{i,Length@BFBranch}]//Expand;
	C=C/.z[i_][j_]:>z[i][j]Total[z[i]/@Range[Length@BFBranch[[i]]]];
	U = Det[A]//Expand;
	If[U===0,Return@{0,0}];(*singular case*)
	F = Expand[B . Adjugate[A] . B - U* C]/.BFSPToSTU;
	G = F+U Sum[z[1][j],{j,Length@BFBranch[[1]]}]^2;
	{U,F,G}//Expand
];


EvaluateSR[family_]:=Module[{R,S,addmatrix},
	R=Table[If[i===j,2Coefficient[FamilyInfo[family,"UF"][[3]],FeynParameters[[i]]^2],Coefficient[FamilyInfo[family,"UF"][[3]],FeynParameters[[i]]FeynParameters[[j]]]]
		,{i,Length@FeynParameters},{j,Length@FeynParameters}];
	addmatrix=Table[Boole[Head[#][[1]]===i]&/@FeynParameters,{i,Length@BFBranch}];
	S=ArrayFlatten[{{0IdentityMatrix[Length@BFBranch],addmatrix},{addmatrix//Transpose,R}}];
	{S,R}
];


genIBPSolutionSector[fam_,sector_]:=Module[{n=Length@sector,secpos,down,up,rescVar,repVar,F,S,T,R,detS,invS,cz,C,derPos,rankS,nullS,gi,sub,sol,index,vec,relation,needpos},
	secpos=Flatten[Position[sector,1]];
	{S,R}=FamilyInfo[fam,"SR"];
	S=S[[Join[Range[Length@BFBranch],secpos+Length@BFBranch],Join[Range[Length@BFBranch],secpos+Length@BFBranch]]];
	detS=Det[S]//Together;
	vec=Table[Slot[secpos[[i]]]KI[fam,Slot[n+1]+2,Array[Slot,n]+UnitVector[n,secpos[[i]]]],{i,Length@secpos}];
	needpos=Table[Count[sector[[BFBranch[[i]]]],1],{i,Length@BFBranch}];
	relation=S[[;;Length@BFBranch,-Length@secpos;;]] . vec+ConstantArray[KI[fam,Slot[n+1]+2,Array[Slot,n]],Length@BFBranch];
	relation=Join[relation,Join@@Table[
		Table[(S[[Length@BFBranch+Total[needpos[[;;i-1]]]+j,Length@BFBranch+1;;]]-S[[Length@BFBranch+Total[needpos[[;;i-1]]]+j+1,Length@BFBranch+1;;]]) . vec-
		(KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[Total[needpos[[;;i-1]]]+j]]]]-KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[Total[needpos[[;;i-1]]]+j+1]]]])
		,{j,needpos[[i]]-1}]
	,{i,Length@BFBranch}]];
	
	If[detS=!=0,
		invS=Inverse[S]//Together;
		cz=invS . Join[ConstantArray[1,Length@BFBranch],ConstantArray[0,Length@secpos]]//Together;
		C=Total[cz[[;;Length@BFBranch]]]//Together;
		If[C=!=0,
			BFSectorType[fam,sector]=1;
			up=(C KI[fam,Slot[n+1]+2,Array[Slot,n]]-Sum[cz[[Length@BFBranch+i]]KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}])/
				(Total[Array[Slot,n]]-Length@BFBranch-FBId-Slot[n+1]);
			down=((Total[Array[Slot,n]]-Length@BFBranch-FBId-Slot[n+1]+2)KI[fam,Slot[n+1]-2,Array[Slot,n]]+
				Sum[cz[[Length@BFBranch+i]]KI[fam,Slot[n+1]-2,Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}])/C;
			BFDimShiftInfo[fam,sector]={down,up};
			BFIBPIdentity[fam,sector]={};
			(*BFEtaDiffInfo[fam,sector]=((2(FBIdelta+Slot[1])-Total[sector]-Length@BFBranch)KI[fam,Slot[1],sector]+
				Sum[cz[[Length@BFBranch+i]]KI[fam,Slot[1]-1,sector-UnitVector[n,secpos[[i]]]],{i,Length@secpos}])/(2eta-C);*)
			sol=invS . Join[Table[-KI[fam,Slot[n+1],Array[Slot,n]],Length@BFBranch],
				Table[KI[fam,Slot[n+1]-2,Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]];
			sol=Association@Table[secpos[[i]]->Collect[sol[[i+Length@BFBranch]]/Slot[secpos[[i]]],_KI,Together],{i,Length@secpos}];
			Return[sol]
			,
			BFSectorType[fam,sector]=2;
			BFIBPIdentity[fam,sector]=relation;
			(*BFEtaDiffInfo[fam,sector]=((2(FBIdelta+Slot[1])-Total[sector]-Length@BFBranch)KI[fam,Slot[1],sector]+
				Sum[cz[[Length@BFBranch+i]]KI[fam,Slot[1]-1,sector-UnitVector[n,secpos[[i]]]],{i,Length@secpos}])/(2eta);*)
			sol=Sum[cz[[Length@BFBranch+i]]KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]/
				(FBId+Slot[n+1]+Length@BFBranch-Total[Array[Slot,n]]);
			Return[sol]
		];
	];
	
	rankS=MatrixRank[S];
	nullS=NullSpace[S]//Together;
	
	If[rankS<Length@S-1,
		BFSectorType[fam,sector]=4;
		relation=Join[relation,Table[
			C=Together[Total[nullS[[m,;;Length@BFBranch]]]];
			C KI[fam,Slot[n+1]+2,Array[Slot,n]]-
				Sum[nullS[[m,Length@BFBranch+i]]KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]
		,{m,Length@nullS}]];
		BFIBPIdentity[fam,sector]=relation;
		sol=Select[nullS,Together[Total[#[[;;Length@BFBranch]]]]=!=0&];
		If[Length@sol>0,
			sol=sol[[1]];
			C=Together[Total[sol[[;;Length@BFBranch]]]];
			sol=Sum[sol[[Length@BFBranch+i]]KI[fam,Slot[n+1]-2,Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]/C
			,
			sol=nullS[[1,Length@BFBranch+1;;]];
			index=FirstPosition[sol,SelectFirst[sol,#=!=0&]][[1]];
			sol=-Sum[sol[[i]]KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[i]]]+UnitVector[n,secpos[[index]]]],{i,DeleteCases[Range[Length@secpos],index]}]/sol[[index]]
		];
		Return[sol];
		,
		sol=nullS[[1]];
		C=Together[Total[sol[[;;Length@BFBranch]]]];
		If[C=!=0,
			BFSectorType[fam,sector]=3;
			BFIBPIdentity[fam,sector]=relation;
			(*BFEtaDiffInfo[fam,sector]=Sum[sol[[Length@BFBranch+i]]KI[fam,Slot[1]-1,sector-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]/(-C);*)
			sol=Sum[sol[[Length@BFBranch+i]]KI[fam,Slot[n+1]-2,Array[Slot,n]-UnitVector[n,secpos[[i]]]],{i,Length@secpos}]/C
			,
			BFSectorType[fam,sector]=4;
			BFIBPIdentity[fam,sector]=relation;
			sol=nullS[[1,Length@BFBranch+1;;]];
			index=FirstPosition[sol,SelectFirst[sol,#=!=0&]][[1]];
			sol=-Sum[sol[[i]]KI[fam,Slot[n+1],Array[Slot,n]-UnitVector[n,secpos[[i]]]+UnitVector[n,secpos[[index]]]],{i,DeleteCases[Range[Length@secpos],index]}]/sol[[index]]
		];
		Return[sol]
	];

];


genSeedSector[sector_,deg_]:=Module[{unit,final={sector},new={sector}},
	unit=Table[If[sector[[i]]===1,UnitVector[Length@sector,i],Nothing],{i,Length@sector}];
	Table[
		new=Join@@Table[Plus[new[[j]],#]&/@unit,{j,Length@new}];
		final=Join[final,new]
	,{n,deg}];
	final
];


genIBP[fam_,deg_]:=Module[{sectors,seeds},
	sectors=BFSectorType[fam]//Keys;
	Join@@Table[
		seeds=genSeedSector[sectors[[i]],deg];
		Join@@Table[Function[Evaluate[BFIBPIdentity[fam,sectors[[i]]]]]@@Append[seeds[[j]],0],{j,Length@seeds}]
	,{i,Length@sectors}]
];


BFIBPSolution=Association[];


BFIBPIdentity=Association[];


BFSectorType=Association[];


BFDimShiftInfo=Association[];


(*BFEtaDiffInfo=Association[];*)


genIBPSolution[fam_]:=Module[{allsectors},
	allsectors=Tuples[{1,0},Length@BFPropagators];
	allsectors=Select[allsectors,Length@DeleteDuplicates[Head/@Pick[FeynParameters,#,1]](*>=Length@BFLoop*)===Length@BFBranch&];
	BFSectorType[fam]=Association[];
	BFDimShiftInfo[fam]=Association[];
	BFIBPIdentity[fam]=Association[];
	(*BFEtaDiffInfo[fam]=Association[];*)
	BFIBPSolution[fam]=Association[Table[allsectors[[i]]->genIBPSolutionSector[fam,allsectors[[i]]],{i,Length@allsectors}]];
];


ClearAll[WJraw];
WJraw[fam_,m_,vs_]:=(*WJraw[fam,m,vs]=*)If[!KeyExistsQ[BFSectorType[fam],JSector[vs]],0,
	Switch[BFSectorType[fam,JSector[vs]],
		1,If[Max@vs<=1,If[m===0,KI[fam,m,vs],
			If[m>0,Function[Evaluate[BFDimShiftInfo[fam,JSector[vs]][[1]]]]@@Append[vs,m]/.KI:>WJraw,
			Function[Evaluate[BFDimShiftInfo[fam,JSector[vs]][[2]]]]@@Append[vs,m]/.KI:>WJraw]],
		Function[Evaluate[BFIBPSolution[fam,JSector[vs],FirstPosition[vs,Max[vs]][[1]]]]]@@Append[vs-UnitVector[Length@vs,FirstPosition[vs,Max[vs]][[1]]],m]/.KI:>WJraw],
		2|3|4,Function[Evaluate[BFIBPSolution[fam,JSector[vs]]]]@@Append[vs,m]/.KI:>WJraw,
		_,0
	]
]


BFCornerRelation=Association[];


genCornerRelation[fam_]:=Module[{ibps,ints,sol,pos},
	ibps=genIBP[fam,1];
	ibps=DeleteDuplicates@DeleteCases[Collect[ibps/.KI:>WJraw,_KI,Together],0];
	If[Length@ibps===0,BFCornerRelation[fam]={},
		ints=SortFIs[Cases[ibps,_KI,Infinity]//DeleteDuplicates];
		sol=CoefficientArrays[ibps,ints]//Last//Normal//RowReduce;
		BFCornerRelation[fam]=Table[
			pos=Position[sol[[i]],1];
			If[Length@pos>0,pos=pos[[1,1]];ints[[pos]]->ints[[pos]]-sol[[i]] . ints,Nothing]
		,{i,Length@sol}]
	]
];


ClearAll[WJ];
WJ[fam_,m_,vs_]:=(*WJ[fam,m,vs]=*)WJraw[fam,m,vs]/.BFCornerRelation[fam];


KMasterInfo=Association[];


genBranchDE[fam_,vars_]:=Block[{master=KMasterInfo[fam],exp,R,secpos,r},
	R=FamilyInfo[fam,"SR"][[2]];
	If[Head@vars===List,
		Table[
			Table[
				secpos=Flatten@Position[master[[ii,-1]],1];
				r=D[R[[secpos,secpos]],var];
				exp=Sum[If[i===j,-r[[i,j]]KI[fam,2,master[[ii,-1]]+2UnitVector[Length@master[[ii,-1]],secpos[[i]]]],
					-r[[i,j]]KI[fam,2,master[[ii,-1]]+UnitVector[Length@master[[ii,-1]],secpos[[i]]]+UnitVector[Length@master[[ii,-1]],secpos[[j]]]]]
					,{i,Length@secpos},{j,i,Length@secpos}];
				exp=exp/.KI:>WJ;
				If[exp===0,ConstantArray[0,Length@master],CoefficientArrays[exp,master]//Last//Normal]
			,{ii,Length@master}]
		,{var,vars}]
		,
		Table[
			secpos=Flatten@Position[master[[ii,-1]],1];
			r=D[R[[secpos,secpos]],vars];
			exp=Sum[If[i===j,-r[[i,j]]KI[fam,2,master[[ii,-1]]+2UnitVector[Length@master[[ii,-1]],secpos[[i]]]],
				-r[[i,j]]KI[fam,2,master[[ii,-1]]+UnitVector[Length@master[[ii,-1]],secpos[[i]]]+UnitVector[Length@master[[ii,-1]],secpos[[j]]]]]
				,{i,Length@secpos},{j,i,Length@secpos}];
			exp=exp/.KI:>WJ;
			If[exp===0,ConstantArray[0,Length@master],CoefficientArrays[exp,master]//Last//Normal]
		,{ii,Length@master}]
	]
];


(*genEtaDE[fam_]:=Module[{fameta=StringJoin[ToString[fam],"eta"],master,topnu},
	topnu=Total[KMasterInfo[fam][[-1,-1]]];
	master=KI[fam,0,#]&/@Keys[Select[BFSectorType[fam],#=!=4&]]//SortFIs//Reverse;
	master=Table[KI[fam,Total[master[[i,-1]]]-topnu,master[[i,-1]]],{i,Length@master}];
	
	{master/.fam:>fameta,Table[
		CoefficientArrays[Function[Evaluate[BFEtaDiffInfo[fam,master[[i,-1]]]]]@@{master[[i,2]]},master]//Last//Normal
	,{i,Length@master}]}
];*)


(*genInfBoundary[fam_]:=Module[{master,topnu},
	master=KI[fam,0,#]&/@Keys[Select[BFSectorType[fam],#=!=4&]]//SortFIs//Reverse;
	topnu=Total[KMasterInfo[fam][[-1,-1]]];
	Table[
		FBIdelta-topnu->(-1)^(Total[master[[i,-1]]])Gamma[topnu-FBIdelta]
	,{i,Length@master}]
];*)


(* ::Section:: *)
(*test*)


(* ::Subsection::Closed:: *)
(*massive quark self energy*)


loopmom={l1,l2};
extmom={p};
pdlist={(l1)^2-msq,(l1+p)^2-msq,(l2)^2-msq,(l2-p)^2-msq,(l1+l2)^2-msq}/.msq->1;
spsRep={p^2->11};
(*branch={2,2,1}*)(*propagator number of each branch*)


BFFamilyDefine["se1",pdlist,loopmom,extmom,spsRep];


{diffGX1,diffGX2,diffGX3}=genBranchDE["se1",{X[1],X[2],X[3]}];


diffGX1[[-1]]//Together


Put[{KMasterInfo["se1"],{diffGX1,diffGX2,diffGX3}},FileNameJoin[{current,"se1bde"}]]


(*reduction*)
(*0 means the dimension is D+0, FBId is D*)
Collect[KI["se1",0,{2,1,1,1,1}]/.KI:>WJ,_KI,Together]


(* ::Subsection:: *)
(*three points massless*)


loopmom={l1,l2};
extmom={p1,p2};
pdlist={(l1)^2,(l1-p1)^2,(l1-p1-p2)^2,(l2-p1-p2)^2,(l2)^2,(l1-l2)^2};
spsRep={p1^2->0,p2^2->0,p1*p2->1/2};
(*branch={2,2,1}*)(*propagator number of each branch*)


BFFamilyDefine["rk",pdlist,loopmom,extmom,spsRep];


BFSectorType


KMasterInfo["rk"]


diffGX1=genBranchDE["rk",X[3]];


diffGX1//Together


diffGX1[[-5]]//Together


diffGX12[[-5]]


diffGX1-diffGX12//Together


(* ::Subsection:: *)
(*cross double box*)


loopmom={l1,l2};
extmom={p1,p2,p3};
pdlist={(l1)^2,(l1+p1)^2,(l1+p1+p2)^2,(l2)^2-mt2,(l2-p3)^2-mt2,(l1+l2)^2-mt2,(l1+l2+p1+p2-p3)^2-mt2}/.mt2->23;
spsRep={p1^2->0,p2^2->0,p3^2->mh2,p1*p2->s/2,p1*p3->(mh2-t)/2,p2*p3->(s+t-mh2)/2}/.{s->101,t->-43,mh2->12};


BFFamilyDefine["cdb",pdlist,loopmom,extmom,spsRep];


BFSectorType["cdb"]


KMasterInfo["cdb"]


diffGX1=genBranchDE["cdb",X[3]];


diffGX1/.{X[1]->1/11,X[2]->3/17,X[3]->7/23}//Together
