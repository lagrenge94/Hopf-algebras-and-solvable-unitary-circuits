(* ::Package:: *)

BeginPackage[ "HATools`"]
 LoadHA::usage="Input name of HA, return HA Data {mul,\[CapitalDelta], \[Eta], \[Epsilon], S, R, Rb}.";
 HAData::usage="HAData[dA, Mul, CoMul, \[Epsilon], S] computes the Hopf agebra data {MulTable, CoMulTable, \[Epsilon]Table, STable} from the linear maps Mul, CoMul, \[Epsilon], S.";
 GroupHA::usage="GroupHA[MulTable,InvTable] cosntructs the group Hopf algebra using group multiplication table and inverses.";
 DirectSum::usage = "DirectSum[\[Rho]1,\[Rho]2] Return the direct sum of Reps \[Rho]1,\[Rho]2.";
 DirectSums::usage = "DirectSum[\[Rho]1,\[Rho]2,...] Return the direct sum of Reps \[Rho]1,\[Rho]2...";
 CheckRep::usage="CheckRep[\[Rho],mul,\[Eta]] Check representation \[Rho] of an algebra (mul,\[Eta]).";
 CheckCorep::usage="CheckCorep[v,\[CapitalDelta],\[Epsilon]]Check corepresentation v of a coalgebra (\[CapitalDelta],\[Epsilon])." ;
 RepToMap::usage="Transform representation \[Rho] to a linear map/matrix, for storage.";
 RepToCorep::usage="RepToCorep[\[Rho],R] use the representation \[Rho] and the R-matrix to construct corepresentation \!\(\*SubscriptBox[\(v\), \(ij\)]\)=(id\[TensorProduct] \!\(\*SubscriptBox[\(\[Rho]\), \(ij\)]\))R.";
 \[Rho]Reg::usage = "\[Rho]Reg[mul] gives the regular representation ";
 \[Rho]Regr::usage = "\[Rho]Regr[mul] gives the regular (right) representation ";
 vRegl::usage = "vRegl[\[CapitalDelta]] gives the regular left corepresentation ";
 vRegr::usage = "vRegr[\[CapitalDelta]] gives the regular right corepresentation ";
 VerifyAlgebra::usage="VerifyAlgebra[mul,\[Eta]] verify unital associative algebra."
 VerifyCoalgebra::usage="VerifyCoalgebra[\[CapitalDelta],\[Epsilon]] verify counital coassociative coalgebra."
 VerifyHAAxioms::usage="VerifyHAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon],S] verify the Hopf algebra axioms. ";
 VerifyWHAAxioms::usage="VerifyWHAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon],S] verify the weak Hopf algebra axioms. ";
 VerifyWeakAntipode::usage=" VerifyWeakAntipode[mul,\[CapitalDelta],\[Eta],\[Epsilon],S] verify the weak antipode axioms.";
 VerifyBAAxioms::usage="VerifyBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]] verify the bialgebra axioms. ";
 VerifyWBAAxioms::usage="VerifyWBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]] verify the weak bialgebra axioms. ";
 VerifyPreBAAxioms::usage="VerifyPreBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]] verify the prebialgebra axioms. ";
 VerifyQuasitriangularity::usage="VerifyQuasitriangularity[mul,\[CapitalDelta],\[Eta],\[Epsilon],R,Rb] verify the quasitriangular bialgebra axioms.";
 VerifyHA::usage="Verify HA axioms."
 VerifyQHA::usage="Verify quasitriangular HA axioms."
 Dual::usage="Construct the dual (co)algebra of a (algebra)coalgebra.";
 DualHA::usage="Construct the dual Hopf algebra";
 DrinfeldDouble::usage="DrinfeldDouble[mul,\[CapitalDelta],\[Eta],\[Epsilon],S] Construct the Drinfeld Double {mulD,\[CapitalDelta]D,\[Eta]D,\[Epsilon]D,SD,RD,RDb} of the Hopf algebra {mul,\[CapitalDelta],\[Eta],\[Epsilon],S}.";
 DHaction\[Rho]H::usage="DHaction\[Rho]H[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,\[Rho]H] constructs a representation (generally not irreducible) of D(H) using the structure data of H and a representation \[Rho]H of H."
 DHaction\[Rho]HSub::usage="DHaction\[Rho]HSub[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,\[Rho]H,V] constructs a subrep of the above rep, where V is a subspace of H^* 
	invariant under the left regular action of H^* and the adjoint action of H."
 DGRep::usage="DGRep[CC,CG,\[Rho]CG] constructs a rep of D(G) from a conjugacy class CC of G and a representation \[Rho]CG of the group Hopf algebra CG."
 DHaction\[Rho]Hd::usage="DHaction\[Rho]Hd[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,\[Rho]Hd] constructs a representation (generally not irreducible) of D(H) using the structure data of H and a representation \[Rho]Hd of the dual of H."
 UniversalRX::usage="UniversalRX[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,Rb] construct the universal RX tensors {UL,UR,VL,VR} (elements of H\[TensorProduct]End(H) ).";
 VerifyUniversalRX::usage="VerifyUniversalRX[mul,R,UL,UR,VL,VR] verify the universal RX equations."
 DrinfeldElement::usage="DrinfeldElement[mul,S,R] computes the Drinfeld element and its inverse {u,uinv}.";
 HaarIntegral::usage=" HaarIntegral[\[CapitalDelta]] compute the Haar Integral of a semi-simple Hopf algebra.";
 RmatQYBE::usage="RmatQYBE[R,\[Rho]1,\[Rho]2] gives the R-matrix satisfying QYBE R12 R13 R23=R23 R13 R12, constructed out of universal R and reps \[Rho]1, \[Rho]2";
 RmatYBE::usage="RmatYBE[R,\[Rho]] gives the R-matrix satisfying YBE R12 R23 R12=R23 R12 R23, constructed out of universal R and rep \[Rho]";
 muln::usage="muln[mul,n] gives the n-fold multiplication map.";
 \[CapitalDelta]n::usage="\[CapitalDelta]n[\[CapitalDelta],n] gives the n-fold comultiplication map.";
 DTwistHA::usage="DTwistHA[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,J,Jb] verifies the DTwist and returns the twisted Hopf Algebra.";
 
  Begin[ "Private`"]
  Needs["TensorTools`"];
  id[x_]:=x;
  SA=SparseArray;
  LoadHA[HAName_]:=Module[{HAfilename,Repfilename},
  HAfilename=FileNameJoin[{NotebookDirectory[],"HALibrary",StringJoin[HAName,"Data.mx"]}];
  Repfilename=FileNameJoin[{NotebookDirectory[],"HALibrary",StringJoin[HAName,"Reps.mx"]}];
  Get[Repfilename];
 Import[HAfilename]
  ];
  muln[mul_,n_]:=Module[{dA,idt},
  If[n==1,Return[mul]];
  dA=Dimensions[mul][[1]];
  idt=Id[dA, SparseArray];
  mul.KP[muln[mul,n-1],idt]
  ];
  \[CapitalDelta]n[\[CapitalDelta]_,n_]:=Module[{dA,idt},
  If[n==1,Return[\[CapitalDelta]]];
  dA=Dimensions[\[CapitalDelta]][[2]];
  idt=Id[dA, SparseArray];
  KP[\[CapitalDelta]n[\[CapitalDelta],n-1],idt].\[CapitalDelta]
  ];
  GroupHA[MulTable_,InvTable_]:=Module[{n,mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,Rb},
  n=Length[MulTable];
  mul=ArrayReshape[SparseArray[Flatten@Table[{MulTable[[i,j]],i,j}->1,{i,1,n},{j,1,n}],{n,n,n}],{n,n^2}];
  \[CapitalDelta]=ArrayReshape[SparseArray[Flatten@Table[{i,i,i}->1,{i,1,n}],{n,n,n}],{n^2,n}];
  \[Eta]=SparseArray[{{1,1}->1},{n,1}];
  \[Epsilon]=ConstantArray[1,{1,n},SparseArray];
  S=SparseArray[Flatten@Table[{InvTable[[i]],i}->1,{i,1,n}],{n,n}];
  R=Rb=SparseArray[{1,1}->1,{n,n}];
  {mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,Rb}
  ];
  
  RepToMap[\[Rho]_, dA_]:=Module[{e},
  e[k_]:=SA@UnitVector[dA,k];
  (*d=Dimensions[\[Rho]@e[1]][[1]];*)
  Transpose[SA@Table[\[Rho][e[k]],{k,1,dA}],{3,1,2}]
  ];
  RepToCorep[\[Rho]_,R_]:=SparseArray[\[Rho].Transpose[R]];
  \[Rho]Reg[mul_]:=Module[{dA},
  dA=Dimensions[mul][[1]];
  Transpose[ArrayReshape[mul,{dA,dA,dA}],2<->3]
  ];
  \[Rho]Regr[mul_]:=Module[{dA},
  dA=Dimensions[mul][[1]];
  Transpose[ArrayReshape[mul,{dA,dA,dA}],1<->2]
  ];
  vRegl[\[CapitalDelta]_]:=Module[{dA},
  (*The regular left corep*)
  dA=Dimensions[\[CapitalDelta]][[2]];
  Transpose[ArrayReshape[\[CapitalDelta],{dA,dA,dA}],1<->3]
  ];
  vRegr[\[CapitalDelta]_]:=Module[{dA},
  (*The regular right corep*)
  dA=Dimensions[\[CapitalDelta]][[2]];
  Transpose[ArrayReshape[\[CapitalDelta],{dA,dA,dA}],2<->3]
  ];
  CheckRep[\[Rho]_,mul_,\[Eta]_]:=Module[{dA=Length[\[Eta]],\[Rho]T,d},
  d=Dimensions[\[Rho]][[1]];
   \[Rho]T=Transpose[\[Rho],2<->3];
   Transpose[ArrayReshape[\[Rho]T.\[Rho]T,{d,dA^2,d}],2<->3]==\[Rho].mul && ArrayReshape[\[Rho].\[Eta],{d,d}]==Id[d]
  ];
  CheckCorep[v_,\[CapitalDelta]_,\[Epsilon]_]:=Module[{dA,vT,d,vAdd},
  {d,d,dA}=Dimensions[v];
   vT=Transpose[v,2<->3];
   vAdd=Transpose[vT,1<->2];
   Transpose[ArrayReshape[vT.vT,{d,dA^2,d}],1<->2]==\[CapitalDelta].vAdd&& ArrayReshape[\[Epsilon].vAdd,{d,d}]==Id[d]
  ];
  HAData[dA_, Mul_, CoMul_, \[Epsilon]_, S_]:=Module[{MulTable, CoMulTable, \[Epsilon]Table, STable,e},
  e[k_]:=SA@UnitVector[dA,k];
  STable=SA@Table[S[e[i]],{i,1,dA}]//Simplify;
  MulTable=SA@Transpose[Table[Mul[e[i],e[j]],{i,1,dA},{j,1,dA}]//Simplify,2<->3];
  CoMulTable=SA@Table[CoMul[e[i]],{i,1,dA}]//Simplify;
  \[Epsilon]Table=SA@Table[\[Epsilon][e[i]],{i,1,dA}]//Simplify;
  {MulTable, CoMulTable, \[Epsilon]Table, STable}
  ];
  DirectSum[\[Rho]1_,\[Rho]2_]:=Module[{dA},
  dA=Dimensions[\[Rho]1][[3]];
  Transpose[Table[ArrayFlatten[{{\[Rho]1[[All,All,k]],0},{0,\[Rho]2[[All,All,k]]}}],{k,1,dA}],{3,1,2}]
  ];
  DirectSums[\[Rho]s__]:=If[Length[{\[Rho]s}]==2,DirectSum[\[Rho]s],DirectSum[{\[Rho]s}[[1]],DirectSums@@({\[Rho]s}[[2;;-1]])]];
 


VerifyAlgebra[mul_,\[Eta]_]:=Module[{dA=Length@\[Eta],idt},
idt=Id[dA, SparseArray];
(mul.KP[mul,idt]==mul.KP[idt,mul])&&(mul.KP[\[Eta],idt]==mul.KP[idt,\[Eta]])
];
VerifyCoalgebra[\[CapitalDelta]_,\[Epsilon]_]:=Module[{dA=Dimensions[\[Epsilon]][[2]],idt},
idt=Id[dA, SparseArray];
(KP[\[CapitalDelta],idt].\[CapitalDelta]==KP[idt,\[CapitalDelta]].\[CapitalDelta])&&(KP[\[Epsilon],idt].\[CapitalDelta]==KP[idt,\[Epsilon]].\[CapitalDelta])
];
Dual[mul_,\[Eta]_]:={Transpose@mul,Transpose@\[Eta]};
DualHA[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=Transpose/@{\[CapitalDelta],mul,\[Epsilon],\[Eta],S};
VerifyPreBAAxioms[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_]:=Module[{dA,idt,m2,\[Tau],Verify\[CapitalDelta]Multiplicative},
dA=Length@\[Eta];
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
m2=KP[mul,mul].KP[idt,\[Tau],idt];
Verify\[CapitalDelta]Multiplicative=m2.KP[\[CapitalDelta],\[CapitalDelta]]== \[CapitalDelta].mul;
VerifyAlgebra[mul,\[Eta]]&&VerifyCoalgebra[\[CapitalDelta],\[Epsilon]]&&Verify\[CapitalDelta]Multiplicative//FullSimplify
];
Verify\[Eta]\[Epsilon][mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_]:=(\[CapitalDelta].\[Eta]==KP[\[Eta],\[Eta]])&&(\[Epsilon].mul== KP[\[Epsilon],\[Epsilon]]);
VerifyBAAxioms[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_]:=VerifyPreBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]]&&Verify\[Eta]\[Epsilon][mul,\[CapitalDelta],\[Eta],\[Epsilon]];
VerifyWBAAxioms[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_]:=Module[{dA=Length@\[Eta],idt,VerifyWeakunit,VerifyWeakcounit,mm,\[CapitalDelta]\[CapitalDelta]},
\[Tau]=\[CapitalPi][dA];
idt=Id[dA, SparseArray];
\[CapitalDelta]\[CapitalDelta]=KP[\[CapitalDelta],idt].\[CapitalDelta];
mm=mul.KP[mul,idt];
VerifyWeakunit=\[CapitalDelta]\[CapitalDelta].\[Eta]==KP[idt,mul.\[Tau],idt].KP[\[CapitalDelta].\[Eta],\[CapitalDelta].\[Eta]]==KP[idt,mul,idt].KP[\[CapitalDelta].\[Eta],\[CapitalDelta].\[Eta]] //FullSimplify;
VerifyWeakcounit=\[Epsilon].mm== KP[\[Epsilon].mul,\[Epsilon].mul].KP[idt,\[CapitalDelta],idt]==  KP[\[Epsilon].mul,\[Epsilon].mul].KP[idt,\[Tau].\[CapitalDelta],idt]//FullSimplify;
VerifyPreBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]]&&VerifyWeakunit&&VerifyWeakcounit
];


VerifyWeakAntipode[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=Module[{dA=Length@\[Eta],\[Epsilon]t,\[Epsilon]s,idt,\[Tau],VerifySAntiAlgebra,VerifySAntiCoalgebra,VerifyWeakSAxioms},
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
\[Epsilon]t=KP[\[Epsilon].mul.\[Tau],idt].KP[idt,\[CapitalDelta].\[Eta]];
\[Epsilon]s=KP[idt,\[Epsilon].mul.\[Tau]].KP[\[CapitalDelta].\[Eta],idt];
VerifySAntiAlgebra=(S.mul==mul.KP[S,S].\[Tau])&&(S.\[Eta]==\[Eta]);
VerifySAntiCoalgebra=(\[CapitalDelta].S==\[Tau].KP[S,S].\[CapitalDelta])&&(\[Epsilon].S==\[Epsilon]);
VerifyWeakSAxioms=(mul.KP[S,idt].\[CapitalDelta]==\[Epsilon]s)&&(mul.KP[idt,S].\[CapitalDelta]==\[Epsilon]t);
VerifyWeakSAxioms&&VerifySAntiAlgebra&&VerifySAntiCoalgebra
];
VerifyWHAAxioms[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=VerifyWBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]]&&VerifyWeakAntipode[mul,\[CapitalDelta],\[Eta],\[Epsilon],S];
VerifyHAAxioms[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=VerifyBAAxioms[mul,\[CapitalDelta],\[Eta],\[Epsilon]]&&VerifyWeakAntipode[mul,\[CapitalDelta],\[Eta],\[Epsilon],S];
VerifyHA[HA_]:=VerifyHAAxioms@@HA[[1;;5]];
VerifyQHA[HA_]:=Module[{mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,Rb},
{mul,\[CapitalDelta],\[Eta],\[Epsilon],S,R,Rb}=HA;
(VerifyHAAxioms@@HA[[1;;5]])&&VerifyQuasitriangularity[mul,\[CapitalDelta],\[Eta],\[Epsilon],R,Rb]
];


VerifyQuasitriangularity[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,R_,Rb_]:=Module[{dA=Length@\[Eta],Rm,Rbm,VerifyRInv,VerifyR\[CapitalDelta]R,VerifyRCE,VerifyRCEL,VerifyRCER,idt,m2,\[Tau]},
m2:=KP[mul,mul].KP[idt,\[Tau],idt];
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
Rm=ArrayReshape[R,{dA^2,1}];
Rbm=ArrayReshape[Rb,{dA^2,1}];
VerifyRInv=FullSimplify[(m2.KP[Rbm,Rm]==\[CapitalDelta].\[Eta]) &&(m2.KP[Rm,Rbm]==\[Tau].\[CapitalDelta].\[Eta])];
VerifyR\[CapitalDelta]R= FullSimplify[m2.KP[\[Tau].\[CapitalDelta],Rm]== m2.KP[Rm,\[CapitalDelta]]];
VerifyRCEL=FullSimplify@(KP[\[CapitalDelta],idt].Rm== KP[idt,KP[idt,mul].KP[\[Tau],idt]].KP[Rm,Rm]);
VerifyRCER=FullSimplify@(KP[idt,\[Tau].\[CapitalDelta]].Rm== KP[KP[mul,idt].KP[idt,\[Tau]],idt].KP[Rm,Rm]);
VerifyRCE=VerifyRCEL&&VerifyRCER;
Print[{VerifyRInv,VerifyR\[CapitalDelta]R,VerifyRCEL,VerifyRCER}];
VerifyRInv&&VerifyR\[CapitalDelta]R&&VerifyRCE
];

(*VerifyDTwist[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,J_,Jb_]:=Module[{dA,idt,m2,\[CapitalDelta]mJ},
dA=Length@\[Eta];
m2:=KP[mul,mul].KP[idt,\[Tau],idt];
idt=Id[dA, SparseArray];
\[CapitalDelta]mJ=m2.KP[\[CapitalDelta],J];
(KP[\[Epsilon],idt].J\[Equal]KP[idt,\[Epsilon]].J\[Equal]\[Eta])&&(KP[\[CapitalDelta]mJ,idt].J\[Equal]KP[idt,\[CapitalDelta]mJ].J)
];*)

DTwistHA[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_,R_,J_,Jb_]:=Module[{dA,idt,m2,\[CapitalDelta]mJ,\[CapitalDelta]J,SJ,RJ,RJm,RJb,RJbm,w,winv,VerifyDTwist,Rm,J21b},
dA=Length@\[Eta];
\[Tau]=\[CapitalPi][dA];
m2=KP[mul,mul].KP[idt,\[Tau],idt];
idt=Id[dA, SparseArray];
\[CapitalDelta]mJ=m2.KP[\[CapitalDelta],J];
VerifyDTwist=(m2.KP[J,Jb]==KP[\[Eta],\[Eta]])&&(KP[\[Epsilon],idt].J==KP[idt,\[Epsilon]].J==\[Eta])&&(KP[\[CapitalDelta]mJ,idt].J==KP[idt,\[CapitalDelta]mJ].J);
Print@FullSimplify@VerifyDTwist;
\[CapitalDelta]J=m2.KP[Jb,\[CapitalDelta]mJ];
Rm=ArrayReshape[R,{dA^2,1}];
J21b=ArrayReshape[Transpose@ArrayReshape[Jb,{dA,dA}],{dA^2,1}];
RJm=m2.KP[J21b,m2.KP[Rm,J]];
RJ=ArrayReshape[RJm,{dA,dA}];
w=mul.KP[S,idt].J;
winv=mul.KP[idt,S].J21b;
SJ=mul.KP[winv,mul.KP[S,w]];
RJbm=KP[SJ,idt].RJm;
RJb=ArrayReshape[RJbm,{dA,dA}];
{mul,\[CapitalDelta]J,\[Eta],\[Epsilon],SJ,RJ,RJb}
];

AF2FA[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=AF2FA[mul,\[CapitalDelta],\[Eta],\[Epsilon],S]=Module[{dA,idt,mm,\[Tau],S\[CapitalDelta]\[CapitalDelta]cop,apart,fpart},
dA=Length@\[Eta];
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
mm=SA@mul.KP[mul,idt];

S\[CapitalDelta]\[CapitalDelta]cop=SA@KP[S,idt,idt].KP[\[Tau].\[CapitalDelta],idt].\[Tau].\[CapitalDelta];
fpart=ArrayReshape[Transpose[ArrayReshape[mm,{dA,dA,dA,dA}],2<->3],{dA,dA,dA^2}];
apart=ArrayReshape[Transpose[ArrayReshape[S\[CapitalDelta]\[CapitalDelta]cop,{dA,dA,dA,dA}],2<->3],{dA^2,dA,dA}];
SA@ArrayReshape[Transpose[fpart.apart,{4,1,2,3}],{dA^2,dA^2}]
];

DHRep[\[Rho]H_,\[Rho]Hd_,af2fa_]:=Module[{\[Rho]Ht,\[Rho]Hdt,LHS,RHS},
{\[Rho]Ht,\[Rho]Hdt}=Transpose[#,2<->3]&/@(SA/@{\[Rho]H,\[Rho]Hd});
LHS=Flatten[\[Rho]Ht.\[Rho]Hdt,{{1},{4},{2,3}}];
RHS=Flatten[\[Rho]Hdt.\[Rho]Ht,{{1},{4},{2,3}}];
(*Print[LHS,RHS.af2fa];*)
(*Print[LHS\[Equal]RHS.af2fa];*)
SA@RHS
];

DHaction\[Rho]H[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_,\[Rho]H_]:=Module[{d\[Rho],dA,af2fa,af2faT,muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd, \[Rho]Hd,\[Rho]H2},
dA=Length@\[Eta];
{muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd}=DualHA[mul,\[CapitalDelta], \[Eta], \[Epsilon], S];
d\[Rho]=Dimensions[\[Rho]H][[1]];
\[Rho]Hd=KP[\[Rho]Reg[muld],ArrayReshape[Id[d\[Rho],SparseArray],{d\[Rho],d\[Rho],1}]];
(*Print[\[Rho]Hd//Dimensions];*)
af2fa=AF2FA[mul,\[CapitalDelta],\[Eta],\[Epsilon],S];
af2faT=Transpose[ArrayReshape[af2fa,{dA,dA,dA,dA}],{2,1,4,3}];
\[Rho]H2=Flatten[\[Rho]H.af2faT,{{3,1},{4,2},{5}}];
(*Print@CheckRep[\[Rho]H2,mul,\[Eta]];
Print@CheckRep[\[Rho]Hd,muld,\[Eta]d];
Print[\[Rho]H2//Dimensions];*)
SA@DHRep[\[Rho]H2,\[Rho]Hd,af2fa]
];

DHaction\[Rho]HSub[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_,\[Rho]H_,V_]:=Module[{d\[Rho],dA,af2fa,af2faT,muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd, \[Rho]Hd,\[Rho]H2,Psub},
dA=Length@\[Eta];
{muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd}=DualHA[mul,\[CapitalDelta], \[Eta], \[Epsilon], S];
d\[Rho]=Dimensions[\[Rho]H][[1]];
Psub=SparseArray[Orthogonalize[V]];
\[Rho]Hd=KP[Transpose[Psub.Transpose[\[Rho]Reg[muld],2<->3].ConjugateTranspose[Psub],2<->3],ArrayReshape[Id[d\[Rho],SparseArray],{d\[Rho],d\[Rho],1}]];

af2fa=AF2FA[mul,\[CapitalDelta],\[Eta],\[Epsilon],S];
af2faT=Transpose[Psub.ArrayReshape[af2fa,{dA,dA,dA,dA}].ConjugateTranspose[Psub],{2,1,4,3}];
\[Rho]H2=Flatten[SA[\[Rho]H.af2faT],{{3,1},{4,2},{5}}];
(*Print[\[Rho]H2," ",\[Rho]Hd];*)
SA@DHRep[\[Rho]H2,\[Rho]Hd,af2fa]
];

DGRep[CC_,CG_,\[Rho]CG_]:=Module[{mul,\[CapitalDelta],\[Eta],\[Epsilon],S,V,dA},
{mul,\[CapitalDelta],\[Eta],\[Epsilon],S}=CG;
dA=Length@\[Eta];
V=Table[UnitVector[dA,k],{k,CC}];
DHaction\[Rho]HSub[mul,\[CapitalDelta],\[Eta],\[Epsilon],S,\[Rho]CG,V]
];

DHaction\[Rho]Hd[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_,\[Rho]Hd_]:=Module[{muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd,dA,fa2af},
dA=Length@\[Eta];
{muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd}=DualHA[mul,\[CapitalDelta], \[Eta], \[Epsilon], S];
fa2af=AF2FA[muld,\[CapitalDelta]d,\[Eta]d,\[Epsilon]d,Sd];
SA[DHaction\[Rho]H[muld,\[CapitalDelta]d,\[Eta]d,\[Epsilon]d,Sd,\[Rho]Hd].fa2af]
];

DrinfeldDouble[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_]:=Module[{dA,idt,\[Tau],muld,\[CapitalDelta]d,\[Eta]d, \[Epsilon]d, Sd,mulD,\[CapitalDelta]D, \[Eta]D, \[Epsilon]D, SD,af2fa,RD,RDb,RDm,RDbm},
dA=Length@\[Eta];
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
af2fa=AF2FA[mul,\[CapitalDelta],\[Eta],\[Epsilon],S];
{muld,\[CapitalDelta]d, \[Eta]d, \[Epsilon]d, Sd}=DualHA[mul,\[CapitalDelta], \[Eta], \[Epsilon], S];
mulD=KP[muld,mul].KP[idt,af2fa,idt];
\[Eta]D=KP[\[Eta]d,\[Eta]];

\[CapitalDelta]D=KP[idt,\[Tau],idt].KP[\[Tau].\[CapitalDelta]d,\[CapitalDelta]];
\[Epsilon]D=KP[\[Epsilon]d,\[Epsilon]];
SD=af2fa.KP[S,Sd].\[Tau];

RDm=KP[\[Eta]d,ArrayReshape[Id[dA,SparseArray],{dA^2,1}],\[Eta]];
RDbm=KP[SD,Id[dA^2,SparseArray]].RDm;
RD=ArrayReshape[RDm,{dA^2,dA^2}];
RDb=ArrayReshape[RDbm,{dA^2,dA^2}];
SA/@{mulD,\[CapitalDelta]D, \[Eta]D, \[Epsilon]D, SD,RD,RDb}
];

UniversalRX[mul_,\[CapitalDelta]_,\[Eta]_,\[Epsilon]_,S_,R_,Rb_]:=Module[{dA,idt,\[Tau],UL,VL,UR,VR,R21,R21m,Rb21,Rb21m},
dA=Length@\[Eta];
idt=Id[dA, SparseArray];
\[Tau]=\[CapitalPi][dA];
R21=Transpose@R;
R21m=ArrayReshape[R21,{dA^2,1}];
Rb21=Transpose@Rb;
Rb21m=ArrayReshape[Rb21,{dA^2,1}];
UL=KP[idt,mul].KP[R21m,idt]; 
UR=KP[idt,mul.\[Tau]].KP[Rb21m,idt];
VL=KP[S,idt].\[Tau].\[CapitalDelta];
VR=\[CapitalDelta];
{UL,UR,VL,VR}
]

VerifyUniversalRX[mul_,R_,UL_,UR_,VL_,VR_]:=Module[{dA,\[Tau],\[Tau]2,idt,m2,TF,RXMul,Rm,CheckRX,CheckRUU,VerifyRUU},
dA=Dimensions[mul][[1]];
idt=Id[dA, SparseArray];
m2=KP[mul,mul].KP[idt,\[Tau],idt];
\[Tau]=\[CapitalPi][dA];
\[Tau]2=\[CapitalPi][dA^2];
Rm=ArrayReshape[R,{dA^2,1}];
TF[U_]:=Transpose[ArrayReshape[U,{dA,dA,dA}],1<->2];
RXMul[U_,V_]:=ArrayReshape[Transpose[ArrayReshape[TF[U].TF[V],{dA,dA^2,dA}],1<->2],{dA^2,dA^2}];
(*CheckUniversalRX=*)
CheckRX={m2.KP[Rm,\[Tau].RXMul[UL,VR]]==RXMul[VR,UL],
m2.KP[Rm,\[Tau].RXMul[UR,VL]]==RXMul[VL,UR],
m2.\[Tau]2.KP[Rm,RXMul[VL,UL]]==\[Tau].RXMul[UL,VL],
m2.\[Tau]2.KP[Rm,RXMul[VR,UR]]==\[Tau].RXMul[UR,VR],
RXMul[UL,UR]==\[Tau].RXMul[UR,UL]&&RXMul[VL,VR]==\[Tau].RXMul[VR,VL]};

VerifyRUU[U_]:=m2.KP[Rm,\[Tau].RXMul[U,U]]==m2.\[Tau]2.KP[Rm,RXMul[U,U]];
(*CheckRUU=And@@VerifyRUU/@{UL,UR,VL,VR};*)
CheckRUU=VerifyRUU/@{UL,UR,VL,VR};
(*CheckRX&&CheckRUU*)
{CheckRX,CheckRUU}
];
DrinfeldElement[mul_,S_,R_]:=Module[{dA,idt,R21m,u,uinv},
dA=Dimensions[mul][[1]];
idt=Id[dA, SparseArray];
R21m=ArrayReshape[Transpose@R,{dA^2,1}];
u=mul.KP[S,idt].R21m;
uinv=mul.KP[idt,S.S].R21m;
{u,uinv}
];
HaarIntegral[\[CapitalDelta]_]:=1/Dimensions[\[CapitalDelta]][[2]] Tr[vRegl[\[CapitalDelta]],Plus,2];
RmatQYBE[R_,\[Rho]1_,\[Rho]2_]:=Module[{d1,d2,dA},
{d1,d1,dA}=Dimensions[\[Rho]1];
{d2,d2,dA}=Dimensions[\[Rho]2];
ArrayReshape[Transpose[\[Rho]1.R.Transpose[\[Rho]2,{2,3,1}],2<->3],{d1 d2,d1 d2}]
];
RmatYBE[R_,\[Rho]_]:=Module[{d,dA},
{d,d,dA}=Dimensions[\[Rho]];
RmatQYBE[R,\[Rho],\[Rho]].\[CapitalPi][d]
];


  End[]

  EndPackage[]

