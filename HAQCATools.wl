(* ::Package:: *)

BeginPackage[ "HAQCATools`"]

  QCATensors::usage = 
	"QCATensors[mul,\[CapitalDelta], \[Eta], \[Epsilon], S, \[Rho], v]: input HA data {mul,\[CapitalDelta], \[Eta], \[Epsilon], S}, a representation \[Rho] and a corepresentation v, 
 output all the basic tensors {Ugate,\[Rho]LTensor,\[Rho]RTensor,vLTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor} necessary for QCA computations.";
  VerifyPentagon::usage=" VerifyPentagon[Ugate,\[Rho]LTensor,\[Rho]RTensor,vLTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor] verifies the pentagon equation.";
  VerifyBoundary\[Eta]\[Epsilon]::usage="VerifyBoundary\[Eta]\[Epsilon][Ugate,\[Rho]LTensor,\[Rho]RTensor,vLTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor] returns {VerifyBA\[Eta]\[Epsilon],VerifypreBA\[Eta]\[Epsilon]}.";
  SWAP::usage="SWAP[d\[Rho],dv] returns the d\[Rho] dv-dimensional swap gate.";
  Corr::usage="Corr[O\[Rho],Ov,\[CapitalPsi],nx,t,\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor] computes the equal time  correlation function O\[Rho](t)Ov(t) for an initial MPS \[CapitalPsi].";
  TMRenyiSmall\[Alpha]::usage="TMRenyiSmall\[Alpha][\[CapitalPsi], \[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor,\[Alpha]] computes all the important TMs {TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary}
  relevant for the computation of Renyi entanglement entropy with finite index \[Alpha]. ";
  RenyiSemifinite::usage="RenyiSemifinite[\[Alpha],t,\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor] computes the bipartite Renyi \[Alpha]-entanglement entropy 
   at time t for a semi-infinite subsystem";
  RenyiFiniteSystemSmall\[Alpha]::usage="RenyiFiniteSystemSmall\[Alpha][\[Alpha],L,t,\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor] computes the Renyi \[Alpha]-entropy for finite subsystem with length L
  (an integer) at time t, for initial state \[CapitalPsi].";
  TMRenyiSmallSystemLarge\[Alpha]::usage="TMRenyiSmallSystemLarge\[Alpha][\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor].";
  RenyiSmallSystemLarge\[Alpha]::usage="RenyiSmallSystemLarge\[Alpha][\[Alpha],L,t,TMs] computes Renyi \[Alpha]-entropy for small system of size L at time t.";
  RenyiEE::usage="RenyiEE[\[Rho],\[Alpha]] computes the Renyi \[Alpha]-entropy of the density matrix \[Rho].";
  RFCrank::usage="RFCrank[U,m] computes the RFC entanglement rank of the m*m square TN generated by the gate U."
  DihedralGate::usage="DihedralGate[n,k] gives the QCA gate constructed from the dihedral group Dn."
  PBCEvol::usage="PBCEvol[U,L] gives the PBC evolution operator Ue.Uo, where L is the system site and U is a 2-site gate."
  Mpm::usage="Mpm[U] computes {Mp,Mm}, the M+ (for the v direction) and M- (for the \[Rho] direction) matrices for the computation of the ST corr for a dual unitary gate U."
  Begin[ "Private`"]
  Needs["TensorTools`"];
  SA=SparseArray;
  SWAP[d\[Rho]_,dv_]:=ArrayReshape[Transpose[ArrayReshape[Id[d\[Rho] dv],{d\[Rho],dv,d\[Rho],dv}]],{d\[Rho] dv,d\[Rho] dv}];
opj[O_,j_,n_]:=Module[{d},
(*Here O is assumed to be a 2-site operator!*)
d=Sqrt[Length[O]];
KP[SparseArray[Id[d^(j-1)]],SparseArray[O],SparseArray[Id[d^(n-j-1)]]]
];
PBCEvol[U_,L_]:=Module[{Shift,Ue,Uo,j,d},
d=Sqrt[Length[U]];
Shift=SparseArray[Id[d^(2L)]];
For[j=1,j<=2L-1,j++,Shift=Shift.opj[\[CapitalPi][d],j,2L]];
Uo=SparseArray[Id[d^(2L)]];
Ue=ConjugateTranspose[Shift].KP[U,Id[d^(2L-2)]].Shift;
(*Print["OK"];*)
For[j=1,j<=L,j++,Uo=Uo.opj[U,2j-1,2L]];
For[j=1,j<=L-1,j++,Ue=Ue.opj[U,2j,2L]];
{Ue.Uo,Ue,Uo}
];

RFCrank[U_,m_]:=Module[{Um,d=Sqrt[Dimensions[U][[1]]],UmR},
Um=SquareTN[U,m];
UmR=Flatten[ArrayReshape[Um,{d^m,d^m,d^m,d^m}],{{1,4},{2,3}}];
MatrixRank@UmR
];
DihedralGate[n_,k_]:=Module[{\[Omega],s,r,P1,P2},
\[Omega]=Exp[(2 \[Pi] I)/n];
       s={{0,1},{1,0}};
       r=DiagonalMatrix[{\[Omega]^k,\[Omega]^-k}];
P1=DiagonalMatrix[{1,0}];
P2=DiagonalMatrix[{0,1}];
\[CapitalPi][2].(KP[s,P1]+KP[r,P2])
];
Mpm[U_]:=Module[{d=Sqrt[Dimensions[U][[1]]],UT,UTD,Mp,Mm,IT},
UT=ArrayReshape[U,{d,d,d,d}];
UTD=KP[UT,Conjugate@UT];
IT=Flatten@Id[d];
Mm=IT.UTD.IT/d;
Mp=IT.Transpose[UTD,{2,1,4,3}].IT/d;
{Mp,Mm}
];
  
  QCATensors[mul_,\[CapitalDelta]_, \[Eta]_, \[Epsilon]_, S_, \[Rho]_, v_] :=
    Module[ {dA,d\[Rho],dv,\[CapitalDelta]Tensorcop,mulTensorop,Ugate,\[Rho]LTensor,\[Rho]RTensor,vLTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor,vT,mulTensor,\[CapitalDelta]Tensor},
    d\[Rho]=Dimensions[\[Rho]][[1]];
    dv=Dimensions[v][[1]];
    dA=Length[\[Eta]];
    \[CapitalDelta]Tensor=ArrayReshape[\[CapitalDelta],{dA,dA,dA}];
    mulTensor=ArrayReshape[mul,{dA,dA,dA}];
    \[CapitalDelta]Tensorcop=Transpose[\[CapitalDelta]Tensor,1<->2];
    mulTensorop=Transpose[mulTensor,2<->3];
    vT=Transpose[v,{2,3,1}];
    
    Ugate=ArrayReshape[Transpose[\[Rho].vT,{2,3,1,4}],{d\[Rho] dv, d\[Rho] dv}];
    \[Rho]LTensor=Transpose[\[Rho].\[CapitalDelta]Tensor,{1,4,3,2}];
    \[Rho]RTensor=Transpose[\[Rho].\[CapitalDelta]Tensorcop,{1,4,2,3}];
    vRTensor=Transpose[mulTensor.vT,{2,3,1,4}];
    vLTensor=Transpose[mulTensorop.vT,{3,2,1,4}];
    
    \[Eta]Tensor=Flatten@\[Eta];
    \[Epsilon]Tensor=Flatten@\[Epsilon];
    {Ugate,\[Rho]LTensor,\[Rho]RTensor,vLTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor}
    ];
  VerifyPentagon[Ugate_,\[Rho]LTensor_,\[Rho]RTensor_,vLTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{dA,d\[Rho],dv,PentagonR,PentagonL},
    d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
    dv=Dimensions[vRTensor][[1]];
    dA=Length[\[Eta]Tensor];
    
  PentagonR=(SWAP[d\[Rho],dv].ContractH[\[Rho]RTensor,vRTensor]== ContractH[vRTensor,\[Rho]RTensor].Ugate);
  PentagonL=(Ugate.ContractH[\[Rho]LTensor,vLTensor]== ContractH[vLTensor,\[Rho]LTensor].SWAP[d\[Rho],dv]);
  (*\[Epsilon]1=KP[\[Epsilon]Tensor,Id[]]*)
  {PentagonL,PentagonR}
  ];
  VerifyBoundary\[Eta]\[Epsilon][Ugate_,\[Rho]LTensor_,\[Rho]RTensor_,vLTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{dA,d\[Rho],dv,\[Epsilon]1,\[Eta]1,VerifyBA\[Eta]\[Epsilon],VerifypreBA\[Eta]\[Epsilon],\[Rho]vContracted,\[Rho]3v3,\[Rho]3v4,\[Rho]4v3},
    d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
    dv=Dimensions[vRTensor][[1]];
    dA=Length[\[Eta]Tensor];
    VerifyBA\[Eta]\[Epsilon]=(V2H[\[Rho]RTensor].\[Eta]Tensor==ArrayReshape[KP[\[Eta]Tensor,Id[d\[Rho]]],{dA,d\[Rho],d\[Rho]}])&&(\[Epsilon]Tensor.V2H[vRTensor]==ArrayReshape[KP[Id[dv],\[Epsilon]Tensor],{dv,dv,dA}]);
    \[Rho]vContracted=ContractH[\[Rho]RTensor,vRTensor];
    
   \[Rho]4v3=(SWAP[d\[Rho],dv].ContractH41[\[Rho]vContracted,\[Eta]Tensor]== ContractH41[ContractH[vRTensor,SwapTensor[d\[Rho],dA]],\[Eta]Tensor].Ugate);
   (*Print[1];*)
   \[Rho]3v4=(SWAP[d\[Rho],dv].ContractH14[\[Epsilon]Tensor,\[Rho]vContracted]== ContractH14[\[Epsilon]Tensor,ContractH[SwapTensor[d\[Rho],dA],\[Rho]RTensor]].Ugate);
  (* Print[2];*)
   \[Rho]3v3=SWAP[d\[Rho],dv].Ugate==(\[Epsilon]Tensor.V2H[\[Rho]vContracted].\[Eta]Tensor);
   (*Print[\[Rho]3v3,\[Rho]3v4,\[Rho]4v3]*);
   VerifypreBA\[Eta]\[Epsilon]=\[Rho]3v3&&\[Rho]3v4&&\[Rho]4v3;
   {VerifyBA\[Eta]\[Epsilon],VerifypreBA\[Eta]\[Epsilon]}
  ];
BoundaryVector[\[Chi]_,\[Eta]_]:=Module[{dA,D\[Psi]},
dA=Length[\[Eta]];
D\[Psi]=Sqrt[Length[\[Chi]]];
Flatten@Transpose[ArrayReshape[KP[\[Chi],Flatten@KP[Conjugate[\[Eta]],\[Eta]]],{D\[Psi],D\[Psi],dA,dA}],{1,4,2,3}]
];

Corr[Ov_?MatrixQ,O\[Rho]_?MatrixQ,\[CapitalPsi]_,nx_,t_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v,nt,dA,d\[Rho],dv,T\[Rho],Tv,TO\[Rho],TOv,\[Eta]\[Eta],\[Epsilon]\[Epsilon],Ttemp},
{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v}=\[CapitalPsi];
(*D\[Psi]= Dimensions[\[Psi]\[Rho]][[2]];*)

nt=2t+nx;
d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
dv=Dimensions[vRTensor][[1]];
dA=Length[\[Eta]Tensor];

\[Eta]\[Eta]=BoundaryVector[\[Chi]R,\[Eta]Tensor];
\[Epsilon]\[Epsilon]=BoundaryVector[\[Chi]L,\[Epsilon]Tensor];

If[t==0,
T\[Rho]=TM\[Psi][\[Psi]\[Rho]];
Tv=TM\[Psi][\[Psi]v];
TO\[Rho]=TMO\[Psi][O\[Rho],\[Psi]\[Rho]];
TOv=TMO\[Psi][Ov,\[Psi]v];
Return[\[Chi]L.TOv.(T\[Rho].Tv)^nx.TO\[Rho].\[Chi]R]];

T\[Rho]=TM\[Psi][ContractV[\[Rho]RTensor,\[Psi]\[Rho]]];
Tv=TM\[Psi][ContractV[vRTensor,\[Psi]v]];
TO\[Rho]=TMO\[Psi][O\[Rho],ContractV[\[Rho]RTensor,\[Psi]\[Rho]]];
TOv=TMO\[Psi][Ov,ContractV[vRTensor,\[Psi]v]];

If[nt>= 2nx+1,
Ttemp=MatrixPow[T\[Rho].Tv,nx];
Return[\[Epsilon]\[Epsilon].Ttemp.TO\[Rho].MatrixPow[Tv.T\[Rho],nt-2nx-1].TOv.Ttemp.\[Eta]\[Eta]];
];

If[nt<= 2nx,
Ttemp=MatrixPow[T\[Rho].Tv,nt-nx-1];
Return[\[Epsilon]\[Epsilon].Ttemp.T\[Rho].TOv.MatrixPow[T\[Rho].Tv,2nx-nt].TO\[Rho].Tv.Ttemp.\[Eta]\[Eta]];
];

];

TMRenyiOneLayer[\[CapitalPsi]_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v,dA,d\[Rho],dv,D\[Psi],T\[Rho],Tv,T\[Rho]b,Tvb,\[Eta]\[Eta],\[Epsilon]\[Epsilon]},
{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v}=\[CapitalPsi];
D\[Psi]= Dimensions[\[Psi]\[Rho]][[2]];

d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
dv=Dimensions[vRTensor][[1]];
dA=Length[\[Eta]Tensor];

\[Eta]\[Eta]=BoundaryVector[\[Chi]R,\[Eta]Tensor];
\[Epsilon]\[Epsilon]=BoundaryVector[\[Chi]L,\[Epsilon]Tensor];

T\[Rho]=TM\[Psi][ContractV[\[Rho]RTensor,\[Psi]\[Rho]]];
Tv=TM\[Psi][ContractV[vRTensor,\[Psi]v]];
T\[Rho]b=SwapX[T\[Rho],dA D\[Psi], dA D\[Psi]];
Tvb=SwapX[Tv,dA D\[Psi], dA D\[Psi]];

{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv,T\[Rho]b,Tvb}
];

TMOt[\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{dA,d\[Rho],dv,D\[Psi],T\[Rho],Tv,T\[Rho]b,Tvb,\[Eta]\[Eta],\[Epsilon]\[Epsilon]},
d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
dv=Dimensions[vRTensor][[1]];
dA=Length[\[Eta]Tensor];


T\[Rho]=TM\[Psi][ContractV[\[Rho]RTensor,\[Psi]\[Rho]]];
Tv=TM\[Psi][ContractV[vRTensor,\[Psi]v]];

{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv}
];


TMRenyiSmall\[Alpha][\[CapitalPsi]_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_,\[Alpha]_]:=TMRenyiSmall\[Alpha][\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor,\[Alpha]]=Module[{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v,dA,d\[Rho],dv,D\[Psi],T\[Rho],Tv,T\[Rho]b,Tvb,TT\[Rho],TTv,TT\[Rho]b,TTvb,\[Eta]\[Eta],\[Epsilon]\[Epsilon],RBoundary,LBoundary},
{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v}=\[CapitalPsi];
D\[Psi]= Dimensions[\[Psi]\[Rho]][[2]];

d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
dv=Dimensions[vRTensor][[1]];
dA=Length[\[Eta]Tensor];

{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv,T\[Rho]b,Tvb}=TMRenyiOneLayer[\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor];

RBoundary=Flatten@KPPow[\[Eta]\[Eta],\[Alpha]];
LBoundary=Flatten@KPPow[\[Epsilon]\[Epsilon],\[Alpha]];

(*{T\[Rho],Tv,T\[Rho]b,Tvb}*)
TT\[Rho]=KPPow[T\[Rho],\[Alpha]];
TTv=KPPow[Tv,\[Alpha]];
TT\[Rho]b=SwapX[KPPow[T\[Rho]b,\[Alpha]],dA D\[Psi],(dA D\[Psi])^(2\[Alpha]-1)];
TTvb=SwapX[KPPow[Tvb,\[Alpha]],dA D\[Psi],(dA D\[Psi])^(2\[Alpha]-1)];
{TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary}
];

RenyiSemifinite[\[Alpha]_,t_,\[CapitalPsi]_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary,nt},
{TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary}=TMRenyiSmall\[Alpha][\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor,\[Alpha]];
nt=2t;
(*LBoundary.MatrixPow[TT\[Rho]b.TTv,nt].RBoundary*)
If[t==0,Return[0.0]];
1/(1-\[Alpha]) Log[LBoundary.MatrixPow[TT\[Rho]b.TTv,nt].RBoundary]
];

RenyiFiniteSystemSmall\[Alpha][\[Alpha]_,L_,t_,\[CapitalPsi]_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=RenyiFiniteSystemSmall\[Alpha][\[Alpha],L,t,\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor]=Module[{TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary,nt,nx},
If[t==0,Return[0.0]];
{TT\[Rho],TTv,TT\[Rho]b,TTvb,LBoundary,RBoundary}=TMRenyiSmall\[Alpha][\[CapitalPsi],\[Rho]RTensor,vRTensor,\[Eta]Tensor,\[Epsilon]Tensor,\[Alpha]];
nx=L;
nt=2t+nx;

If[nt>= 2nx,
1/(1-\[Alpha]) Log[LBoundary.MatrixPow[TT\[Rho]b.TTv,nx].MatrixPow[TT\[Rho].TTv,2t-nx].MatrixPow[TT\[Rho].TTvb,nx].RBoundary],
1/(1-\[Alpha]) Log[LBoundary.MatrixPow[TT\[Rho]b.TTv,2t].MatrixPow[TT\[Rho]b.TTvb,nx-2t].MatrixPow[TT\[Rho].TTvb,2t].RBoundary]
]
];

RenyiEE[\[Rho]_,\[Alpha]_]:=1/(1-\[Alpha]) Log[Tr[MatrixPow[\[Rho],\[Alpha]]]];
(*TMIndexTranspose[n_]:={Table[i,{i,1,n}],Table[i,{i,n+1,2n}]}//Transpose//Flatten;*)
DMTensor2Mat[DMTensor_]:=Module[{n},
n=ArrayDepth[DMTensor]/2;
Flatten[DMTensor,{Range[1,2n-1,2],Range[2,2n,2]}]
];

TMRenyiSmallSystemLarge\[Alpha][\[CapitalPsi]_,\[Rho]RTensor_,vRTensor_,\[Eta]Tensor_,\[Epsilon]Tensor_]:=Module[{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v,dA,d\[Rho],dv,D\[Psi],T\[Rho],Tv,T\[Rho]b,Tvb,\[Eta]\[Eta],\[Epsilon]\[Epsilon],\[Rho]\[Psi]\[Rho],v\[Psi]v,K\[Rho],Kv},
{\[Chi]L,\[Chi]R,\[Psi]\[Rho],\[Psi]v}=\[CapitalPsi];
D\[Psi]= Dimensions[\[Psi]\[Rho]][[2]];

d\[Rho]=Dimensions[\[Rho]RTensor][[1]];
dv=Dimensions[vRTensor][[1]];
dA=Length[\[Eta]Tensor];

\[Eta]\[Eta]=BoundaryVector[\[Chi]R,\[Eta]Tensor];
\[Epsilon]\[Epsilon]=BoundaryVector[\[Chi]L,\[Epsilon]Tensor];

\[Rho]\[Psi]\[Rho]=ContractV[\[Rho]RTensor,\[Psi]\[Rho]];
v\[Psi]v=ContractV[vRTensor,\[Psi]v];
(*Print@\[Rho]\[Psi]\[Rho];*)
T\[Rho]=TM\[Psi][\[Rho]\[Psi]\[Rho]];
Tv=TM\[Psi][v\[Psi]v];

K\[Rho]=Transpose[ArrayReshape[KP[Conjugate[\[Rho]\[Psi]\[Rho]],\[Rho]\[Psi]\[Rho]],{d\[Rho],d\[Rho],dA^2,dA^2}],{2,3,1,4}];
Kv=Transpose[ArrayReshape[KP[Conjugate[v\[Psi]v],v\[Psi]v],{dv,dv,dA^2,dA^2}],{2,3,1,4}];

{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv,K\[Rho],Kv}
];
RenyiSmallSystemLarge\[Alpha][\[Alpha]_,L_,t_,TMs_]:=RenyiSmallSystemLarge\[Alpha][\[Alpha],L,t,TMs]=Module[{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv,K\[Rho],Kv,nx,nt,DMTensor,DM},

If[t==0,Return[0.0]];
{\[Eta]\[Eta],\[Epsilon]\[Epsilon],T\[Rho],Tv,K\[Rho],Kv}=TMs;
nx=L;
nt=2t+nx;

DMTensor=
If[nt>= 2nx,
\[Epsilon]\[Epsilon].TensorPow[K\[Rho].Tv,nx].MatrixPow[T\[Rho].Tv,2t-nx].TensorPow[T\[Rho].Kv,nx].\[Eta]\[Eta],
\[Epsilon]\[Epsilon].TensorPow[K\[Rho].Tv,2t].TensorPow[K\[Rho].Kv,nx-2t].TensorPow[T\[Rho].Kv,2t].\[Eta]\[Eta]
];
DM=DMTensor2Mat[DMTensor];
(*{DM,RenyiEE[DM,\[Alpha]]}*)
RenyiEE[DM,\[Alpha]]
];





  End[]

  EndPackage[]

