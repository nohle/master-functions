(* ::Package:: *)

(* 
Many of the functions related to graph generation, graph isomorphisms, 
pure YM Feynman rules, momentum conservation, relabeling, as well as
many of the structure definitions are thanks to Tristan Dennen. 
Much of this code builds off of his great foundation. Many thanks 
to Tristan Dennen, Scott Davies and Zvi Bern for code and/or the 
physical insights therein.
*)


(* 
Note: MonsterSolve is a propriety function that has been removed. 
Contact me for more information. It only affects some color functions. 
Many of the results are already saved though.
*)


(* ::Subtitle::Closed:: *)
(****************************************)
(*Public Initializations*)
(****************************************)


BeginPackage["MasterFunctions`"];
ClearAll["MasterFunctions`*"];
(* 
- Since we use Mathematica 8 often, I needed to demand that Mathematica use the 
  Combinatorica package instead of its native functions for this to run on mma10, 
  especially in the isomorphism code. Will updates make these faster??
*)
Needs["Combinatorica`"];
(* 
- For ContractionCycles[] in Tensor Reduction -> Coefficients from Tensor Contractions
- Mma 8 on OSX 10.10.5 with XCode 6 had some issues compiling in C
  before upgrading to XCode 7.
	- Error: "ld: library not found for -lgcc_s.10.4
			  nclang: error: linker command failed with exit code 1 (use -v to see invocation)"
  So if you have issues, comment out CompilationTarget -> "C" in ContractionCycles[]. 
  Then it will compile to "WVM".
*)
Needs["CCompilerDriver`"];


(* Needed to find saved results in subdirectories *)
MASTERFUNCTIONSDIRECTORY = $HomeDirectory<>"/Desktop/MasterFunctions/";
AppendTo[$Path, MASTERFUNCTIONSDIRECTORY];


(* We utilize parallel processing, so launch all kernels *)
(* I have 2 cores, but mma8 will launch 4 due to hyperthreading
	mma10 will only launch 2. Mma8 seems faster though. *)
LaunchKernels[];
(* For ParallelTable, which has a habit of severe misbehavior! Beware! *)
distributedContexts = {"MasterFunctions`Private`","MasterFunctions`"};


$MinPrecision = 100;
(*Used in Cut Numerics*)
rationalizePrecision = 50;
(*Used in RandComp*)
randSize = 100;


(* HARD CODED *)
(* Non-kinematic symbols that might appear in an ansatz or Feynman rules *)
(* Used in analytical cutting *)
SYMBOLSLIST = {Ds,Df,c[__],a[__]};
(* Used in tensor reduction to simplify things *)
ZEROREPS = {kk[i_,i_]:>0, ke[i_,i_]:>0};


(* Optional printing. Useful for debugging. Annoying otherwise. *)
PRINTQ = False;
PrintQ[expr_] := If[PRINTQ,Print[expr]];


{Ds,ke,le,kk,lk,ll,ee,EE,kE,lE,eE,s,t,d,l,p,\[Epsilon],\[Epsilon]2,\[Mu],\[Eta],\[CapitalGamma],Df,fColor,RG,BG,R,S,AngleBracket,bg,a,b,c,\[Beta],\[Alpha],tree,times,h0,F,COLOR};

(****************************************************)
(* MomConsN, Relab, Propagate, GetVars *)
(****************************************************)
MomConsN::usage = "MomConsN[# of legs, expression, in terms of s[i,j,...]'s = False]";
Relab::usage = "Relab[list of replacement rules, expression]";
Propagate::usage = "Propagate[expression with p[_][\[Mu][__]], \[Epsilon][_][\[Mu][__], \[Eta][\[Mu][__],\[Mu][__]], ...]";
GetVars::usage = "GetVars[expression, {type to exclude, eg. a[__] or ke[1,2]}]";
(****************************************************)
(* Graph Information *)
(****************************************************)
DisplayGraph::usage = "";
GraphInfo::usage = "";
MomRules::usage = "";
ChooseLoopMom::usage = "";
CutSolution::usage = "";
GetCutReps::usage = "";
(****************************************************)
(* Graph Isomorphisms *)
(****************************************************)
BuildGraph::usage = "";
EdgewiseAutomorphisms::usage = "";
EdgewiseIsomorphism::usage = "";
LabelwiseAutomorphisms::usage = "";
LabelwiseIsomorphism::usage = "";
PreIsoInvariant::usage = "";
isIso::usage = "";
FindIso::usage = "";
(****************************************************)
(* Ansatz Builder *)
(****************************************************)
BuildAnsatz::usage = "";
(****************************************************)
(* Linear Equation Solver *)
(****************************************************)
MonsterSolve2::usage = "";
(****************************************************)
(* General Numerics *)
(****************************************************)
VEC::usage = "";
POL::usage = "";
POL2::usage = "";
g::usage = "";
Mink::usage = "";
sbl::usage = "";
sbr::usage = "";
sb::usage = "";
abl::usage = "";
abr::usage = "";
ab::usage = "";
SetVec::usage = "";
GenPol::usage = "";
GenPol2::usage = "";
NDot::usage = "";
NFreeIndex::usage = "";
CutNumerics::usage = "";
CutNumericsSoft::usage = "";
(****************************************************)
(* Sewing Technology *)
(****************************************************)
H::usage = "";
HnoProp::usage = "";
MakeGravProj::usage = "";
NMakeGravProj::usage = "";
NMatNoLR::usage = "";
NMat::usage = "";
NMatH::usage = "";
ContractFull::usage = "";
(****************************************************)
(* Tensor Reduction *)
(****************************************************)
TensorReduce::usage = "";
(****************************************************)
(* Graph and Jacobi Generation *)
(****************************************************)
GetGraphBasis::usage = "";
GetGraphJacobis::usage = "";
(****************************************************)
(* fabc Color Basis *)
(****************************************************)
GetSymBasis::usage = "";
GetSymReps::usage = "";
GetDDMBasis::usage = "";
GetDDMReps::usage = "";
GetCOBasis::usage = "";
(****************************************************)
(* YM Feynman Rules *)
(****************************************************)
FeynmanYM::usage = "Individual graphs";
GetFeynYM::usage = "All n-point, L-loop topologies";
(****************************************************)
(* NLSM Feynman Rules *)
(****************************************************)
FeynNLSM::usage = "";
(****************************************************)
(* Tree Amplitudes *)
(****************************************************)
Atree::usage = "";
NLSMtree::usage = "";
Tree3::usage = "";
Mtree3::usage = "";
Mtree4::usage = "";
Mtree5::usage = "";
Mtree::usage = "";
StripPol4::usage = "";
StripPol::usage = "";
StripPol4LR::usage = "";
StripPolLR::usage = "";
(****************************************************)
(* Numerical Sewing *)
(****************************************************)
SewYM::usage = "";
SewGDA::usage = "";
SewGD::usage = "";
SewG::usage = "";
(****************************************************)
(* Analytical Cutting *)
(****************************************************)
GetCutIsos::usage = "";
CutYM::usage = "";
CutYMColor::usage = "";
CutGravity::usage = "";


(* ::Subtitle::Closed:: *)
(****************************************)
(*Functions*)
(****************************************)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Structure Definitions*)


(* Many structures thanks to Tristan Dennen *)


SetAttributes[{R,S,AngleBracket,\[Eta],ee,EE,s,t,kk,ll},Orderless];

\[Eta] /: \[Eta][\[Mu][j1__],\[Mu][j3__]] * a__[\[Mu][j3__]] := a[\[Mu][j1]];
\[Eta] /: \[Eta][\[Mu][j1__],\[Mu][j3__]] * \[Eta][\[Mu][j2__],\[Mu][j3__]] := \[Eta][\[Mu][j2],\[Mu][j1]];
\[Eta] /: \[Eta][\[Mu][j1__],\[Mu][j1__]] := Ds;
\[Eta] /: \[Eta][\[Mu][j1__],\[Mu][j2__]]^2 := Ds;
\[CapitalGamma] /: \[CapitalGamma][lorentzList1_, {\[Mu][j1__],\[Mu][j2__]}] * \[CapitalGamma][lorentzList2_,{\[Mu][j2__],\[Mu][j3__]}] := \[CapitalGamma][Join[lorentzList1,lorentzList2],{\[Mu][j1],\[Mu][j3]}];
\[CapitalGamma] /: \[CapitalGamma][lorentzList_, {\[Mu][j1__],\[Mu][j1__]}] := GammaTrace[lorentzList];
AngleBracket /: AngleBracket[a_+b_,c_] := AngleBracket[a,c]+AngleBracket[b,c];
AngleBracket /: AngleBracket[a_ b_p,c_] := a AngleBracket[b,c];
AngleBracket /: AngleBracket[a_ b_ \[Epsilon],c_] := a AngleBracket[b,c];
AngleBracket /: AngleBracket[a_ b_ \[Epsilon]2,c_] := a AngleBracket[b,c];
AngleBracket /: AngleBracket[a_ b_l,c_] := a AngleBracket[b,c];
AngleBracket /: AngleBracket[\[Epsilon][i_],\[Epsilon][j_]] := ee[i,j];
AngleBracket /: AngleBracket[\[Epsilon]2[i_],\[Epsilon]2[j_]] := EE[i,j];
AngleBracket /: AngleBracket[\[Epsilon][i_],\[Epsilon]2[j_]] := eE[i,j];
AngleBracket /: AngleBracket[\[Epsilon]2[i_],\[Epsilon][j_]] := eE[j,i];
AngleBracket /: AngleBracket[p[i_],\[Epsilon][j_]] := ke[i,j];
AngleBracket /: AngleBracket[p[i_],\[Epsilon]2[j_]] := kE[i,j];
AngleBracket /: AngleBracket[p[i_],\[Epsilon][j_]] := ke[i,j];
AngleBracket /: AngleBracket[p[i_],p[j_]] := kk[i,j];
AngleBracket /: AngleBracket[l[i_],p[j_]] := lk[i,j];
AngleBracket /: AngleBracket[l[i_],l[j_]] := ll[i,j];
AngleBracket /: AngleBracket[l[i_],\[Epsilon][j_]] := le[i,j];
AngleBracket /: AngleBracket[l[i_],\[Epsilon]2[j_]] := lE[i,j];
S /: S[a_] := AngleBracket[a,a];
(*By convention since R[a]=R[a^2] multiply by -1 if first term is negative*)
R /: R[a_] := R[-a] /; Or[Count[{a[[1]]},Times[-1,_]] > 0, Count[{a},Times[-1,_]]>0];

(*Aesthetics*)
s /: MakeBoxes[s[i__],TraditionalForm]:= ToBoxes[Interpretation[Subscript["s",StringJoin@@ToString/@List[i]],s[i]]];
R /: MakeBoxes[R[i_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["R",ToString[i]],R[i]]];
p /: MakeBoxes[p[j_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["p",ToString[j]],p[j] ]];
l /: MakeBoxes[l[j_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["l",ToString[j]],l[j]]];
\[Epsilon] /: MakeBoxes[\[Epsilon][j_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["\[Epsilon]",ToString[j]],\[Epsilon][j]]];
\[Epsilon]2 /: MakeBoxes[\[Epsilon]2[j_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["\[Epsilon]2",ToString[j]],\[Epsilon]2[j]]];
\[Mu] /: MakeBoxes[\[Mu][j_],TraditionalForm]:= ToBoxes[Interpretation[Subscript["\[Mu]",ToString[j]],\[Mu][j]]];
\[Mu] /: MakeBoxes[\[Mu][j_,i_],TraditionalForm]:= ToBoxes[Interpretation[Subsuperscript["\[Mu]",ToString[j],ToString[i]],\[Mu][j,i]]];
\[Eta] /: MakeBoxes[\[Eta][\[Mu][j1__],\[Mu][j2__]],TraditionalForm] := ToBoxes[Interpretation[Superscript["\[Eta]",Row[{\[Mu][j1],\[Mu][j2]}]],\[Eta][\[Mu][j1],\[Mu][j2]]],TraditionalForm];

(*RG /: Format[RG[a_,b_,c_]]:= "RG["<>ToString[a]<>",<>]";*)


(* ::Section::Closed:: *)
(*MomConsN, Relab, Propagate, GetVars*)


(* ::Subsection::Closed:: *)
(*MomConsN*)


(*ClearAll[MomConsN];*)
(* sQ = False: Change s's to kk's *)
MomConsN[n_,args_,sQ_:False] := 
	Module[
		{expr = args,rep,sourceVars,targetVars,pos},
		sourceVars = Variables[expr];
		targetVars = sourceVars;
		targetVars = targetVars /. {
			s[a__]:>AngleBracket[Plus@@(p/@{a}),Plus@@(p/@{a})],
			R[a_]:>AngleBracket[a,a]}; (* changes s and R to kk *)
		targetVars = targetVars //. {
			ke[n,i_]:>-Sum[ke[j,i],{j,1,n-1}],
			kE[n,i_]:>-Sum[kE[j,i],{j,1,n-1}],
			kk[n,i_]:> -Sum[kk[j,i],{j,1,n-1}],
			lk[i_,n]:> -Sum[lk[i,j],{j,1,n-1}]
		}; (* removes p[n] *)
		rep = Join[
			{ke[j1_,j1_]:>0,
				ke[n-1,n]:>-Sum[ke[i,n],{i,1,n-2}],kE[j1_,j1_]:>0,
				kE[n-1,n]:>-Sum[kE[i,n],{i,1,n-2}]},
			Flatten@Table[kk[i,j] -> (s@@Range[i,j] - s@@Range[i+1,j] - s@@Range[i,j-1] + s@@Range[i+1,j-1])/2,
				{j,2,n-1},{i,1,j-1}],
			Flatten@Table[kk[i,i]->s[i],{i,1,n-1}]
		];
		rep = rep //. {s[]->0, s[a_]:>0, s@@Range[1,n-1]:>0};
		targetVars = targetVars //. rep;
		(*Change s[i,j]:>2kk[i,j]*)
		If[!sQ,
			targetVars = (targetVars /. s[a__]:>AngleBracket[Plus@@(p/@{a}),Plus@@(p/@{a})])/.Table[kk[i,i]->0,{i,1,n}] //Expand;
		];
		pos = Flatten@Position[Table[SameQ[sourceVars[[i]],targetVars[[i]]],{i,Length[sourceVars]}],False,Heads->False];
		sourceVars = sourceVars[[pos]];
		targetVars = targetVars[[pos]];
		rep = Table[sourceVars[[i]]->targetVars[[i]],{i,Length[sourceVars]}];
		expr //. rep
	];


(* ::Subsection::Closed:: *)
(*Relab*)


(* Example: Relab[{p[1]->p[2],\[Epsilon][1]->\[Epsilon][2],...},ke[1,2]] *)
ClearAll[tAngleBracket];
tAngleBracket /: tAngleBracket[a_+b_,c_] := tAngleBracket[a,c]+tAngleBracket[b,c];
tAngleBracket /: tAngleBracket[a_ b_p,c_] := a tAngleBracket[b,c];
tAngleBracket /: tAngleBracket[a_ b_ \[Epsilon],c_] := a tAngleBracket[b,c];
tAngleBracket /: tAngleBracket[a_ b_ \[Epsilon]2,c_] := a tAngleBracket[b,c];
tAngleBracket /: tAngleBracket[a_ b_l,c_] := a tAngleBracket[b,c];

(*Clear[Relab];*)
Relab[arg2_,arg3_] := 
	Module[
		{reps = arg2,expr = arg3,sourceVars,targetVars,pos},
		sourceVars = Variables[expr];
		targetVars = sourceVars /. {
			ee[i_,j_]:>tAngleBracket[\[Epsilon][i],\[Epsilon][j]],
			ke[i_,j_]:>tAngleBracket[p[i],\[Epsilon][j]],
			le[i_,j_]:>tAngleBracket[l[i],\[Epsilon][j]],
			EE[i_,j_]:>tAngleBracket[\[Epsilon]2[i],\[Epsilon]2[j]],
			kE[i_,j_]:>tAngleBracket[p[i],\[Epsilon]2[j]],
			lE[i_,j_]:>tAngleBracket[l[i],\[Epsilon]2[j]],
			eE[i_,j_]:>tAngleBracket[\[Epsilon][i],\[Epsilon]2[j]],
			kk[i_,j_]:>tAngleBracket[p[i],p[j]],
			lk[i_,j_]:>tAngleBracket[l[i],p[j]],
			ll[i_,j_]:>tAngleBracket[l[i],l[j]],
			s[a__] :> tAngleBracket[Plus@@(p/@{a}),Plus@@(p/@{a})]};
		targetVars = targetVars /. reps;
		targetVars = targetVars /. tAngleBracket->AngleBracket;
		
		pos = Flatten@Position[Table[SameQ[sourceVars[[i]],targetVars[[i]]],{i,Length[sourceVars]}],False,Heads->False];
		sourceVars = sourceVars[[pos]];
		targetVars = targetVars[[pos]];
		reps = Table[sourceVars[[i]]->targetVars[[i]],{i,Length[sourceVars]}];
		expr = expr /. reps;
		expr
	];


(* ::Subsection::Closed:: *)
(*Propagate*)


(*Clear[Propagate];*)
Propagate[args_] := 
	Module[
		{expr = Expand[args]},
		expr = expr //. {a_[\[Mu][j1__]]b_[\[Mu][j1__]] :> AngleBracket[a,b], a_[\[Mu][j1__]]^2 :> AngleBracket[a,a]};
		expr
	];


(* ::Subsection::Closed:: *)
(*GetVars*)


(*Clear[GetVars];*)
GetVars[expr_,removeList_:{}]:=Module[{list=removeList,vars},
	vars=Variables[expr];
	While[Length[list]>0,
		vars=DeleteCases[vars,list[[1]]];
		list=Delete[list,1];
	];
	Sort[vars]
];


(*Example*)
(*
GetVars[numList[[1]],{Ds,a1[_],ke[1,2],lk[6,1]}]
*)


(* ::Section::Closed:: *)
(*Linear Equation Solver*)


(* From Henrik and/or Radu and/or John Joseph and/or Zvi *)


(* ::Subsection:: *)
(*MonsterSolve2*)


MonsterSolve2[eqns_,chunksize_,minsize_]:=Module[{},
	Print["MonsterSolve is a propriety function. Contact me for more information."];
];


(* ::Section::Closed:: *)
(*Graph Information*)


(* ::Subsection::Closed:: *)
(*DisplayGraph*)


DisplayGraph[verts_List] := Module[{g = BuildGraph[verts]},
	If[$VersionNumber>=8 && Length[Union[Sort/@g[[2]]]]==Length[g[[2]]],
		System`Graph[Labeled@@@({DirectedEdge@@@g[[2]],g[[3]] /. a_[i_]:>Subscript[a,i]}\[Transpose]),GraphStyle->"DynamicHighlight",GraphLayout->"SpringElectricalEmbedding"]
	(*	System`Graph[Labeled@@@({DirectedEdge@@@g[[2]],Range[Length[g[[2]]]]}\[Transpose]),GraphStyle->"DynamicHighlight"]*)
	,
		GraphPlot[{Rule@@@g[[2]],g[[3]]}\[Transpose],
			VertexLabeling->False,
			DirectedEdges->True,
			MultiedgeStyle->1,
			EdgeLabeling->True(*Automatic*),
			Method->"SpringEmbedding",
			ImageSize->Small,
			PlotStyle->PointSize[Medium]]
	]
];


(* ::Subsection::Closed:: *)
(*GraphInfo*)


(*Clear[GraphInfo];*)
GraphInfo[g_List]:= Module[{extEdges,intEdges,verts},	
	extEdges = Flatten[Abs[Select[g, Length[#]==1&]]];
	verts = Select[g, Length[#]>1&];
	intEdges = Complement[Flatten[Abs[verts]],extEdges];
	{extEdges,intEdges,verts}
];


(* ::Subsection::Closed:: *)
(*MomRules*)


(*Clear[MomRules];*)
MomRules[g_List,loopChoices_:{}]:=Module[{extEdges,intEdges,verts,eqs,solveFor,soln},
	{extEdges,intEdges,verts} = GraphInfo[g];
	eqs = Plus@@((p/@#)/.p[i_]:>Sign[i]p[Abs[i]])==0&/@verts;
	solveFor = Join[{p[extEdges[[-1]]]},p/@Complement[intEdges,loopChoices]];
	soln = Solve[eqs,solveFor][[1]];//Quiet;
	soln
];


(* ::Subsection::Closed:: *)
(*ChooseLoopMom*)


(*Out: {{loopEdge1, loopEdge2, ... }, {maxPower1, maxPower2, ...}}*)
ChooseLoopMom[g_]:=Module[
{intVerts,extEdges,intEdges,edges,edgePos,nPoint,
momrules,momcons,indep,loopNum,
mat,subs,indepSubs,
candidates,choice,loopmom,
min,small,minSum,smaller},

	(**************************************************************)
	(*Graph info*)
	{extEdges,intEdges,intVerts} = GraphInfo[g];
	edges = Join[extEdges,intEdges];
	edgePos[edge_]:=Position[edges,edge][[1,1]];
	nPoint = Length[extEdges];
	(**************************************************************)

	(**************************************************************)
	(*Momentum Conservation*)
	momrules = MomRules[g];
	momcons = (p/@edges)/.momrules;
	indep = GetVars[momcons,{}];
	loopNum = Length[indep] - (nPoint - 1);
	(**************************************************************)

	(**************************************************************)
	(*Find all independent sets of loop momenta*)
	mat=Table[Coefficient[momcons[[i]],indep[[j]]],{i,1,Length[edges]},{j,1,Length[indep]}];
	subs=Join[extEdges[[1;;Length[extEdges]-1]],#]&/@Subsets[Table[i,{i,intEdges}],{loopNum}];
	(*independent loop momenta*)
	indepSubs=subs[[Position[Table[Det[mat[[edgePos/@(subs[[i]])]]]!=0,{i,1,Length[subs]}],True][[All,1]]]];
	(**************************************************************)

	(**************************************************************)
	(*Find n-gons for each possible set of loop legs*)
	candidates=Range[Length[indepSubs]];
	For[choice=1,choice<=Length[indepSubs],choice++,
		(*If[Mod[choice,200]===0,Print[choice]];*)
	
		momcons=Solve[momrules/.Rule->Equal,p/@Complement[edges,indepSubs[[choice]]]];
		If[momcons==={},
			candidates[[choice]]=Null;,

			momcons=momcons[[1]];
			momcons=(p/@intEdges)/.momcons;
			loopmom=Select[GetVars[momcons,{}]/.p[i_]:>i,#>nPoint&];

			candidates[[choice]]=
			{loopmom,
			Length[Position[momcons,#][[All,1]]]&/@loopmom};
		];
	];
	(**************************************************************)

	(**************************************************************)
	(*Select loop momenta with good power counting*)
	min=Min[Max/@(candidates[[All,2]])];
	small=Select[candidates,Max[#[[2]]]===min&];

	minSum=Min[Plus@@#&/@small[[All,2]]];
	smaller=Select[candidates,Plus@@(#[[2]])===minSum&];
	(**************************************************************)

	smaller[[1]]
];


(* ::Subsection::Closed:: *)
(*CutSolution (All exposed internal edges are put on-shell)*)


CutSolution[verts_List, loopEdges0_List:{}] := Module[
{loopEdgesRecommended,loopEdges,extEdges,n,intEdges,edges,intVerts,cutEdges,edgePos,
momConsRep,momList,onShellEdges,onShellEqs,toElimVars,cutSoln,kks,kes,kin,reps,rule,x},

	(*************************************)
	(*Set Loop Momenta*)
	(*************************************)
	(* Make sure we have enough loop momenta *)
	(* Still, if the user enters the correct number of loop momenta
		which are not independent, then bad things will happend *)
	loopEdgesRecommended = ChooseLoopMom[verts][[1]];
	If[ Length[loopEdges0] != Length[loopEdgesRecommended],
		loopEdges = loopEdgesRecommended;,
		loopEdges = loopEdges0;
	];
	
	(*************************************)
	(*Graph Information*)
	(*************************************)
	{extEdges,intEdges,intVerts} = GraphInfo[verts];
	extEdges = Sort[extEdges];
	n = Length[extEdges];
	If[ extEdges != Range[n],
		Print["The external edges are not labelled 1,2,...,n, which the function 'CutSolution' demands."];
		Return;
	];
	intEdges = Sort[intEdges];
	edges = Join[extEdges,intEdges];
	cutEdges = intEdges;
	(*Since we don't expect every edge to appear in our contact graphs, this function will tell you the position of the edge number*)
	(*eg. edges={1,2,3,4,5,7,9,10,11}, edgePos[7]=6*)
	edgePos[edge_] := Position[edges,edge][[1,1]];

	(*************************************)
	(*Momentum Conservation*)
	(*************************************)
	momConsRep = MomRules[verts, loopEdges];
	momList = (p/@edges)/.momConsRep;

	(*************************************)
	(*On-Shell*)
	(*************************************)
	onShellEdges = Sort[Join[extEdges,cutEdges]];
	onShellEqs = (AngleBracket[#,#]&/@(momList[[edgePos/@onShellEdges]]));
	onShellEqs = Thread[Expand[MomConsN[n,onShellEqs]]==0];
	toElimVars = Flatten[Join[
				{kk[1,n-1]},(*convention of MomConsN to eliminate this*)
				Table[kk[i,i],{i,1,n-1}],
				Table[kk[loopEdges[[i]],j],{i,1,Length[loopEdges]},{j,1,n-1}],
				Table[kk[loopEdges[[i]],loopEdges[[j]]],{i,1,Length[loopEdges]},{j,i,Length[loopEdges]}]
			]];
	Quiet[cutSoln = ToRules[Reduce[onShellEqs,toElimVars]];];
	
	(*************************************)
	(*Kinematic Variable Replacements*)
	(*************************************)
	kks = Flatten[Table[kk[edges[[i]],edges[[j]]],{i,1,Length[edges]},{j,i,Length[edges]}]];
	reps = rule[
				#, 
				Expand[
					MomConsN[n,
						Relab[momConsRep, #]
					]/. cutSoln
				]
			]&/@kks;
	reps = DeleteCases[reps,rule[x_,x_]]/.rule->Rule;

	reps
];


(* ::Subsection::Closed:: *)
(*GetCutReps (Careful with reference momenta replacements. Specialized!)*)


GetCutReps[verts_List, loopEdges_List:{}, projLegChoices_List:{}]:=Module[
{extEdges,intEdges,intVerts,n,edges,cutEdges,
kkReps,momReps,rule,x,sReps,Rreps,eeReps,kes,keReps,
refs,keCantAppear,kkCantAppear,refCantAppear,kkRefReps,cutReps},

	(*************************************)
	(*Graph Information*)
	(*************************************)
	{extEdges,intEdges,intVerts} = GraphInfo[verts];
	extEdges = Sort[extEdges];
	n = Length[extEdges];
	If[ extEdges != Range[n],
		Print["The external edges are not labeled 1,2,...,n, which the function 'CutSolution' demands."];
		Return;
	];
	intEdges = Sort[intEdges];
	edges = Join[extEdges,intEdges];
	cutEdges = intEdges;

	(*************************************)
	(*Cut and Mom Cons Replacements*)
	(*************************************)
	kkReps = CutSolution[verts,loopEdges];
	momReps = MomRules[verts, loopEdges];

	(*************************************)
	(*Additional Replacements*)
	(*************************************)
	sReps = rule[#, 
				Expand[MomConsN[n,#]/.kkReps]
			]&/@(s@@@Subsets[Range[n],{2,n}]);
	sReps = DeleteCases[sReps,rule[x_,x_]]/.rule->Rule;
	
	Rreps=DeleteDuplicates[Table[
			If[MemberQ[cutEdges,i],
				{rule[R[p[i]], 1(*I*)], rule[R[Expand[p[i]/.momReps]], 1(*I*)]},
				rule[R[p[i]], R[Expand[p[i]/.momReps]]]
			],
		{i,intEdges}
		]//Flatten];
	Rreps = DeleteCases[Rreps,rule[x_,x_]]/.rule->Rule;

	eeReps = Table[ee[i,i]->0,{i,extEdges}];
	eeReps = Join[eeReps, eeReps/.ee->EE];

	kes = Flatten[Table[ke[edges[[i]],extEdges[[j]]],{i,1,Length[edges]},{j,1,Length[extEdges]}]];
	keReps = rule[
				#, 
				Expand[
					MomConsN[n,
						Relab[momReps, #]
					]
				]
			]&/@kes;
	keReps = DeleteCases[keReps,rule[x_,x_]]/.rule->Rule;
	keReps = Join[keReps, keReps/.ke->kE];

	(*reference momenta terms that can't cancel light-cone denominator*)
	refs = ToExpression["ref"<>ToString[#]]&/@projLegChoices;
	keCantAppear = Flatten[Table[ke[#,extEdges[[i]]]->0,{i,1,Length[extEdges]}]&/@refs];
	keCantAppear = Join[keCantAppear,keCantAppear/.{ke->kE}];
	kkCantAppear = Join[
					Flatten[Table[kk[#,extEdges[[i]]]->0,{i,1,Length[extEdges]}]&/@refs],
					Thread[kk@@#&/@DeleteDuplicates[Sort/@Tuples[refs,2]]->0],
					Thread[Apply[kk,Complement[Tuples[{projLegChoices,refs}],{#,ToExpression["ref"<>ToString[#]]}&/@projLegChoices],1]->0]
				];
	refCantAppear=Join[keCantAppear,kkCantAppear];
	kkRefReps=Join@@Table[kk[#,refs[[i]]]->(Relab[momReps,kk[#,refs[[i]]]]//.refCantAppear)&/@Complement[intEdges,projLegChoices],{i,1,Length[refs]}];

	(*************************************)
	(*Combine and Return*)
	(*************************************)
	cutReps = Dispatch[Expand[Join[kkReps,sReps,keReps,eeReps,Rreps,refCantAppear,kkRefReps]]];

	cutReps
];


(* ::Section::Closed:: *)
(*Graph Isomorphisms*)


(* Many functions written by Tristan Dennen *)


(* ::Subsection::Closed:: *)
(*BuildGraph*)


(*ClearAll[BuildGraph];*)
BuildGraph[vertices_List]:=
	Module[
		{edges,labels},
		labels = Sort[DeleteDuplicates[Abs/@Flatten[vertices]]];
		edges = Table[{Position[vertices,i][[1,1]],Position[vertices,-i][[1,1]]},{i,labels}];
		labels = p/@labels;
		BG[vertices,edges,labels]
	];


(* ::Subsection::Closed:: *)
(*EdgewiseAutomorphisms*)


(*ClearAll[EdgewiseAutomorphisms];*)
EdgewiseAutomorphisms[vertices_List] := 
	Module[
		{g1,g2,
			e1,e2,v1,v2,
			bubbles,
			automorphisms,
			sortedPermutedE1,permutedE1,signaturePermutedE1,sortedE2,signatureE2,
			tmp,len,ind,target,source,i,
			labels,repsToEdgeNumbers,
			signatures,edgewiseAutomorphisms},
		g1 = BuildGraph[vertices];
		g2 = g1;
		v1 = g1[[1]]; v2 = g2[[1]];
		e1 = g1[[2]]; e2 = g2[[2]];
		bubbles = { Length[e1] != Length@Union[Sort/@e1],Length[e2]!=Length@Union[Sort/@e2] };
		(* Need to use Combinatorica to make a list of ALL automorphisms. mma8 only gives one. *)
		automorphisms = Combinatorica`Automorphisms[Combinatorica`FromUnorderedPairs[e1]];
		If[Length[automorphisms]==0,Return[{{},{}}]];
		automorphisms = Rule@@@({Range[Length[#]],#}\[Transpose])&/@automorphisms;
		permutedE1 = e1 /. automorphisms;
		sortedPermutedE1 = Map[Sort,permutedE1,{2}];
		sortedE2 = Map[Sort,e2,{1}];
		edgewiseAutomorphisms = {};
		For[i=1,i<=Length[automorphisms],i++,
			tmp = {Position[sortedE2,#]&/@sortedPermutedE1[[i]]};
	(****************************************************************************************************************)
	(* Bubbles need extra care: *)
	(****************************************************************************************************************)
			While[Max[(len = Map[Length,tmp,{2}])]>1, (* Permutation isn't one-to-one *)
				ind = Position[len,2][[1]]; (* Find the first instance of a sortedPermutedE1 edge that goes two places *)
				target = tmp[[ind[[1]],ind[[2]]]]; (* Positions in sortedE2 *)
				source = Position[tmp[[ind[[1]]]],target];
				tmp = ReplacePart[tmp,ind[[1]]->
					Sequence[ReplacePart[tmp[[ind[[1]]]],{source[[1]]->{target[[1]]},source[[2]]->{target[[2]]}}],
					ReplacePart[tmp[[ind[[1]]]],{source[[1]]->{target[[2]]},source[[2]]->{target[[1]]}}]]]; (* Split into the two possible permutations *)
			];
	(****************************************************************************************************************)
	(* Figure out signs of permuted edges *)
	(****************************************************************************************************************)
			tmp = Flatten/@tmp;
			signaturePermutedE1 = Map[Signature,permutedE1[[i]],{1}];
			signatureE2 = Map[Signature,e2[[#]]&/@tmp,{2}];
			tmp = Table[tmp[[j]]*signatureE2[[j]]*signaturePermutedE1,{j,Length[tmp]}];
			edgewiseAutomorphisms = Join[edgewiseAutomorphisms,tmp];
		];
	(****************************************************************************************************************)
	(* Determine the signature of the permutation. That is, g1->g2 or g1->-g2 under this permutation? *)
	(****************************************************************************************************************)
		labels=Sort[DeleteDuplicates[Abs[Flatten[v1]]]];
		repsToEdgeNumbers=Thread[labels->Range[Length[labels]]];
		signatures=Table[
			(Times@@(Signature/@Map[Sign[#]*edgewiseAutomorphisms[[i1,Abs[#]/.repsToEdgeNumbers]]&,v1,{2}])) *
			(Times@@(Signature/@v2)),
			{i1,1,Length[edgewiseAutomorphisms]}
		];
		{signatures,edgewiseAutomorphisms}
	];

EdgewiseAutomorphisms[bg_BG] := EdgewiseAutomorphisms[bg[[1]]];


(* ::Subsection::Closed:: *)
(*EdgewiseIsomorphism*)


(*ClearAll[EdgewiseIsomorphism];*)
EdgewiseIsomorphism[vertices1_List,vertices2_List] := 
	Module[
		{g1,g2,
			e1,e2,v1,v2,
			bubbles,
			isomorphisms,
			sortedPermutedE1,permutedE1,signaturePermutedE1,sortedE2,signatureE2,
			tmp,len,ind,target,source,i,
			labels,repsToEdgeNumbers,
			signatures,edgewiseIsomorphisms},
		g1 = BuildGraph[vertices1];
		g2 = BuildGraph[vertices2];
		v1 = g1[[1]]; v2 = g2[[1]];
		e1 = g1[[2]]; e2 = g2[[2]];
		bubbles = { Length[e1] != Length@Union[Sort/@e1],Length[e2]!=Length@Union[Sort/@e2] };
		If[$VersionNumber===8 && Nor@@bubbles, (* Do the faster check in newer versions. Doesn't cost anything. *)
			isomorphisms = FindGraphIsomorphism[System`Graph[UndirectedEdge@@@e1],System`Graph[UndirectedEdge@@@e2]];
			isomorphisms = Range[Length[isomorphisms]] /. isomorphisms;,
			isomorphisms = Combinatorica`Isomorphism[Combinatorica`FromUnorderedPairs[e1],Combinatorica`FromUnorderedPairs[e2]];
		];
		(* Need to use Combinatorica to make a list of ALL isomorphisms. mma8 only gives one isomorphism. *)
		If[Length[isomorphisms]==0,Return[{{},{}}]];
		isomorphisms = {isomorphisms};
		isomorphisms = Rule@@@({Range[Length[#]],#}\[Transpose])&/@isomorphisms;
		permutedE1 = e1 /. isomorphisms;
		sortedPermutedE1 = Map[Sort,permutedE1,{2}];
		sortedE2 = Map[Sort,e2,{1}];
		edgewiseIsomorphisms = {};
		For[i=1,i<=Length[isomorphisms],i++,
			tmp = {Position[sortedE2,#]&/@sortedPermutedE1[[i]]};
	(****************************************************************************************************************)
	(* This is good enough most of the time, but we have to deal with bubbles: *)
	(****************************************************************************************************************)
			While[Max[(len = Map[Length,tmp,{2}])]>1, (* Permutation isn't one-to-one *)
				ind = Position[len,2][[1]]; (* Find the first instance of a sortedPermutedE1 edge that goes two places *)
				target = tmp[[ind[[1]],ind[[2]]]]; (* Positions in sortedE2 *)
				source = Position[tmp[[ind[[1]]]],target];
				tmp = ReplacePart[tmp,ind[[1]]->
					Sequence[ReplacePart[tmp[[ind[[1]]]],{source[[1]]->{target[[1]]},source[[2]]->{target[[2]]}}],
					ReplacePart[tmp[[ind[[1]]]],{source[[1]]->{target[[2]]},source[[2]]->{target[[1]]}}]]]; (* Split into the two possible permutations *)
			];
	(****************************************************************************************************************)
	(* Now need to figure out signs of permuted edges *)
	(****************************************************************************************************************)
			tmp = Flatten/@tmp;
			signaturePermutedE1 = Map[Signature,permutedE1[[i]],{1}];
			signatureE2 = Map[Signature,e2[[#]]&/@tmp,{2}];
			tmp = Table[tmp[[j]]*signatureE2[[j]]*signaturePermutedE1,{j,Length[tmp]}];
			edgewiseIsomorphisms = Join[edgewiseIsomorphisms,tmp];
		];
	(****************************************************************************************************************)
	(* Determine the signature of the permutation. That is, g1->g2 or g1->-g2 under this permutation? *)
	
	(* Map[Sign[#]*edgewiseIsomorphisms[[i1,Abs[#]]]&,v1,{2}] makes a new set of vertices by taking the *)
	(* original set and replacing the edges by the rules of edgewiseIsomorphisms *)
	(****************************************************************************************************************)
		labels=Sort[DeleteDuplicates[Abs[Flatten[v1]]]];
		repsToEdgeNumbers=Thread[labels->Range[Length[labels]]];
		signatures=Table[
			(Times@@(Signature/@Map[Sign[#]*edgewiseIsomorphisms[[i1,Abs[#]/.repsToEdgeNumbers]]&,v1,{2}])) *
			(Times@@(Signature/@v2)),
			{i1,1,Length[edgewiseIsomorphisms]}
		];
		{signatures,edgewiseIsomorphisms}
	];

EdgewiseIsomorphism[bg1_BG,bg2_BG] := EdgewiseIsomorphism[bg1[[1]],bg2[[1]]];


(* ::Subsection::Closed:: *)
(*LabelwiseAutomorphisms*)


(*ClearAll[LabelwiseAutomorphisms];*)
LabelwiseAutomorphisms[vertices_List] := 
	Module[
		{signature,iso,ordering},
		{signature,iso} = EdgewiseAutomorphisms[vertices];
		If[Length[iso]==0,Return[{{},{}}]];
		iso = LabelwiseIsomorphism[vertices,vertices,#]&/@iso;
		ordering = Ordering[iso];
		iso = iso[[ordering]];
		signature = signature[[ordering]];
		{signature,iso}
	];

LabelwiseAutomorphisms[bg_BG] := LabelwiseAutomorphisms[bg[[1]]];


(* ::Subsection::Closed:: *)
(*LabelwiseIsomorphism*)


(*ClearAll[LabelwiseIsomorphism];*)
LabelwiseIsomorphism[verts1_List,verts2_List] := 
	Module[
		{signature,iso,ordering},
		
		{signature,iso} = EdgewiseIsomorphism[verts1,verts2];
		If[Length[iso]==0,Return[{{},{}}]];
		iso = LabelwiseIsomorphism[verts1,verts2,#]&/@iso;
		ordering = Ordering[iso];
		iso = iso[[ordering]];
		signature = signature[[ordering]];
		{signature,iso}
	];
LabelwiseIsomorphism[verts1_List,verts2_List,arg3_List] :=
	Module[
		{source = BuildGraph[verts1], target = BuildGraph[verts2], eAut = arg3,
			autSigns,autPerms,sourceLabels,sourceExternalEdges,n,e,
			pRelabelings,epsRelabelings,eps2Relabelings,lRelabelings},
		autSigns = Sign[eAut];
		autPerms = Abs[eAut];
		sourceLabels = source[[3]];
	(* Identify which edges are external. *)
		sourceExternalEdges = Abs@Flatten@Cases[source[[1]],a_ /; Length[a]==1];
		n = Length[sourceExternalEdges];
		e = Length[source[[2]]];
	(* make list of relabels *)
		pRelabelings = Table[source[[3,i]] -> autSigns[[i]] target[[3, autPerms[[i]] ]],{i,1,n}];
		epsRelabelings = Table[(source[[3,i]]/.p->\[Epsilon]) -> (target[[3, autPerms[[i]] ]]/.p->\[Epsilon]),{i,1,n}];
		lRelabelings = Table[source[[3,i]] -> autSigns[[i]] target[[3, autPerms[[i]] ]],{i,n+1,e}];
		Join[pRelabelings,lRelabelings,epsRelabelings]
	];

(* For the old BG structure *)
LabelwiseIsomorphism[bg1_BG,bg2_BG] := LabelwiseIsomorphism[bg1[[1]],bg2[[1]]];
LabelwiseIsomorphism[bg1_BG,bg2_BG,arg3_List] := LabelwiseIsomorphism[bg1[[1]],bg2[[1]],arg3];


(* ::Subsection::Closed:: *)
(*Pre-Isomorphism Invariants*)


(*From Tristan's BCJFunctions*)
PreIsoInvariant[graph0_Graph] := 
	Module[
		{g = graph0,degrees,minDegree,numMin,posMin,invariant},
		invariant = "";
		While[V[g]>0,
			degrees = Length/@Combinatorica`ToAdjacencyLists[g];
			minDegree = Min[degrees];
			posMin = Flatten@Position[degrees,minDegree];
			numMin = Length[posMin];
			invariant = StringJoin[invariant, Table[ToString[minDegree],{numMin}],"3"];
			g = Combinatorica`DeleteVertices[g,posMin];
		];
		invariant
	];

PreIsoInvariant[arg_BG] := PreIsoInvariant[Combinatorica`FromUnorderedPairs[arg[[2]]]];


(* ::Subsection::Closed:: *)
(*isIso*)
(*(Input: vertex list,  basis vertex list)*)


(*Clear[isIso];*)
isIso[graph_List, basis_List]:=Module[
{ib,labels,edges,graphUnOrd,ibGraphList,baseInvariants,myInvariant,toCheck,i,isoQ},

	ib = BuildGraph/@basis;

	(*Make unordered graphs from vertices list structure*)
	labels = Sort[DeleteDuplicates[Abs/@Flatten[graph]]];
	edges = Table[{Position[graph,i][[1,1]],Position[graph,-i][[1,1]]},{i,labels}];
	graphUnOrd = Combinatorica`FromUnorderedPairs[edges];

	(*The unordered graphs for the INTEGRALBASIS*)
	ibGraphList = Table[Combinatorica`FromUnorderedPairs[ib[[i,2]]],{i,1,Length[ib]}];
	
	baseInvariants = PreIsoInvariant[#]&/@ibGraphList;
	myInvariant = PreIsoInvariant[graphUnOrd];
	toCheck = Flatten@Position[baseInvariants,myInvariant,{1},Heads->False];

	i = 1;
	isoQ = False;
	While[And[isoQ === False, i <= Length[toCheck]],
		isoQ = Combinatorica`IsomorphicQ[graphUnOrd,ibGraphList[[toCheck[[i]]]]];
		i = i+1;
	];
	(* If we found it, also return the index of the graph it's isomorphic to*)
	Return[
		If[isoQ===False,
			{isoQ,{}},
			{isoQ, toCheck[[i-1]]}
		]
	];
];


(* ::Subsection::Closed:: *)
(*FindIso*)
(*(Input: vertex list,  basis vertex list)*)


Clear[InverseReps];
InverseReps[arg_List] :=
	Module[
		{reps = arg,sourceLabels,targetLabels,sourceVars,targetVars,
			relabelings,eqns,solns},
		targetLabels = reps[[All,1]];
		sourceLabels = reps[[All,2]];
		sourceVars = Variables[sourceLabels];
		targetVars = Variables[targetLabels];
(*		If[Length[sourceVars]!=Length[targetVars],Message[InverseReps::NonInvertibleSystem];];*)
		relabelings = Table[sourceVars[[i1]]->\[Beta][i1],{i1,Length[sourceVars]}];
		eqns = Table[(sourceLabels[[i1]] /. relabelings) == targetLabels[[i1]],{i1,Length[targetLabels]}];
		solns = ToRules@Reduce[eqns,Array[\[Beta],Length[sourceVars]]];
		solns = Cases[solns,a_/;Head[a[[1]]]===\[Beta]];
		relabelings = relabelings /. solns;
		relabelings
	];


(*Find the basis element that is isomorphic to our input vertex list*)
(*InverseQ = False: returns relabeling of basis to graph*)
(*InverseQ = True: returns relabeling of graph to basis*)
(*Clear[FindIso];*)
FindIso[graph_List, basis_List, InverseQ_:False]:=Module[
{isoQ,isoIndex,isoStruc,isoSig,iso},

	{isoQ,isoIndex} = isIso[graph,basis];

	(****************************************************************)
	(*If it can't find isomorphism, the return {0,0,{}}*)
	If[isoQ===False,
		{0,0,{}},
	
		(*Find the actual isomorphism*)
		iso = Transpose[LabelwiseIsomorphism[basis[[isoIndex]],graph]][[1]];
		isoSig = iso[[1]];
		iso = Sort[iso[[2]]/.l->p];
		If[InverseQ === True,
			iso = InverseReps[iso];
		];
		{isoIndex,isoSig,iso}
	]
];


(* ::Section::Closed:: *)
(*Graph and Jacobi Generation*)


(* A lot of code from Tristan Dennen *)


(* ::Subsection::Closed:: *)
(*Jacobi*)


Default[Jacobi,3]=MINGIRTH;
Jacobi[arg1_BG,arg2_Integer,mingirth_.] := 
	Module[
		{g1 = arg1,
			g2 = arg1,
			g3 = arg1,
			e=arg2,
			v,edges,m,momenta},
		v = Join[Position[g1[[1]],e],Position[g1[[1]],-e]];
		m = Max@Flatten@Abs[g1[[1]]];
	(* Cycle the two vertices until the Jacobi edge is the first entry in each. *)
		edges = Flatten[{RotateLeft[g1[[1]][[v[[1,1]]]],v[[1,2]]-1],
			RotateLeft[g1[[1]][[v[[2,1]]]],v[[2,2]]-1]}];
		momenta = Sign[edges] g1[[3,Abs[edges]]];
	(* Build the two other graphs by hand. *)
		g2[[1]][[v\[Transpose][[1]]]] = {edges[[{1,2,5}]],edges[[{4,6,3}]]};
		g3[[1]][[v\[Transpose][[1]]]] = {edges[[{1,2,6}]],edges[[{4,3,5}]]};

		(*Screws up signs in GenerateJacobis I think*)	
		(*g2[[1]][[v\[Transpose][[1]]]] = {edges[[{1,2,5}]],edges[[{4,3,6}]]}; (*u-channel - Changed 10/12/2012 - Josh Nohle. Using this, BCJ is g1-g2-g3=0 and g3 is color-ordered. (g2 might not be for nonplanar.)*)*)
		(*g3[[1]][[v\[Transpose][[1]]]] = {edges[[{1,6,2}]],edges[[{4,3,5}]]}; (*s-channel - Changed 10/12/2012 - Josh Nohle*)*)
		
		g2[[2]] = Table[{Position[g2[[1]],i][[1,1]],Position[g2[[1]],-i][[1,1]]},{i,m}];
		g3[[2]] = Table[{Position[g3[[1]],i][[1,1]],Position[g3[[1]],-i][[1,1]]},{i,m}];

	(* 
		We will eventually want to solve this Jacobi identity for one of the three graphs. 
		To do this, it is important to maintain invertibility of the relabeling from the base graphs.
		To this end, we need to make sure that the momentum on the original Jacobi leg appears in all three graphs. 
		This momentum assignment is necessary unless we change how GenerateJacobis[] works. 
	*)
		g2[[3,e]] = momenta[[1]]+momenta[[3]]-momenta[[5]];
		g3[[3,e]] = momenta[[1]]+momenta[[3]]-momenta[[6]];

		If[mingirth==0,
			{g1,g2,g3},
			Cases[{g1,g2,g3},a_ /; Girth[Combinatorica`FromUnorderedPairs[a[[2]]]]>=mingirth,1]
		]
	];


(* ::Subsection::Closed:: *)
(*SortRG*)


SortRG[expr_] := expr /. a_RG :>SortRG[a];
SortRG[arg_RG] :=
	Module[
		{ind = arg[[1]],sig = arg[[2]],reps = Expand[arg[[3]]],ord},
		reps = Sort/@reps;
		ord = Ordering[reps];
		reps = reps[[ord]];
		sig = sig[[ord]];
		RG[ind,sig,reps]
	];


(* ::Subsection::Closed:: *)
(*PlanarGraphQ*)


ClearAll[PlanarGraphQQ];
PlanarGraphQQ[arg_BG] := 
	Module[
		{g = arg[[2]],externalV},
		externalV = Cases[Tally[Flatten[g]],a_ /; a[[2]]==1,{1},Heads->False][[All,1]];
		g = g /. {a_ /; MemberQ[externalV,a] :> externalV[[1]]};
		PlanarQ@Combinatorica`FromUnorderedPairs[g]
	];


(* ::Subsection::Closed:: *)
(*GenerateJacobis*)


ClearAll[BuildGraphGJ];
(* Old version of build graph. *)
(* The new version is more robust with regard to labels. *)
(* They should give the same result as long as the labels go from 1 to Max[Flatten[vertices]] *)
BuildGraphGJ[vertices_]:=
	Module[
		{m,edges,labels},
		m = Max@Flatten[vertices];
		edges = Table[{Position[vertices,i][[1,1]],Position[vertices,-i][[1,1]]},{i,m}];
		labels = AssignLabelsNoMomCons[vertices];
		BG[vertices,edges,labels]
	];


ClearAll[AssignLabels,AssignLabelsNoMomCons];
AssignLabels[args_] := 
	Module[
		{g=args,
			edges,v,localL,n,m,i,
			externalV,internalV,externalE,internalE,
			momcons,pos,sign,soln},
		m = Max@Flatten[g];
		edges = Table[{Position[g,-i][[1,1]],Position[g,i][[1,1]]},{i,m}];
		v = Length/@g;
		externalV = Flatten@Position[v,1]; (* External vertices *)
		internalV = Flatten@Position[v,3]; (* Internal vertices *)
		externalE = Abs[g[[externalV]]];
		internalE = Complement[Range[m],externalE];
	
		(* Set up momentum conservation constaints *)
		momcons = (Plus@@(localL/@(g[[#]])) == 0) &/@internalV;
	
		(* Label external lines *)
		localL[arg_] /; (arg < 0) := -localL[-arg]; 
		localL[arg_] /; (MemberQ[externalV,arg]) := p[arg]; 
		localL[arg_] := l[arg];
	
		soln = Reduce[momcons,localL/@internalE] /. {Equal->Rule,And->List};
		localL/@Range[m] /. soln
	];
AssignLabelsNoMomCons[args_] := 
	Module[
		{g=args,
			edges,v,e,localL,n,m,i,
			externalV,internalV,externalE,internalE,
			momcons,pos,sign,soln},
		m = Max@Flatten[g];
		edges = Table[{Position[g,-i][[1,1]],Position[g,i][[1,1]]},{i,m}];
		e = Length[edges];
		v = Length/@g;
		externalV = Flatten@Position[v,1]; (* External vertices *)
		externalE = Flatten@Abs[g[[externalV]]];
		Table[If[MemberQ[externalE,i],p[i],l[i]],{i,1,e}]
	];


ClearAll[GenerateJacobis];
GenerateJacobis[n_,L_,minGirth_] := 
	Module[
		{i1=1,i2,i3,triples = {{}},i=0,j, 
			me,permutedMe,mySignature,mySignatures,
			you,permutedYou,yourIndex,isoSignatures,
			edgeAutomorphisms,labelAutomorphisms,
			allYourMomenta,yourTargetLabels,yourReplacedLabels,momentumRelabelings,eqns,solns,\[Beta],
			eaut,restrictedEAut,toCheck,baseInvariants,myInvariant,
			selfAntisymmetric,positions,yourEdge,
			uniqueEdges,alreadyChecked,iso,restrictedIso,testCondition},
	(*****************************************************************************************)
	(* Build the starting diagram. *)
	(*****************************************************************************************)
		BASISJACOBIS = {};
		If[n==3&&L==1,
			INTEGRALBASIS = BuildGraphGJ/@{{{-1},{-2},{-3},{-5,1,4},{-4,2,6},{-6,3,5}}};
			MINGIRTH = Max[2,minGirth];
			,
			If[n==3 && L==2,
				INTEGRALBASIS = BuildGraphGJ/@{{{-1},{-2},{-3},{1,4,5},{-4,6,7},{-5,8,9},{-6,2,-8},{-7,3,-9}}};
				MINGIRTH = Max[2,minGirth];,
				Which[
					L==0,
						INTEGRALBASIS = BuildGraphGJ/@{Join[Table[{-i},{i,1,n}],
							{{1,2,n+1}},
						Table[{-n-i,i+2,n+i+1},{i,1,n-4}],
								{{-2n+3,n-1,n}}]};
						MINGIRTH = 0;,
					L>=1,
						INTEGRALBASIS = BuildGraphGJ/@{Join[Table[{-i},{i,1,n}],
							{{n+2L+1,2,n+1}},
							{{n+L+1,1,-n-2L-1}},
							Table[{n+2L+j,-n-j+1,n+j},{j,2,L}],
							Table[{n+L+j,-n-L-j+1,-n-2L-j},{j,2,L}],
							{{-n-L,3,n+3L+1}},
							Table[{-n-3L-j,j+3,n+3L+j+1},{j,1,n-4}],
							{{-n-3L-(n-3),n,-n-2L}}]}; 
						MINGIRTH = Max[2,minGirth];,
					True,
						Break[];
				];
			];
		];
		PLANAR = {True};
		edgeAutomorphisms = {EdgewiseAutomorphisms[INTEGRALBASIS[[1]]]};
		labelAutomorphisms = {LabelwiseAutomorphisms[INTEGRALBASIS[[1]]]};
		alreadyChecked = {{}};
		baseInvariants = {PreIsoInvariant[INTEGRALBASIS[[1]]]};
		PrintTemporary["Working on graph ",Dynamic[i1]," of ",Dynamic[Length[INTEGRALBASIS]],"."];
		For[i1=1,i1<=Length[INTEGRALBASIS],i1++,
	(*****************************************************************************************)
	(* Set up initial variables *)
	(*****************************************************************************************)
			me = INTEGRALBASIS[[i1]];
	(*****************************************************************************************)
	(* Generate a Jacobi identity on each distinct internal leg, *) 
	(* unless it has already been involved in another Jacobi *)
	(* This avoids writing every Jacobi three times. *)
	(*****************************************************************************************)
			{mySignatures,eaut} = edgeAutomorphisms[[i1]];
			uniqueEdges = Sort/@(Abs[eaut]\[Transpose][[n+1;;-1]]) //DeleteDuplicates; 
			uniqueEdges = uniqueEdges[[All,1]];
			uniqueEdges = Complement[uniqueEdges,alreadyChecked[[i1]]];
			triples = Table[Jacobi[me,i],{i,uniqueEdges}]; 
			positions = triples; (* Initialize array to the correct size *)
	(*****************************************************************************************)
	(* Get momentum relabelings and assemble RG structure *)
	(* Note: To do this, "me" has to be unaltered by Jacobi. *)
	(* Otherwise, "me" would have to be treated like the other two graphs here. *)
	(*****************************************************************************************)
			momentumRelabelings = labelAutomorphisms[[i1,2]];
			positions[[All,1]] = Table[RG[bg[i1],mySignatures,momentumRelabelings],{Length[triples]}]; 
			For[i=1,i<=Length[triples],i++,
				For[j=2,j<=Length[triples[[i]]],j++,
	(*****************************************************************************************)
	(* Task: Identify each new graph with a base graph in INTEGRALBASIS, *)
	(* or append the graph to INTEGRALBASIS, *)
	(*****************************************************************************************)
					me = triples[[i,j]];
					yourIndex = 0;
	(*****************************************************************************************)
	(* Check equality of a graph invariant before wasting time with Isomorphism. *)
	(*****************************************************************************************)
					myInvariant = PreIsoInvariant[me];
					toCheck = Flatten@Position[baseInvariants,myInvariant,{1},Heads->False];
					For[i2=1,i2<=Length[toCheck],i2++,
						{isoSignatures,iso} = EdgewiseIsomorphism[INTEGRALBASIS[[ toCheck[[i2]] ]],me];
						If[Length[iso]>0,
							yourIndex = toCheck[[i2]];
	(* Need to determine which edge was reached. It is edge uniqueEdges[[i]] in "me". What edge is it in the base graph? *)
							momentumRelabelings = LabelwiseIsomorphism[INTEGRALBASIS[[yourIndex]],me,#]&/@iso;
							yourEdge = Position[Abs[iso[[1]]],uniqueEdges[[i]],{1},Heads->False ][[1,1]];
							yourEdge = edgeAutomorphisms[[yourIndex,2,All,yourEdge]]; 
							yourEdge = DeleteDuplicates[Abs[yourEdge]];
							alreadyChecked[[yourIndex]] = Union[alreadyChecked[[yourIndex]],yourEdge];
							Break[];
						];
					]; 
	(*****************************************************************************************)
	(* If the new graph is not in the INTEGRALBASIS, add the graph to the list and store its properties. *)
	(*****************************************************************************************)
					If[yourIndex == 0,
						INTEGRALBASIS = Append[INTEGRALBASIS,BuildGraphGJ[me[[1]]]]; 
						yourIndex = Length[INTEGRALBASIS];
						{isoSignatures,iso} = EdgewiseAutomorphisms[INTEGRALBASIS[[yourIndex]]];
						AppendTo[edgeAutomorphisms,{isoSignatures,iso}];
						AppendTo[labelAutomorphisms,
							{
								edgeAutomorphisms[[yourIndex,1]],
								LabelwiseIsomorphism[
									INTEGRALBASIS[[yourIndex]],
									INTEGRALBASIS[[yourIndex]],#]&/@edgeAutomorphisms[[yourIndex,2]
								]
							}
						];
						yourEdge = uniqueEdges[[i]];
						yourEdge = edgeAutomorphisms[[yourIndex,2,All,yourEdge]]; 
						yourEdge = DeleteDuplicates[Abs[yourEdge]];
						AppendTo[alreadyChecked,yourEdge];
						AppendTo[PLANAR, PlanarGraphQQ[ INTEGRALBASIS[[yourIndex]] ] ];
						AppendTo[baseInvariants,myInvariant];
						isoSignatures = {1};
						momentumRelabelings = 
							{LabelwiseIsomorphism[
								INTEGRALBASIS[[ yourIndex ]], 
								me, 
								Range[ Length[ me[[2]] ] ] 
							]};
					];	
	(*****************************************************************************************)
	(* Assemble RG structure. *)
	(*****************************************************************************************)
					isoSignatures = labelAutomorphisms[[yourIndex,1]] * isoSignatures[[1]];
					momentumRelabelings = labelAutomorphisms[[yourIndex,2]] /. Rule[a_,b_]:>Rule[a,b/.momentumRelabelings[[1]]];
					positions[[i,j]] = SortRG[RG[bg[yourIndex],isoSignatures,momentumRelabelings]];
				];
			];
			BASISJACOBIS = Join[BASISJACOBIS,positions];
		];
		i1 = Length[INTEGRALBASIS];
		BASISJACOBIS = Total/@BASISJACOBIS;
		DEPENDENCIES = Range[Length[INTEGRALBASIS]];
		(*Print["Data saved in global variables INTEGRALBASIS, BASISJACOBIS, DEPENDENCIES and PLANAR."];*)
	];


(* ::Subsection::Closed:: *)
(*User Interface: GetGraphBasis*)


(*ClearAll[GetGraphBasis];*)
GetGraphBasis[n_Integer, L_Integer, minGirth_Integer:0] := Module[
{directoryGB,directoryGJ,name,suffix,fileNameGB,fileNameGJ,gb,gj},
	Clear[INTEGRALBASIS,BASISJACOBIS];
	(*Search if saved*)
	directoryGB = MASTERFUNCTIONSDIRECTORY<>"aux/graphBasis/";
	directoryGJ = MASTERFUNCTIONSDIRECTORY<>"aux/graphJacobis/";
	name = ToString[n]<>"."<>ToString[L]<>"."<>ToString[minGirth];
	suffix = ".m";
	fileNameGB = directoryGB<>name<>suffix;
	fileNameGJ = directoryGJ<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading graph basis..."];
	PrintQ[AbsoluteTiming[Quiet[gb = Get[fileNameGB];];][[1]]];
	If[ gb === $Failed,
		PrintQ["Can't find saved graph basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating graphs and Jacobi relations..."];
		PrintQ[AbsoluteTiming[GenerateJacobis[n,L,minGirth];][[1]]];
		gb = INTEGRALBASIS[[All,1]];
		gj = BASISJACOBIS;
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameGB];];
		Put[gb,fileNameGB];
		Quiet[DeleteFile[fileNameGJ];];
		Put[gj,fileNameGJ];
	];
	Clear[INTEGRALBASIS,BASISJACOBIS];
	Return[gb];
];


(* ::Subsection::Closed:: *)
(*User Interface: GetGraphJacobis*)


(*ClearAll[GetGraphJacobis];*)
GetGraphJacobis[n_Integer, L_Integer, minGirth_Integer:0] := Module[
{directoryGB,directoryGJ,name,suffix,fileNameGB,fileNameGJ,gb,gj},
	Clear[INTEGRALBASIS,BASISJACOBIS];
	(*Search if saved*)
	directoryGB = MASTERFUNCTIONSDIRECTORY<>"aux/graphBasis/";
	directoryGJ = MASTERFUNCTIONSDIRECTORY<>"aux/graphJacobis/";
	name = ToString[n]<>"."<>ToString[L]<>"."<>ToString[minGirth];
	suffix = ".m";
	fileNameGB = directoryGB<>name<>suffix;
	fileNameGJ = directoryGJ<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading graph basis..."];
	PrintQ[AbsoluteTiming[Quiet[gj = Get[fileNameGJ];];][[1]]];
	If[ gj === $Failed,
		PrintQ["Can't find saved graph basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating graphs and Jacobi relations..."];
		PrintQ[AbsoluteTiming[GenerateJacobis[n,L,minGirth];][[1]]];
		gb = INTEGRALBASIS[[All,1]];
		gj = BASISJACOBIS;
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameGB];];
		Put[gb,fileNameGB];
		Quiet[DeleteFile[fileNameGJ];];
		Put[gj,fileNameGJ];
	];
	Clear[INTEGRALBASIS,BASISJACOBIS];
	Return[gj];
];


(* ::Section::Closed:: *)
(*fabc Color Basis*)


(* Overkill for most purposes, but it finds everything you could want to know about color. *)
(* I wouldn't go beyond n = 8 though *)


(* 
There are more efficient algorithms to generate "color-ordered basis,"
i.e. orderings that contribute to color-ordered amplitudes.
See http://arxiv.org/pdf/1506.00974.pdf for instance
*)


(* Note: Had to take one ParallelTable out of GenerateSymBasis, otherwise it forgets what perms is. *)


(* ::Subsection::Closed:: *)
(*Generate Symmetry Basis*)


ClearAll[GenerateSymBasis];
GenerateSymBasis[n_Integer] := Module[
{symBasis,symReps,ib,i,graph,ea,ordEC,signs,perms,permSigns,tab,pos,reps,sig},
	symBasis = {};
	symReps = {};
	ib = GetGraphBasis[n,0];
	For[ i = 1, i <= Length[ib], i++,
		graph = ib[[i]];
		(*Get automorphisms of graph i*)
		ea = (Transpose[EdgewiseAutomorphisms[graph]]);
		(* Ordering equivalence class *)
		ordEC = ea[[All,2]][[All,1;;n]];
		signs = ea[[All,1]];
		perms = Permutations[Range[n]];
		(* Gather permutations by symmetry equivalence class, then take first representatives *)
		perms = GatherBy[
					perms,
					Sort[
						Table[#[[ordEC[[i]]]],{i,1,Length[ordEC]}]
					][[1]]&
				];
		permSigns = ParallelTable[
						signs[[Flatten[ 
									Position[
										Table[perms[[j,1]][[ordEC[[k]]]],{k,1,Length[ordEC]}],#
									]&/@(perms[[j]]) 
								]
							]],
					{j,1,Length[perms]}
					,DistributedContexts :> distributedContexts
					];
		(* Take the first representative from each equivalence class, and add it to the symBasis *)
		AppendTo[ symBasis, fColor[i,#]&/@(perms[[All,1]]) ];
		(* Make replacement rules to put into symmetry basis *)
		reps = (*Parallel*)Table[
					Table[ 
						fColor[i,perms[[j,k]]] -> ( permSigns[[j,1]] * permSigns[[j,k]] * fColor[i,perms[[j,1]]] ),
						{k,1,Length[perms[[j]]]}
					],
				{j,1,Length[perms]}
				(*,DistributedContexts :> distributedContexts*)
				];
		AppendTo[ symReps, reps ];
	];
	symBasis = Flatten[symBasis];
	symReps = Flatten[symReps];

	{symBasis,symReps}
];


(*ClearAll[GetSymBasis];*)
GetSymBasis[n_Integer] := Module[
{directorySB,directorySR,name,suffix,fileNameSB,fileNameSR,
symBasis,symReps},
	(*Search if saved*)
	directorySB = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/symBasis/";
	directorySR = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/symReps/";
	name = ToString[n];
	suffix = ".m";
	fileNameSB = directorySB<>name<>suffix;
	fileNameSR = directorySR<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading color symmetry basis..."];
	PrintQ[AbsoluteTiming[Quiet[symBasis = Get[fileNameSB];];][[1]]];
	If[ symBasis === $Failed,
		PrintQ["Can't find saved color symmetry basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating color symmetry basis and replacements..."];
		{symBasis,symReps} = GenerateSymBasis[n];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameSB];];
		(*Put[symBasis,fileNameSB];*)
		Quiet[DeleteFile[fileNameSR];];
		(*Put[symReps,fileNameSR];*)
	];

	symBasis
];


(*ClearAll[GetSymReps];*)
GetSymReps[n_Integer] := Module[
{directorySB,directorySR,name,suffix,fileNameSB,fileNameSR,
symBasis,symReps},
	(*Search if saved*)
	directorySB = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/symBasis/";
	directorySR = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/symReps/";
	name = ToString[n];
	suffix = ".m";
	fileNameSB = directorySB<>name<>suffix;
	fileNameSR = directorySR<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading color symmetry replacement rules..."];
	PrintQ[AbsoluteTiming[Quiet[symReps = Get[fileNameSR];];][[1]]];
	If[ symReps === $Failed,
		PrintQ["Can't find saved color symmetry basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating color symmetry basis and replacements..."];
		{symBasis,symReps} = GenerateSymBasis[n];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameSB];];
		Put[symBasis,fileNameSB];
		Quiet[DeleteFile[fileNameSR];];
		Put[symReps,fileNameSR];
	];

	symReps
];


(* ::Subsection::Closed:: *)
(*Generate DDM Basis Using Jacobis*)


(* DDM Basis: *)
(* See Eq.(2.1) and Fig. 1 of http://arxiv.org/pdf/hep-ph/9910563.pdf *)
(* Final basis is fColor[1, {1, ...<(n-2)! perms>..., n}] *)


ClearAll[GenerateDDMBasis];
GenerateDDMBasis[n_Integer] := Module[
{bj,symBasis,symReps,eqs,j,jacob,perms,ddmBasis,aRepBase,vars,aRepOther,msReps,msInvReps,soln,ib,ddmReps},
	bj = GetGraphJacobis[n,0];
	symBasis = GetSymBasis[n];
	symReps = GetSymReps[n];
	perms = Permutations[Range[n]];
	perms = Thread[Table[a[i],{i,1,n}] -> #]&/@perms;
	(* We look at all permutations of labels, then delete duplicate equations *)
	(* Is there a better way? My first attempt just looking at symBasis orderings failed *)
	(* Still, this GatherBy is relatively fast. The bottleneck appears to be in generating the symBasis/Reps *)
	eqs = ParallelTable[
			jacob = Table[{bj[[j,k,1]]/.bg[x_]:>x,bj[[j,k,2,1]],bj[[j,k,3,1]]},{k,1,Length[bj[[j]]]}];
			jacob = Plus@@(#[[2]]*fColor[#[[1]],(Table[\[Epsilon][m],{m,1,n}]/.#[[3]])/.\[Epsilon]->a]&/@jacob);
			jacob = (jacob /. #)&/@perms;
			jacob = Expand[jacob /. Dispatch[symReps]];
			jacob = GatherBy[jacob,Sort[{#,Expand[-#]}]&];
			jacob = jacob[[All,1]];
			jacob,
		{j,1,Length[bj]},
		DistributedContexts -> distributedContexts
		];
	eqs = Expand[Flatten[eqs]];

	(* Solve in terms of the standard color factors *)
	(* See Eq.(2.1) or Fig. 1 of http://arxiv.org/pdf/hep-ph/9910563.pdf *)
	(* GenerateJacobis makes the first diagram the ladder diagram of Fig. 1, so this is kosher *)
	(* The following is formatting to make MonsterSolve happy *)
	ddmBasis = fColor[1,Join[{1},#,{n}]]&/@Permutations[Table[i,{i,2,n-1}]];
	aRepBase = Thread[ddmBasis -> Table[a[1,i],{i,1,Length[ddmBasis]}]];
	vars = GetVars[eqs,ddmBasis];
	aRepOther = Table[vars[[i]] -> a[0,i],{i,1,Length[vars]}];
	msReps = Join[aRepBase,aRepOther];
	msInvReps = msReps /. Rule[a_,b_]:>Rule[b,a];
	Print["Solving "<>ToString[Length[eqs]]<>" equations..."];
	soln = MonsterSolve2[eqs/.Dispatch[msReps],100,20]/. Dispatch[msInvReps];

	ib = GetGraphBasis[n,0];
	(* All color factors *)
	ddmReps = Flatten[ Table[ fColor[i,#],{i,1,Length[ib]} ]&/@Permutations[Range[n]] ];
	(* All color replacement rules under symmetry and Jacobis *)
	ddmReps = Thread[ ddmReps -> (( ddmReps /. Dispatch[symReps] ) /. Dispatch[soln] ) ];

	{ddmBasis,ddmReps}
];


(*ClearAll[GetDDMBasis];*)
GetDDMBasis[n_Integer] := Module[
{directoryB,directoryR,name,suffix,fileNameB,fileNameR,
ddmBasis,ddmReps},
	(*Search if saved*)
	directoryB = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/DDMBasis/";
	directoryR = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/DDMReps/";
	name = ToString[n];
	suffix = ".m";
	fileNameB = directoryB<>name<>suffix;
	fileNameR = directoryR<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading DDM basis..."];
	PrintQ[AbsoluteTiming[Quiet[ddmBasis = Get[fileNameB];];][[1]]];
	If[ ddmBasis === $Failed,
		PrintQ["Can't find saved DDM basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating DDM basis and replacements..."];
		{ddmBasis,ddmReps} = GenerateDDMBasis[n];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameB];];
		Put[ddmBasis,fileNameB];
		Quiet[DeleteFile[fileNameR];];
		Put[ddmReps,fileNameR];
	];

	ddmBasis
];


(*ClearAll[GetDDMReps];*)
GetDDMReps[n_Integer] := Module[
{directoryB,directoryR,name,suffix,fileNameB,fileNameR,
ddmBasis,ddmReps},
	(*Search if saved*)
	directoryB = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/DDMBasis/";
	directoryR = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/DDMReps/";
	name = ToString[n];
	suffix = ".m";
	fileNameB = directoryB<>name<>suffix;
	fileNameR = directoryR<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading DDM replacement rules..."];
	PrintQ[AbsoluteTiming[Quiet[ddmReps = Get[fileNameR];];][[1]]];
	If[ ddmReps === $Failed,
		PrintQ["Can't find saved DDM basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating DDM basis and replacements..."];
		{ddmBasis,ddmReps} = GenerateDDMBasis[n];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileNameB];];
		Put[ddmBasis,fileNameB];
		Quiet[DeleteFile[fileNameR];];
		Put[ddmReps,fileNameR];
	];

	ddmReps
];


(* ::Subsection::Closed:: *)
(*Generate Color-Ordered Basis for (12 ... n) Ordering*)
(*(Maybe "Basis" is a bad word here)*)


(* Find graph orderings that contribute to fColor[1, {1,2,...,n}] *)


ClearAll[GenerateCOBasis];
GenerateCOBasis[n_Integer] := Module[
{symBasis,ddmBasis,ddmReps,mask,maskRep,signs,coBasis},
	symBasis = GetSymBasis[n];
	ddmBasis = GetDDMBasis[n];
	ddmReps = GetDDMReps[n];
	(*Get reps that change f[1,{permutation}] into contribution to f[1,{1,2,...,n-1,n}]*)
	mask = ConstantArray[0,Length[ddmBasis]];
	mask[[1]] = 1;
	maskRep = Thread[ddmBasis->mask];
	signs = ( symBasis /. Dispatch[ddmReps] ) /. Dispatch[maskRep];
	coBasis = List@@(Plus@@(signs * symBasis));

	coBasis
];


(*ClearAll[GetCOBasis];*)
GetCOBasis[n_Integer] := Module[
{directory,name,suffix,fileName,
coBasis},
	(*Search if saved*)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/colorTree/COBasis/";
	name = ToString[n];
	suffix = ".m";
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading color-ordered basis..."];
	PrintQ[AbsoluteTiming[Quiet[coBasis = Get[fileName];];][[1]]];
	If[ coBasis === $Failed,
		PrintQ["Can't find saved color-ordered basis."];
		(*If not saved build it and save it*)
		PrintQ["Generating color-ordered basis..."];
		coBasis = GenerateCOBasis[n];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[coBasis,fileName];
	];

	coBasis
];


(* ::Section::Closed:: *)
(*Ansatz Builder*)


(* This is old and probably pretty bad *)
(* It's messy and unreadable at least *)
(* ...but it seems to work *)
(* Should rewrite *)
(* Can we immediately incorporate symmetry constraints? *)


(* ::Subsection::Closed:: *)
(*BuildAnsatz*)


(*Clear[BuildAnsatz];*)
BuildAnsatz[verts_,coeff_,N_:0] := Module[
{n,legChoices,maxPowers,numVert,internalVertices,externalEdges,internalEdges,
epList,eps,f,epTerms,epTotMom,epLoopMom,epsGoodAlonePos,epsNoGoodPos,epsGoodAlone,momList,moms,g,ansTerms,sum,
momrules,propsList},
	
	{legChoices,maxPowers}=ChooseLoopMom[verts];
	{externalEdges,internalEdges,internalVertices} = GraphInfo[verts];
	n=Length[externalEdges];
	numVert=Length[internalVertices];
    
	(*SUSY reduces your maximum powers of loop momentua*)
	maxPowers=maxPowers-Switch[N,0,0,1,2,2,2,4,4,_,0];

	(**********Build \[Epsilon] terms*********************)
	epList=Join[Flatten@Table[If[i<j,{ee[i,j]},{}], {i,1,n}, {j,1,n}], Flatten@Table[ke[i,j], {i,1,n}, {j,1,n}], Flatten@Table[le[i,j], {i,legChoices}, {j,1,n}]];
	eps=Variables[Map[MomConsN[n,#]&,epList]];

	f[Range[n],x___] := {{x}}; (*Get all n \[Epsilon]'s*)
	f[{___,x_,___,x_,___},___] ={}; (*Don't duplicate \[Epsilon]'s*)
	(*Stops when we get all n or when we try to apply on an empty list*) 
	f[{x___},a___,b_]:= Join@@(f[
		Sort[Flatten[Join[{x},#/.{ee[e1_,e2_]:>{e1,e2},ke[k_,e_]:>{e},le[l_,e_]:>{e}}]]],
		a ,b,#]&
		/@ eps[[Position[eps,b][[1]][[1]]+1;;-1]]);

	epTerms=Join@@Flatten[f[#/.{ee[e1_,e2_]:>{e1,e2},ke[k_,e_]:>{e},le[l_,e_]:>{e}},#]&/@eps,{2}];
	
	epTotMom=Exponent[Times@@(#/.{ee[e1_,e2_]:>1,ke[k_,e_]:>\[Alpha],le[l_,e_]:>\[Alpha]}),\[Alpha]]&/@epTerms;
	epLoopMom=Table[Count[#,legChoices[[i]],\[Infinity]],{i,1,Length[legChoices]}]&/@epTerms;
	(*********************************************)
	
	(****Pick out the complete terms and drop the disqualified terms****)
	epsNoGoodPos=Join@@Table[
					If[Or@@
						Join[Table[epLoopMom[[j]][[i]]>maxPowers[[i]],{i,1,Length[maxPowers]}],
						{epTotMom[[j]]>numVert}],
						{j},
						{}
					],
				{j,1,Length[epTerms]}];

	epsGoodAlonePos=Join@@Position[epTotMom,numVert];

	epsGoodAlone=epTerms[[Complement[epsGoodAlonePos,epsNoGoodPos]]];
	epTerms = Delete[epTerms, Map[{#}&,Union[epsNoGoodPos,epsGoodAlonePos]]];
	epTotMom = Delete[epTotMom, Map[{#}&,Union[epsNoGoodPos,epsGoodAlonePos]]];
	epLoopMom = Delete[epLoopMom, Map[{#}&,Union[epsNoGoodPos,epsGoodAlonePos]]];
	(*********************************************************************)

	(************Fill out rest of momentum*******)
	momList=Join[
					Flatten@Table[If[i<=j,{ll[i,j]},{}], {i,legChoices}, {j,legChoices}],  
					Flatten@Table[lk[i,j], {i,legChoices}, {j,1,n}],
					Flatten@Table[If[i<=j,{kk[i,j]},{}], {i,1,n}, {j,1,n}]
				];
	moms=Variables[Map[MomConsN[n,#]&,momList]];

	g[{l___},{pow__},x__]:={}/;Or@@Table[{l}[[i]]>maxPowers[[i]],{i,1,Length[maxPowers]}];
	g[{l___},{pow__},x__]:={}/;pow>numVert;
	g[{l___},{pow__},x__]:={{x}}/;pow===numVert;
	g[{l___},{pow__},x__]:=Join@@(g[
		Table[{l}[[i]]+Count[List@@#,legChoices[[i]]],{i,1,Length[legChoices]}],
		{pow+2},
		x,#]&
		/@ moms[[Position[moms,{x}[[-1]]][[1]][[1]];;-1]]);

	ansTerms=Join@@Table[
		Map[Join[epTerms[[j]],#1]&,
			Join@@Flatten[
				g[
					Table[Count[List@@#,legChoices[[i]]]+epLoopMom[[j,i]],{i,1,Length[legChoices]}],
					{epTotMom[[j]]+2},
					#
				]
			&/@ moms,{2}]
		],
	{j,1,Length[epTerms]}];

	ansTerms = Join[epsGoodAlone,ansTerms];

	sum=Plus@@Table[coeff[i]*Times@@(ansTerms[[i]]),{i,1,Length[ansTerms]}];

	(*Build props*)
	momrules = MomRules[verts];
	propsList = R/@((p/@internalEdges)/.momrules);

	{sum,propsList}
];


(* ::Section::Closed:: *)
(*YM Feynman Rules*)


(* Modifications to code originally by Tristan Dennen *)


(* Scalar Feynman rules haven't been tested yet! *)
(* Theories involving explict scalar computation are commented out in AssignFields *)
(* That being said we can get N = 4 in D = 4 from N = 1 with Ds = 10, Df = 8 *)


(* ::Subsection::Closed:: *)
(*Edges Dictionary*)


(* Dictionary between edges to a sequential labeling *)
ClearAll[dictEdges]
dictEdges[vertices_List] := Module[
{externalEdges,internalEdges,edges},
	{externalEdges,internalEdges} = GraphInfo[vertices][[1;;2]];
	edges = Join[externalEdges,internalEdges];
	Thread[edges->Range[Length[edges]]]
];


(* ::Subsection::Closed:: *)
(*A (Feynman vertex definitions)*)


(*See Peskin and Schroeder Appendix*)
ClearAll[A];
A[args__] := 
	Module[
		{vertexType = {args}[[1]],
			index = {args}[[2]],
			mom = {args}[[3]],
			amplitude,
			ord},
		(* - indicates incoming ?? Momenta conventionally outgoing. {_,_,_} is read clockwise *)
		vertexType = vertexType //. {-"A"->"A", -"B"->"B", -"\[Phi]"->"-\[Phi]", -"\[Psi]"->"-\[Psi]", -"c"->"-c", -"\[Chi]"->"-\[Chi]"};
		ord = Ordering[vertexType];
		vertexType = vertexType[[ord]];
		index = index[[ord]];
		mom = mom[[ord]];
	
		Switch[vertexType,
			{"A","A","A"},
				amplitude = 1/(I Sqrt[2]) * (
					(mom[[1]] - mom[[2]])[ \[Mu][ index[[3]] ] ] \[Eta][\[Mu][index[[1]]],\[Mu][index[[2]]]] +
					(mom[[2]] - mom[[3]])[ \[Mu][ index[[1]] ] ] \[Eta][\[Mu][index[[2]]],\[Mu][index[[3]]]] + 
					(mom[[3]] - mom[[1]])[ \[Mu][ index[[2]] ] ] \[Eta][\[Mu][index[[3]]],\[Mu][index[[1]]]] );,
			{"A","A","B"},
				amplitude =  1/(I 2)  (
					\[Eta][ \[Mu][index[[1]]], \[Mu][index[[3]],1] ] \[Eta][ \[Mu][index[[2]]], \[Mu][index[[3]],2] ] - 
					\[Eta][ \[Mu][index[[2]]], \[Mu][index[[3]],1] ] \[Eta][ \[Mu][index[[1]]], \[Mu][index[[3]],2] ]);,
			{"A","-c","c"}, 
				amplitude =  1/(I Sqrt[2]) (mom[[3]])[ \[Mu][ index[[1]] ] ];,
			(*Contract with the incoming arrow. The Subscript[matrix [\[Gamma]^\[Mu]]^out, in] *)
			{"A","-\[Psi]","\[Psi]"}, (* Put I in denom because it works *)
				amplitude = 1/(I Sqrt[2]) \[CapitalGamma][{\[Mu][index[[1]]]}, {\[Mu][index[[3]],2], \[Mu][index[[2]],1]}];, (* 9/11/2012 - Josh Nohle *)
			{"A","-\[Phi]","\[Phi]"},
				amplitude = 1/(I Sqrt[2]) (mom[[2]] + mom[[3]])[\[Mu][index[[1]]]];
		];
		amplitude * Signature[ord]
	]; 


(* ::Subsection::Closed:: *)
(*Propagators*)


(*See Peskin and Schroeder Appendix*)
(* Added the edge arguement necessary for the \[CapitalGamma]s - 9/13/2012 - Josh Nohle *)
ClearAll[Propagator];
Propagator[label_, field_, edge_] := Module[
{},	
	Switch[field,
		"A", -I/R[label],
		"B", -I*S[label]/R[label],
		-"c", I/R[label],
		"c", I/R[label],
		-"\[Psi]", (I/R[label])*\[CapitalGamma][{\[Mu][edge,3]},{\[Mu][edge,1],\[Mu][edge,2]}]*(-label[\[Mu][edge,3]]),(* 9/12/2012 - Josh Nohle *)
		"\[Psi]", (I/R[label])*\[CapitalGamma][{\[Mu][edge,3]},{\[Mu][edge,1],\[Mu][edge,2]}]*(label[\[Mu][edge,3]]),(* 9/12/2012 - Josh Nohle *)
		-"\[Phi]", I/R[label],
		"\[Phi]", I/R[label]		
	]
];
ClearAll[Propagators];
Propagators[vertices_List, fieldAssignments_List] := Module[
{internalEdges, dict},
	internalEdges = GraphInfo[vertices][[2]];
	dict = dictEdges[ vertices ];
	Times@@Map[ Propagator[p[#], fieldAssignments[[#/.dict]], #]&, internalEdges ]
];


(* ::Subsection::Closed:: *)
(*ValidVertexQ (Allowed vertices)*)


ClearAll[ValidVertexQ];
ValidVertexQ[arg_] := 
	Module[
		{fieldAssignments = arg //. {-a_ :> "-" <> a},
			validVertices = {"A"^3, "A"^2 "B", "A" "-c" "c", "A" "-\[Psi]" "\[Psi]", "A" "-\[Phi]" "\[Phi]"} },
		fieldAssignments = fieldAssignments /. {"-A"->"A","-B"->"B"};
		fieldAssignments = Times@@fieldAssignments;
		Or@@( (#===1)&/@Denominator[validVertices/fieldAssignments] )
	];


(* ::Subsection::Closed:: *)
(*CountGhostLoops*)


ClearAll[CountGhostLoops];
CountGhostLoops[vertices_List, fieldAssignments_List] := Module[
{edges, ghostLegs, cc},
		edges = BuildGraph[vertices][[2]];
		ghostLegs = Position[fieldAssignments,"c"][[All,1]];
		edges = edges[[ghostLegs]];
		cc = Combinatorica`ConnectedComponents@Combinatorica`FromOrderedPairs[edges];
		DeleteCases[cc, a_/;Length[a]==1] //Length
	];


(* ::Subsection::Closed:: *)
(*CountFermionLoops*)


ClearAll[CountFermionLoops];
CountFermionLoops[vertices_List, fieldAssignments_List] := Module[
{edges, fermLegs, cc},
		edges = BuildGraph[vertices][[2]];
		fermLegs = Position[fieldAssignments,"\[Psi]"][[All,1]];
		edges = edges[[fermLegs]];
		cc = Combinatorica`ConnectedComponents@Combinatorica`FromOrderedPairs[edges];
		DeleteCases[cc,a_/;Length[a]==1] //Length
	];


(* ::Subsection::Closed:: *)
(*AssignFields*)


ClearAll[AssignFields];
AssignFields[vertices_List, N_] := Module[
{externalEdges, internalEdges, edges, dict, fields, fieldAssignments, i1, i2, toCheck},
	
	(* Set the field of each external propagator to "A" *)
	{externalEdges,internalEdges} = GraphInfo[vertices][[1;;2]];
	edges = Join[externalEdges,internalEdges];
	(* Map edges to 1,2,3,...,number of edges *)
	(* Don't assume edges have sequential numbering *)
	dict = dictEdges[ vertices ];
	(* a minus sign indicates particle is in the opposite sense to momentum direction *)		
	fields = Switch[N,
		-1, {"A","c",-"c"}, (*Just used to find all of the loops in a graph*)
		0, {"A","B","c",-"c"},
		1, {"A","B","c",-"c","\[Psi]",-"\[Psi]"}, (*Vector multiplet*)
		(*2, {"A","B","c",-"c","\[Psi]",-"\[Psi]","\[Phi]"},
		4, {"A","B","c",-"c","\[Psi]",-"\[Psi]","\[Phi]"},*)
		_, Message["Not a valid/supported theory. Defaulting to N=0."];{"A","B","c",-"c"}
	
	];
	fieldAssignments = {Table[1,{Length[edges]}]};	
	fieldAssignments[[1,externalEdges/.dict]] = "A";

	(* One by one, add field to each internal edge and check for validity *)
	For[i1 = 1, i1 <= Length[internalEdges], i1++,
		fieldAssignments = Reap[
			For[i2 = 1, i2 <= Length[fields], i2++,
				Sow[ReplacePart[#,(internalEdges[[i1]]/.dict)->fields[[i2]] ]&/@fieldAssignments];
			];
		][[2,1]];
		fieldAssignments = Flatten[fieldAssignments,1];
		(* Check for validity *)
		toCheck = {Position[vertices,internalEdges[[i1]]][[1,1]],Position[vertices,-internalEdges[[i1]]][[1,1]]};
		toCheck = vertices[[ toCheck ]];(* Vertices to check *)
		fieldAssignments = fieldAssignments[[Flatten@Position[
			Table[And@@ValidVertexQ/@(Sign[#]*fieldAssignments[[i3,Abs[#]/.dict]]&/@toCheck),
				{i3,Length[fieldAssignments]}],
			True
		]]];
	];
	fieldAssignments
];


(* ::Subsection::Closed:: *)
(*GammaTrace*)


(* 
This is old and probably bad. 
Tensor reduction is similar but my code is much faster there, so we should probably use that.
*)


ClearAll[GammaTrace,SplitByTwosG];
SplitByTwosG[arg_List]:=  
	Catch[Block[
		{list=arg,rlist = Rest[arg],pairs,set},
		If[Length[list]==2,Throw[{{list}}];];
		pairs = Table[Join[{{list[[1]],list[[i]]}},#]&/@SplitByTwosG[Drop[rlist,{i-1}]],{i,2,Length[list]}];
		Throw[Flatten[pairs,1]];
	set
	]];

GammaTrace[arg__]:=
	Module[{list,set,ans,temp},
		list=arg//.{a___,\[Mu][b__],c___,\[Mu][b__],d___}:>{a,\[Mu][b],c,\[Mu][b,temp],d};
		set = SplitByTwosG[list];
		(*In D dimensions, the dimension of the Clifford Algebra representation is 2^\[LeftFloor]D/2\[RightFloor]...sorta*)
		ans=Plus@@Table[
			(*Df are the number of fermion degrees of freedom*)
			Df*Signature[list]*Signature[Flatten[set[[i]]]]*Times@@Apply[\[Eta],set[[i]],1],{i,1,Length[set]}];
		ans//.\[Mu][a__,temp]:>\[Mu][a]
	];


(* ::Subsection::Closed:: *)
(*FeynmanYM*)


(* Given a cubic graph, calculate from Feynman rules, including all valid propagating fields *)
(*ClearAll[FeynmanYM];*)
FeynmanYM[vertices_List, N_:0] := Module[
{internalVertices,externalEdges,internalEdges,dict,
fieldAssignments,ghostLoops,fermLoops,propagators,amplitude},

	{externalEdges,internalEdges,internalVertices} = GraphInfo[vertices];
	dict = dictEdges[ vertices ];

	(* Still need to deal with different number of fields in theories of different N.*)
	fieldAssignments = AssignFields[vertices, N];
	ghostLoops = Table[CountGhostLoops[vertices, fa], {fa,fieldAssignments}];
	fermLoops = Table[CountFermionLoops[vertices, fa], {fa,fieldAssignments}];
	propagators = Table[Propagators[vertices, fa], {fa,fieldAssignments}];
	propagators = propagators * (-1)^ghostLoops * (-1(*nf*))^fermLoops; (*Df in gamma trace tracks number of fermion loops *)

	(* Construct amplitude *)
	amplitude = 
		Table[
			Product[
				A[ Sign[iv] * fa[[ Abs[iv]/.dict ]], 
					Abs[iv], 
					Sign[iv] * p/@(Abs[iv])
				],
				{iv,internalVertices}
			],
			{fa,fieldAssignments}
		];
	(* Pull out the denominators *)
	amplitude *= (propagators/.R[x__]:>1);
	propagators = (*Times@@*)(R[p[#]]&/@internalEdges);
	amplitude *= Product[\[Epsilon][e][\[Mu][e]],{e,externalEdges}];

	{Plus@@(Propagate/@amplitude),propagators}
];


(*
(*Example:*)
(* Note that we can have weird edge labels now *)
graph = GetGraphBasis[4,1][[1]];
graph = graph /. {2->10,-2->-10, 5->13,-5->-13, 1->9, -1->-9}
GraphInfo[graph]
DisplayGraph[ graph ]
FeynmanYM[graph,1]
*)


(* ::Subsection::Closed:: *)
(*User Interface: GetFeynYM*)


(*ClearAll[GetFeynYM];*)
GetFeynYM[n_Integer, L_Integer, minGirth_Integer:0, N_Integer:0] := Module[
{directory,name,suffix,fileName,graphs,feyn,loopMom,momRules},
	(*Search if saved*)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/feynYM/";
	name = ToString[n]<>"."<>ToString[L]<>"."<>ToString[minGirth]<>"."<>ToString[N];
	suffix = ".m";
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading Feynman expressions..."];
	PrintQ[AbsoluteTiming[Quiet[feyn = Get[fileName];];][[1]]];
	If[ feyn === $Failed,
		PrintQ["Can't find saved Feynman expressions."];
		(*If not saved build it and save it*)

		PrintQ["Generating Feynman expressions..."];
		graphs = GetGraphBasis[n,L,minGirth];
		PrintQ[AbsoluteTiming[feyn = Map[FeynmanYM[#,N]&,graphs];][[1]]];
		loopMom = ChooseLoopMom/@graphs; (* Choose best loop momenta*)
		(* Remove dependent momenta*)
		momRules = Table[ MomRules[graphs[[i]],loopMom[[i,1]]], {i,1,Length[graphs]} ];
		feyn = Table[ Relab[ momRules[[i]], feyn[[i]] ], {i,1,Length[feyn]} ];
		(* Appy on-shell, transversality, etc, too!! *)
		(* Might want to remove this for some purposes! *)
		(* Keep R[] propagator structure intact though *)
		feyn = Table[ {MomConsN[n, feyn[[i,1]]], feyn[[i,2]]/.momRules[[i]] }, {i,1,Length[feyn]}];

		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[feyn,fileName];
	];
	Return[feyn];
];


(* ::Section::Closed:: *)
(*NLSM Feynman Rules*)


(* ::Subsection::Closed:: *)
(*Feynman Rules*)


(*The c1, c2 method is not reliable because the 
  prefactors are not multiplicative! See 10-point*)
(*ClearAll[feynNLSM];*)
FeynNLSM[g_]:=Module[
{validGraphs,kinReps,nums,propPrefactors,expansionCoefs,invProps,num,props,loopMom,momRules},
	validGraphs = getValidGraphs[g];
	kinReps = getKinReps[g];
	nums = (evalValidGraph[g,#]&/@validGraphs)/.kinReps;
	propPrefactors = getPropPrefactor[g,#]&/@validGraphs;
	expansionCoefs = getCoefficients[g,#]&/@validGraphs;
	invProps = (getInvPropagators[g,#]&/@validGraphs)/.kinReps;
	num = Plus@@( invProps * propPrefactors * expansionCoefs * nums );
	num = Expand[num];

	loopMom = ChooseLoopMom[g]; (* Choose best loop momenta*)
	momRules = MomRules[g,loopMom[[1]]];(* Remove dependent momenta*)
	props = ExpandAll[getRPropagators[g]/.momRules];

	{num,props}
];


ClearAll[getPropagators];
getPropagators[g_]:=Module[{extEdges,verts,intEdges},
	{extEdges,intEdges,verts} = GraphInfo[g];
	kk[#,#]&/@intEdges
];
ClearAll[getRPropagators];
getRPropagators[g_]:=Module[{extEdges,verts,intEdges},
	{extEdges,intEdges,verts} = GraphInfo[g];
	R[p[#]]&/@intEdges
];


(****Find valid graphs with valid vertices****)
ClearAll[allMemberQ];
allMemberQ[test_List, valid_List, rep_List]:=Module[{len,i},
	len = Length[test];
	For[i = 1, i <= len, i++,
		If[!MemberQ[valid,Sort[test[[i]]/.rep]],
			Return[False];
		];
	];
	Return[True];
];

ClearAll[getValidGraphs];
getValidGraphs[g_List]:=Module[{validVerts,extEdges,verts,intEdges,extReps,test,tups,validGraphs,rep,i},
	validVerts = Sort[{{"B","\[Phi]","\[Phi]"},{"A","B","\[Phi]"}}];

	{extEdges,intEdges,verts} = GraphInfo[g];
	extReps = Thread[extEdges->"\[Phi]"];

	test = Abs[verts] /. extReps;
	tups = Tuples[{"A","B","\[Phi]"},Length[intEdges]];
	validGraphs = {};
	For[i = 1, i <= Length[tups], i++,
		rep = Thread[intEdges->(tups[[i]])];
		If[allMemberQ[test,validVerts,rep],
			AppendTo[validGraphs,rep];
		];
	];
	validGraphs
];
(***********************************************)


(******** Find appropriate prefactors and inverse propagators *********)
ClearAll[particleEdges];
particleEdges[g_List,validGraph_List,partName_String] := Module[{extEdges,intEdges,verts,pos,edges},
	{extEdges,intEdges,verts} = GraphInfo[g];
	pos = Position[validGraph,partName];
	If[pos != {}, 
		pos = pos[[All,1]];
		edges = validGraph[[#,1]]&/@pos;,

		edges = {};
	];
	edges
];

ClearAll[getInvPropagators];
getInvPropagators[g_List,validGraph_List]:=Module[{},
	Times@@Complement[getPropagators[g],kk[#,#]&/@particleEdges[g,validGraph,"\[Phi]"]]
];

ClearAll[getPropPrefactor];
getPropPrefactor[g_List, validGraph_List]:= Module[{extEdges,verts,intEdges,pe,be,ae,notPhiEdges,propFactor,preFactors},
	{extEdges,intEdges,verts} = GraphInfo[g];
	pe = particleEdges[g,validGraph,"\[Phi]"];
	be = particleEdges[g,validGraph,"B"];
	ae = particleEdges[g,validGraph,"A"];
	notPhiEdges = Complement[intEdges, pe];
	(*Does the \[Phi] prop go like I/R ? Other's go like -I??*)
	propFactor = I^Length[pe] * (-I)^Length[be] * (-I)^Length[ae];
	
	propFactor
];
(***********************************************)


ClearAll[getCoefficients];
getCoefficients[g_List,validGraph_List]:=Module[
{extEdges,intEdges,verts,extReps,vertTypes,numB\[Phi]\[Phi],posAB\[Phi],
 numAB\[Phi],pe,not\[Phi]Edges,compressReps,ABstrings,B\[Phi]\[Phi]Coef,AB\[Phi]Coef},
	{extEdges,intEdges,verts} = GraphInfo[g];
	extReps = Thread[extEdges->"\[Phi]"];
	vertTypes = Sort/@(Abs[verts]/.Join[extReps,validGraph]);
	numB\[Phi]\[Phi] = Count[vertTypes,Sort[{"B","\[Phi]","\[Phi]"}]];
	posAB\[Phi] = Flatten[Position[vertTypes,Sort[{"A","B","\[Phi]"}]]];
	numAB\[Phi]= Length[posAB\[Phi]];

	pe = particleEdges[g,validGraph,"\[Phi]"];
	not\[Phi]Edges = Complement[intEdges, pe];
	compressReps = Join[
		({a___,{w___, #, x___},b___,{y___,-#,z___},c___}:>{a,{w,z,y,x},b,c}&/@not\[Phi]Edges),
		({a___,{w___,-#, x___},b___,{y___, #,z___},c___}:>{a,{w,z,y,x},b,c}&/@not\[Phi]Edges),
		({a___,{w___,-#, x___, #, y___},b___}:>{{"ABLoop","ABLoop"}}&/@not\[Phi]Edges)(*Loop which isn't allowed. Only Phi loops allowed.*),
		({a___,{w___, #, x___,-#, y___},b___}:>{{"ABLoop","ABLoop"}}&/@not\[Phi]Edges)(*Loop which isn't allowed. Only Phi loops allowed.*)
	];
	ABstrings = verts[[posAB\[Phi]]]//.compressReps;

	B\[Phi]\[Phi]Coef = (I/(Sqrt[3]*F));
	AB\[Phi]Coef[n_]:=If[Or[OddQ[n],n===0], 0, (-1)^(n/2) * (3 * 2^(n+3)) / ((n+4)! * F^n)];
	B\[Phi]\[Phi]Coef^numB\[Phi]\[Phi] * Times@@(AB\[Phi]Coef/@((Length[#]-2)&/@ABstrings))
];


ClearAll[evalVertex];
evalVertex[mom_List,vertex_List] := Module[{ord,sig,vertexOrd,momOrd},
	ord = Ordering[vertex];
	sig = Signature[ord];
	vertexOrd = vertex[[ord]];
	momOrd = mom[[ord]];
	Switch[vertexOrd,
		{"B","\[Phi]","\[Phi]"},
			I * (sig *
				 (Sign[momOrd[[3]]] * I*p[Abs[momOrd[[3]]]][\[Mu][Abs[momOrd[[1]]]]]
				- Sign[momOrd[[2]]] * I*p[Abs[momOrd[[2]]]][\[Mu][Abs[momOrd[[1]]]]])
			),
		{"A","B","\[Phi]"},
			I * sig * \[Eta][\[Mu][Abs[momOrd[[1]]]],\[Mu][Abs[momOrd[[2]]]]],
		_,
			0
	]
];
ClearAll[evalValidGraph];
evalValidGraph[g_List, validGraph_List]:=Module[{extEdges,verts,intEdges,extReps,mom,vertex},
	{extEdges,intEdges,verts} = GraphInfo[g];
	extReps = Thread[extEdges->"\[Phi]"];
	Propagate[
		Product[
			mom = verts[[i]];
			vertex = Abs[mom]/.Join[extReps,validGraph];
			evalVertex[mom,vertex],
		{i,1,Length[verts]}]
	]
];


Clear[getKinReps];
getKinReps[g_List]:= Module[{extEdges,intEdges,verts,tups,momrules,kin},
	{extEdges,intEdges,verts} = GraphInfo[g];
	tups = Tuples[Join[extEdges,intEdges],2];
	momrules = MomRules[g];
	kin = Join[kk@@@tups,s@@@tups];
	(* True argument means use s's *)
	Dispatch[Thread[kin->ExpandAll[MomConsN[Length[extEdges],ExpandAll[Relab[momrules,kin]],True]]]]
];


(* ::Section::Closed:: *)
(*Tree Amplitudes*)


(* ::Subsection::Closed:: *)
(*GetCOTreeAmp*)


ClearAll[GetCOTreeAmp];
(* sepPropsQ = True keeps propagators separate *)
GetCOTreeAmp[n_Integer, sepPropsQ_:False] := Module[
{directory,name,suffix,fileName,
amp},
	(*Search if saved*)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/ampTreeYM/";
	name = If[!sepPropsQ, ToString[n], ToString[n]<>".sepProps"];
	suffix = ".m";
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading color-ordered tree amplitude..."];
	PrintQ[AbsoluteTiming[Quiet[amp = Get[fileName];];][[1]]];
	If[ amp === $Failed,
		PrintQ["Can't find saved color-ordered tree amplitude."];
		(*If not saved build it and save it*)
		PrintQ["Generating color-ordered tree amplitude..."];
		amp = GenerateCOTreeAmp[n, sepPropsQ];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[amp,fileName];
	];

	amp
];


(* Gives color-ordered partial amplitude Atree(1,2,...,n) *)
ClearAll[GenerateCOTreeAmp];
(* sepPropsQ = True keeps propagators separate *)
GenerateCOTreeAmp[n_Integer, sepPropsQ_:False]:= Module[
{co,signs,feyn,re,amp,
props,propTypes,pows,max,lcd,numFactors},

	If[n < 4, 
		Print["n must be greater 4."];
		Return[];
	];
	co = GetCOBasis[n];
	signs = co/.fColor[__]:>1;
	co = signs*co;
	feyn = GetFeynYM[n,0];
	re[ord_List] := Join[
						Thread[ p/@Range[n] -> p/@ord ],
						Thread[ \[Epsilon]/@Range[n] -> \[Epsilon]/@ord ]
					];
	amp = ({#,1}&/@signs) * ( 
				co /. fColor[i_,ord_] :> {
					Expand[MomConsN[n, Relab[re[ord],feyn[[i,1]]], True]], 
					feyn[[i,2]]/.R[x__]:>R[Expand[Relab[re[ord],x]/.p[n]->-Sum[p[j],{j,1,n-1}]]]
					}
			);
	(* Put over common denominator *)
	props = amp[[All,2]];
	propTypes = GetVars[props];
	pows = Table[ Count[props[[i]], #]&/@propTypes, {i,1,Length[props]} ];
	max = (Max /@ Transpose[pows]);
	lcd = propTypes ^ (Max /@ Transpose[pows]);
	numFactors = Factor[MomConsN[n, 
					Table[ Times@@( propTypes^(max - pows[[i]]) ), {i,1,Length[props]} ]/.R->S,
					True
				]];
	amp = Expand[ Plus@@(amp[[All,1]] * numFactors) ];
	
	If[sepPropsQ,
		{amp, lcd},
		
		lcd = Factor[MomConsN[n,(Times@@lcd)/.R->S, True]];
		amp/lcd
	]
];


(* ::Subsection::Closed:: *)
(*GetNLSMCOAmp*)


ClearAll[GetNLSMCOAmp];
(* sepPropsQ = True keeps propagators separate *)
GetNLSMCOTreeAmp[n_Integer, sepPropsQ_:False] := Module[
{directory,name,suffix,fileName,
amp},
	(*Search if saved*)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/ampTreeNLSM/";
	name = If[!sepPropsQ, ToString[n], ToString[n]<>".sepProps"];
	suffix = ".m";
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading color-ordered NLSM tree amplitude..."];
	PrintQ[AbsoluteTiming[Quiet[amp = Get[fileName];];][[1]]];
	If[ amp === $Failed,
		PrintQ["Can't find saved color-ordered NLSM tree amplitude."];
		(*If not saved, build it and save it*)
		PrintQ["Generating color-ordered NLSM tree amplitude..."];
		amp = GenerateNLSMCOTreeAmp[n, sepPropsQ];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[amp,fileName];
	];

	amp
];


(* Gives color-ordered partial amplitude Atree(1,2,...,n) *)
ClearAll[GenerateNLSMCOTreeAmp];
(* sepPropsQ = True keeps propagators separate *)
GenerateNLSMCOTreeAmp[n_Integer, sepPropsQ_:False]:= Module[
{co,signs,graphs,feyn,re,vars,iGraph,reps,ord,amp,
props,propTypes,pows,max,lcd,numFactors},

	If[n < 4, 
		Print["n must be greater or equal to 4."];
		Return[];
	];
	If[OddQ[n],
		Return[0];
	];
	co = GetCOBasis[n];
	signs = co/.fColor[__]:>1;
	co = signs*co;
	graphs = GetGraphBasis[n,0];
	feyn = FeynNLSM/@graphs;
	re[ord_List] := Thread[ p/@Range[n] -> p/@ord ];
	vars = GetVars[feyn];
	amp = Table[
			{iGraph, ord} = List@@(co[[i]]); (* fColor[iGraph, ord] *)
			reps = Dispatch[Thread[vars -> Expand[MomConsN[n,Relab[re[ord],vars],True]]]];
			{signs[[i]],1} * {Expand[feyn[[iGraph,1]]/.reps],
							feyn[[iGraph,2]]/.R[x__]:>R[Expand[Relab[re[ord],x]/.p[n]->-Sum[p[j],{j,1,n-1}]]]},
		{i,1,Length[co]}];
	amp = DeleteCases[amp,{0,__}];

	(* Put over common denominator *)
	props = amp[[All,2]];
	propTypes = GetVars[props];
	pows = Table[ Count[props[[i]], #]&/@propTypes, {i,1,Length[props]} ];
	max = (Max /@ Transpose[pows]);
	lcd = propTypes ^ (Max /@ Transpose[pows]);
	numFactors = Factor[MomConsN[n, 
					Table[ Times@@( propTypes^(max - pows[[i]]) ), {i,1,Length[props]} ]/.R[x_]:>S[x],
					True
				]];
	amp = Expand[ Plus@@(amp[[All,1]] * numFactors) ];

	If[sepPropsQ,
		{amp, lcd},
		
		lcd = Factor[MomConsN[n,(Times@@lcd)/.R[x_]:>S[x],True]];
		amp/lcd
	]
];


(* ::Subsection::Closed:: *)
(*User Interface: Atree*)


(* Color-ordered YM tree amplitude *)
(*ClearAll[Atree];*)
(*Use negative arguements for incoming momentum*)
(*eg. legList={1,2,-6,3,-5,7}*)
(* sepPropsQ = True keeps propagators separate *)
Atree[legList_List, sepPropsQ_:False]:= Module[
{n,amp,rep},
	n = Length[legList];
	amp = GetCOTreeAmp[n,sepPropsQ]; (* Gives Atree(1,2,...,n) *)
	(* Don't bother relabeling if the legList is (1,2,...,n) ordering*)
	If[ legList === Range[n],
		Return[amp]
	];
	rep = Join[
		Thread[p/@Range[n] -> Sign[legList]*p/@Abs[legList]],
		Thread[\[Epsilon]/@Range[n] -> \[Epsilon]/@Abs[legList]]
	];

	(* s's will get changed to kk's, but we should kill k^s since we assume on-shell *)
	Relab[rep, amp]/.kk[i_,i_]:>0
];

Atree[n_Integer, sepPropsQ_:False] := Atree[Range[n], sepPropsQ];


(* ::Subsection::Closed:: *)
(*User Interface: NLSMtree*)


(* Color-ordered NLSM tree amplitude *)
NLSMtree[legList_List, sepPropsQ_:False]:= Module[
{n,amp,rep},
	n = Length[legList];
	amp = GetNLSMCOTreeAmp[n,sepPropsQ]; (* Gives NLSMtree(1,2,...,n) *)
	(* Don't bother relabeling if the legList is (1,2,...,n) ordering*)
	If[ legList === Range[n],
		Return[amp]
	];
	rep = Thread[p/@Range[n] -> Sign[legList]*p/@Abs[legList]];

	(* s's may get changed to kk's, but we should kill k^s since we assume on-shell *)
	Relab[rep, amp]/.kk[i_,i_]:>0
];

NLSMtree[n_Integer, sepPropsQ_:False] := NLSMtree[Range[n], sepPropsQ];


(* ::Subsection::Closed:: *)
(*User Interface: Atree for 3-point and Common Cases*)


(*Tool to build trees*)
(*ClearAll[Tree3]*)
(*Use negative arguements for incoming momentum*)
(*indexLabel is used to distinguish multiple copies of the three vertex*)
Tree3[vertList_List,indexLabel_]:=Module[{nPoints,input,mom,index,ext,extPol,vert3,col},
	nPoints=4;
	
	input=vertList;
	mom=p/@input;
	index=Abs/@input;

	ext=Select[index,#<=nPoints&];
	extPol = Times@@(((\[Epsilon][#][\[Mu][#,indexLabel]])&/@ext)/.{\[Epsilon][x_]:>\[Epsilon][Abs[x]]/;Negative[x]});
	
	(*Dont apply mom cons or anything (expect simple stuff) since we are using numerics*)
	vert3=1/(I Sqrt[2]) * (
		(mom[[1]][ \[Mu][ index[[3]],indexLabel ] ] - (mom[[2]][ \[Mu][ index[[3]],indexLabel ] ])) \[Eta][\[Mu][ index[[1]],indexLabel],\[Mu][ index[[2]],indexLabel]] +
		(mom[[2]][ \[Mu][ index[[1]],indexLabel ] ] - (mom[[3]][ \[Mu][ index[[1]],indexLabel ] ])) \[Eta][\[Mu][ index[[2]],indexLabel],\[Mu][ index[[3]],indexLabel]] + 
		(mom[[3]][ \[Mu][ index[[2]],indexLabel ] ] - (mom[[1]][ \[Mu][ index[[2]],indexLabel ] ])) \[Eta][\[Mu][ index[[3]],indexLabel],\[Mu][ index[[1]],indexLabel]] );

	Propagate[extPol*vert3/.{p[x_][y__]:>-p[Abs[x]][y]/;Negative[x]}]
];


(* These function names are used in old code *)
ClearAll[Atree4,Atree5];
Atree4[legList_List] := Atree[legList];
Atree5[legList_List] := Atree[legList];


(* ::Subsection::Closed:: *)
(*User Interface: Gravity Tree Amplitudes from KLT *)
(*(3-, 4-, and 5-point)*)


(*Tool to build gravity 3-vertex*)
(*Clear[Mtree3]*)
(*Use negative arguements for incoming momentum*)
(*indexLabelL, indexLabelR are used to distinguish the gravity double copies*)
Mtree3[vertList_List,indexLabelL_,indexLabelR_]:=Module[{nPoints,input,mom,index,ext,extPolL,extPolR,vert3},
	nPoints=4;
	
	input=vertList;
	mom=p/@input;
	index=Abs/@input;

	ext=Select[index,#<=nPoints&];
	extPolL = Times@@(((\[Epsilon][#][\[Mu][#,indexLabelL]])&/@ext)/.{\[Epsilon][x_]:>\[Epsilon][Abs[x]]/;Negative[x]});
	extPolR = Times@@(((\[Epsilon]2[#][\[Mu][#,indexLabelR]])&/@ext)/.{\[Epsilon]2[x_]:>\[Epsilon][Abs[x]]/;Negative[x]});
	
	(*Dont apply mom cons or anything (expect simple stuff) since we are using numerics*)
	vert3[indexLabel_]:=1/(I Sqrt[2]) * (
		(mom[[1]][ \[Mu][ index[[3]],indexLabel ] ] - (mom[[2]][ \[Mu][ index[[3]],indexLabel ] ])) \[Eta][\[Mu][ index[[1]],indexLabel],\[Mu][ index[[2]],indexLabel]] +
		(mom[[2]][ \[Mu][ index[[1]],indexLabel ] ] - (mom[[3]][ \[Mu][ index[[1]],indexLabel ] ])) \[Eta][\[Mu][ index[[2]],indexLabel],\[Mu][ index[[3]],indexLabel]] + 
		(mom[[3]][ \[Mu][ index[[2]],indexLabel ] ] - (mom[[1]][ \[Mu][ index[[2]],indexLabel ] ])) \[Eta][\[Mu][ index[[3]],indexLabel],\[Mu][ index[[1]],indexLabel]] );

	I times[Propagate[extPolL*vert3[indexLabelL]/.{p[x_][y__]:>-p[Abs[x]][y]/;Negative[x]}],Propagate[extPolR*vert3[indexLabelR]/.{p[x_][y__]:>-p[Abs[x]][y]/;Negative[x]}]]
];


(*Tool to build tree amplitude*)
(*Clear[Mtree4]*)
(*Use negative arguements for incoming momentum*)
(*eg. legList={1,2,-5,7}*)
Mtree4[legList_List]:=Module[{a1234,a1243,R0},

	(*Not actually A(1,2,3,4) but A(legList[[1]],legList[[2]],...)*)
	a1234=Atree4[legList];
	a1243=Atree4[legList[[{1,2,4,3}]]];
	(*Distinguish copies*)
	a1243=a1243/.Dispatch[{ee->EE,ke->kE,le->lE}];

	R0[i_,j_]:=R[Sign[legList[[i]]]p[Abs[legList[[i]]]]+Sign[legList[[j]]]p[Abs[legList[[j]]]]]//.{R[x__]:>AngleBracket[x,x],kk[x_,x_]:>0};

	(*Don't immediately multiply since tensor size would double. Want to keep L and R YM copies separate*)
	-I*R0[1,2]*times[a1234,a1243]
];
Mtree[4] := Mtree4[Range[4]];


(*Tool to build tree amplitude*)
(*Clear[Mtree5]*)
(*Use negative arguements for incoming momentum*)
(*eg. legList={3,4,-5,-5,-7}*)
Mtree5[legList_List]:=Module[{a12345,a21435,a13245,a31425,R0},

	(*Not actually A(1,2,3,4,5) but A(legList[[1]],legList[[2]],...)*)
	a12345=Atree5[legList];
	a21435=Atree5[legList[[{2,1,4,3,5}]]];
	a13245=Atree5[legList[[{1,3,2,4,5}]]];
	a31425=Atree5[legList[[{3,1,4,2,5}]]];

	(*Distinguish copies*)
	a21435=a21435/.Dispatch[{ee->EE,ke->kE,le->lE}];
	a31425=a31425/.Dispatch[{ee->EE,ke->kE,le->lE}];

	R0[i_,j_]:=R[Sign[legList[[i]]]p[Abs[legList[[i]]]]+Sign[legList[[j]]]p[Abs[legList[[j]]]]]//.{R[x__]:>AngleBracket[x,x],kk[x_,x_]:>0};

	(*Don't immediately multiply since tensor size would double. Want to keep L and R YM copies separate*)
	I(R0[1,2]*R0[3,4]*times[a12345,a21435] + R0[1,3]*R0[2,4]*times[a13245,a31425])
];
Mtree[5] := Mtree5[Range[5]];


(* ::Subsection::Closed:: *)
(*User Interface: Strip Internal Polarization*)
(*(Used for sewing trees to form 4-point loops)*)


(*strip polarization for the internal legs, \[Epsilon][i>4]*)
(*Use an extra index label that denotes what tree 
amplitude we are working with since for gravity we 
will have multiple copies of say \[Mu][5] (4 copies at one loop)*)
(*ClearAll[StripPol4]*)
StripPol4[expr_,indexLabel_]:=Module[{},
	StripPol[4,expr,indexLabel]
];


StripPol[n_,expr_,indexLabel_]:=Module[{stripRep},
	stripRep=
		{ke[i_,j_]:>p[i][\[Mu][j,indexLabel]]/;Abs[j]>n,
		ee[i_,j_]:>\[Eta][\[Mu][i,indexLabel],\[Mu][j,indexLabel]]/;(Abs[i]>n && Abs[j]>n),
		ee[i_,j_]:>\[Epsilon][j][\[Mu][i,indexLabel]]/;(Abs[i]>n && Abs[j]<=n)};
	
	expr /. stripRep
];


(*strip polarization for the internal legs, \[Epsilon][i>4]*)
(*Use an extra index label that denotes what tree 
amplitude we are working with since for gravity we 
will have multiple copies of say \[Mu][5] (4 copies at one loop)*)
(*ClearAll[StripPol4LR];*)
StripPol4LR[expr_,indexLabelL_,indexLabelR_]:=Module[{nPoints,stripRep},
	StripPolLR[4,expr,indexLabelL,indexLabelR]
];


StripPolLR[n_Integer,expr_,indexLabelL_,indexLabelR_]:=Module[{stripRep},
	nPoints=n;

	stripRep={
		ke[i_,j_]:>p[i][\[Mu][j,indexLabelL]]/;Abs[j]>n,
		ee[i_,j_]:>\[Eta][\[Mu][i,indexLabelL],\[Mu][j,indexLabelL]]/;(Abs[i]>n && Abs[j]>n),
		ee[i_,j_]:>\[Epsilon][j][\[Mu][i,indexLabelL]]/;(Abs[i]>n && Abs[j]<=n),
		kE[i_,j_]:>p[i][\[Mu][j,indexLabelR]]/;Abs[j]>n,
		EE[i_,j_]:>\[Eta][\[Mu][i,indexLabelR],\[Mu][j,indexLabelR]]/;(Abs[i]>n && Abs[j]>n),
		EE[i_,j_]:>\[Epsilon]2[j][\[Mu][i,indexLabelR]]/;(Abs[i]>n && Abs[j]<=n),
		eE[i_,j_]:>\[Eta][\[Mu][i,indexLabelL],\[Mu][j,indexLabelR]]/;(Abs[i]>n && Abs[j]>n),
		eE[i_,j_]:>\[Epsilon]2[j][\[Mu][i,indexLabelL]]/;(Abs[i]>n && Abs[j]<=n),
		eE[i_,j_]:>\[Epsilon][i][\[Mu][j,indexLabelR]]/;(Abs[j]>n && Abs[i]<=n)
	};
	
	expr /. stripRep
];


(* ::Section::Closed:: *)
(*Tensor Reduction*)


(* See "TensorReduction.pdf" in aux/tensorReduction for example and explanation *)


(* ::Subsection::Closed:: *)
(*Split by Twos*)


ClearAll[GetSBT];
GetSBT[n_Integer]:= Block[
{ran,directory,name,suffix,fileName,sbt},
	If[OddQ[n], Print["Odd number of momenta. Returning NULL..."];Return[];];

	ran = Range[n];
	(*Search if saved*)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/splitByTwos/";
	name = ToString[n];
	suffix = ".m";
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading split-by-twos info..."];
	PrintQ[AbsoluteTiming[Quiet[sbt = Get[fileName];];][[1]]];
	If[ sbt === $Failed,
		PrintQ["Can't find cached split-by-twos."];
		(*If not saved build it and save it*)
		PrintQ["Splitting by twos..."];
		PrintQ[AbsoluteTiming[sbt = SplitByTwos[ran];][[1]]];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[sbt,fileName];
	];
	Return[sbt];
];


(* 2-14 points take < 5s, 16-point takes ~ 1 min *)
(* From Tristan Dennen *)
ClearAll[SplitByTwos];
SplitByTwos[arg_List]:=  
	Catch[Block[
		{list=arg,rlist = Rest[arg],pairs,set},
		If[Length[list]==2,Throw[{{list}}];];
		pairs = Table[Join[{{list[[1]],list[[i]]}},#]&/@SplitByTwos[Drop[rlist,{i-1}]],{i,2,Length[list]}];
		Throw[Flatten[pairs,1]];
	set
	]];


(* ::Subsection::Closed:: *)
(*Tensor Structures that Obey Symmetries*)


(*Example: input = {3,3,2} ~ pppqqqrr *)
(* Should sort input labels so, 2.3.3 & 3.3.2 aren't separately generated *)
ClearAll[GetTensorStructs];
GetTensorStructs[input_List]:=Block[
{n,ran,ids,directory,suffix,name,fileName,sbt,tens,i},
	n = Plus@@input;
	ran = Range[n];
	ids = {};
	For[i = 1, i <= Length[input], i++,
		AppendTo[ids,ran[[1;;input[[i]]]]];
		ran = Drop[ran,input[[i]]];
	];
	ran = Range[n];

	If[OddQ[n], PrintQ["Odd number of momenta. Returning NULL..."];Return[];];

	(* Basically, splitting products of \[Eta]'s into equivalence classes based on symmetry *)
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/tensorStructs/";
	suffix = ".m";
	name = StringJoin@@Riffle[ToString/@input,"."];
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading tensors that obey input symmetry..."];
	PrintQ[AbsoluteTiming[Quiet[tens = Get[fileName];];][[1]]];
	If[ tens === $Failed,
		PrintQ["Can't find cached tensor structures."];
		(*If not saved build it and save it*)
		(**************************************)
		(************* SplitByTwos ***********)
		(**************************************)
		sbt = GetSBT[n];
		(**************************************)
		(**** Tensor symmetry structures *****)
		(**************************************)
		PrintQ["Building tensors (2 steps)..."];
		PrintQ[AbsoluteTiming[tens = GetIdenticalStructs[sbt, ids];][[1]]];
		PrintQ[AbsoluteTiming[tens = ParallelTable[(sbt[[tens[[i,j]]]]),{i,1,Length[tens]},{j,1,Length[tens[[i]]]},DistributedContexts :> distributedContexts];][[1]]];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[tens,fileName];
	];
	Return[tens];
];


ClearAll[GetIdenticalStructs];
GetIdenticalStructs[sbt_List, identifications_List] := Block[
{reps,dupPos},
	If[ identifications === {},
		Return[{}];
	];
	reps = (Sort/@identifications);
	reps = Flatten[Thread[Drop[#,1]->#[[1]]]&/@reps];

	dupPos = (Sort[Sort/@#])&/@(sbt//.reps);
	dupPos = PositionDuplicates[dupPos];
	Return[dupPos];
];


ClearAll[PositionDuplicates];
PositionDuplicates[list_] := GatherBy[Range@Length[list], list[[#]] &];


(* ::Subsection::Closed:: *)
(*Scalar Products from Tensor Structures*)


(* Each representative from a tensor structure contracts the input momenta in a unique way *)
(* Example: input = {3,3,2} ~ pppqqqrr *)
(* Should sort input labels so, 2.3.3 & 3.3.2 aren't separately generated *)
ClearAll[GetScalarProducts];
GetScalarProducts[input_List] := Block[
{n,ran,ids,directory,suffix,name,fileName,sps,reps,i},
	n = Plus@@input;
	ran = Range[n];
	ids = {};
	For[i = 1, i <= Length[input], i++,
		AppendTo[ids,ran[[1;;input[[i]]]]];
		ran = Drop[ran,input[[i]]];
	];
	ran = Range[n];

	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/scalarProducts/";
	suffix = ".m";
	name = StringJoin@@Riffle[ToString/@input,"."];
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading scalar products..."];
	PrintQ[AbsoluteTiming[Quiet[sps = Get[fileName];];][[1]]];
	If[ sps === $Failed,
		PrintQ["Can't find cached scalar products."];
		(**************************************)
		(**** Tensor symmetry structures *****)
		(**************************************)
		sps = GetTensorStructs[input];
		(* Only need one representative *)
		sps = sps[[All,1]];
		(**************************************)
		(********** Scalar Products **********)
		(**************************************)
		PrintQ["Building scalar products..."];
		reps = #->Position[ids,#][[1,1]]&/@Flatten[ids];
		PrintQ[AbsoluteTiming[sps = (Sort[Sort/@#])&/@(sps/.reps);][[1]]];
	
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[sps,fileName];
	];
	Return[sps];
];


(* ::Subsection::Closed:: *)
(*Coefficients from Tensor Contractions*)


ClearAll[GetContractions];
(* Out: Sum of tens[[i]] being contracted by any representative of tens[[j]] *)
(* All representatives should give same contraction (if tens[[i]] is summed) *)
(* Example: input = {3,3,2} ~ pppqqqrr *)
(* Should sort input labels so, 2.3.3 & 3.3.2 aren't separately generated *)
GetContractions[input_List] := Block[
{directory,suffix,name,fileName,coefs,tens},

	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/tensorContractions/";
	suffix = ".m";
	name = StringJoin@@Riffle[ToString/@input,"."];
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading tensor contractions..."];
	PrintQ[AbsoluteTiming[Quiet[coefs = Get[fileName];];][[1]]];
	If[ coefs === $Failed,
		PrintQ["Can't find cached tensor contractions."];
		(**************************************)
		(**** Tensor symmetry structures *****)
		(**************************************)
		tens = GetTensorStructs[input];
		(**************************************)
		(****** Contraction of tensors *******)
		(**************************************)
		PrintQ["Finding all relevant contractions to build equations..."];
		PrintQ[AbsoluteTiming[
			coefs = 
				ParallelTable[
					Sum[ 
						Ds^ContractionCycles[tens[[i,a]],tens[[j,1]]]
					,{a,1,Length[tens[[i]]]}]
				,{i,1,Length[tens]},{j,1,Length[tens]}
				,DistributedContexts :> distributedContexts
				];
		][[1]]];
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[coefs,fileName];
	];
	Return[coefs];
];


(* This is the workhorse that runs thousands or millions of times: 
		~(# of tensor structures)*(length of tensor structure) *)
(* Need speed with good constants *)
(* Problem is same as finding connected components in undirected graphs *)
(* Note: Took out error checking for speed *)
ClearAll[ContractionCycles];
ContractionCycles = Compile[{
{pairList, _Integer, 2},
{pairsToContract, _Integer, 2}},
Block[{edges,len,nodes,i,n1,n2,numCC,next,e1,e2},
	(*If[ !(Sort[Flatten[pairList]] === Sort[Flatten[pairsToContract]]),
		Print["Need complete and valid list of contractions."];
		numCC = -1;,*)

		edges = Join[pairList,pairsToContract];
		len = Length[edges];
		(* numNodes === numEdges since every vertex has 2 edges *)
		(* {outEdge1, outEdge2, visited?} *)
		nodes = ConstantArray[{0,0,0},len];
		For[i = 1, i <= len, i++,
			n1 = edges[[i,1]];
			n2 = edges[[i,2]];
			If[nodes[[n1,1]]===0, nodes[[n1,1]] = n2, nodes[[n1,2]] = n2];
			If[nodes[[n2,1]]===0, nodes[[n2,1]] = n1, nodes[[n2,2]] = n1];
		];

		(* Number of connected components *)
		numCC = 0;
		For[i = 1, i <= len, i++,
			(* If unvisited, sorta BFS *) 
			(* Go around circle since only two edges per vertex *)
			next = i;
			If[ nodes[[next,3]] === 0,
				numCC++;
				nodes[[next,3]] = 1;
				e1 = nodes[[next,1]];
				e2 = nodes[[next,2]];
				While[Or[nodes[[e1,3]]===0,nodes[[e2,3]]===0],
					next = If[nodes[[e1,3]]===0, e1, e2];
					nodes[[next,3]] = 1;
					e1 = nodes[[next,1]];
					e2 = nodes[[next,2]];
				];
			 ];
		];
	(*];*)
numCC
],
CompilationTarget -> "C",
RuntimeAttributes -> Listable
];


(*(****Over 10x's slower using Mathematica functions******)
ClearAll[Contract2];
Contract2 = Compile[{
{pairList, _Integer, 2},
{pairsToContract, _Integer, 2}},
	If[ !(Sort[Flatten[pairList]] === Sort[Flatten[pairsToContract]]),
		Print["Need complete and valid list of contractions."];
		-1,

		Length[ConnectedComponents[FromUnorderedPairs[Join[pairList,pairsToContract]]]]
	],
	RuntimeAttributes -> Listable
];
*)


(* ::Subsection::Closed:: *)
(*Equation Building*)


(* Each scalar product has a sum of tensor structures in front of it *)
(* This function generates equations for the coefficients in front of each structure for each scalar product *)
(* Number of equations and coefficients will be (# of tensor structures)^2 *)
(* Example: input = {3,3,2} ~ pppqqqrr *)
(* Should sort input labels so, 2.3.3 & 3.3.2 aren't separately generated *)
ClearAll[GetEqs];
GetEqs[input_List] := Block[
{directory,suffix,name,fileName,eqs,contractCoefs},
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/equations/";
	suffix = ".m";
	name = StringJoin@@Riffle[ToString/@input,"."];
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading equations..."];
	PrintQ[AbsoluteTiming[Quiet[eqs = Get[fileName];];][[1]]];
	If[ eqs === $Failed,
		PrintQ["Can't find cached equations."];
		(**************************************)
		(******** Tensor Contractions ********)
		(**************************************)
		contractCoefs = GetContractions[input];
		(**************************************)
		(************* Equations **************)
		(**************************************)
		PrintQ["Building equations..."];
		PrintQ[AbsoluteTiming[
		eqs = Table[
				Table[
					Sum[a[j,i]*contractCoefs[[i,k]],{i,1,Length[contractCoefs]}]== KroneckerDelta[j,k],
				{j,1,Length[contractCoefs]}],
			{k,1,Length[contractCoefs]}];][[1]]];
		eqs = Transpose[eqs];
	
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[eqs,fileName];
	];
	Return[eqs];
];


(* ::Subsection::Closed:: *)
(*Solver*)


(* Example: input = {3,3,2} ~ pppqqqrr *)
ClearAll[GetSoln];
GetSoln[input_List] := Block[
{directory,suffix,name,fileName,soln,eqs,slen,contractCoefs},
	directory = MASTERFUNCTIONSDIRECTORY<>"aux/tensorReduction/solutions/";
	suffix = ".m";
	name = StringJoin@@Riffle[ToString/@input,"."];
	fileName = directory<>name<>suffix;
	SetDirectory[NotebookDirectory[]];
	PrintQ["Loading solutions..."];
	PrintQ[AbsoluteTiming[Quiet[soln = Get[fileName];];][[1]]];
	If[ soln === $Failed,
		PrintQ["Can't find cached solutions."];
		(**************************************)
		(******** Tensor Contractions ********)
		(**************************************)
		eqs = GetEqs[input];
		(**************************************)
		(************* Equations **************)
		(**************************************)
		slen = ToString[Length[eqs]];
		PrintQ["Solving "<>slen<>"x"<>slen<>" equations..."];
		PrintQ[AbsoluteTiming[
		soln = Factor[
			ParallelTable[
				Solve[eqs[[i]],GetVars[eqs[[i,All,1]],{Ds}]][[1]],
			{i,1,Length[eqs]}
			,DistributedContexts :> distributedContexts]
		];][[1]]];
	
		SetDirectory[NotebookDirectory[]];
		Quiet[DeleteFile[fileName];];
		Put[soln,fileName];
	];
	Return[soln];
];


(* ::Subsection::Closed:: *)
(*User Interface: TensorReduce*)


(* Can handle a sum of inputs. We send each to TensorReduce0. *)
(* ClearAll[TensorReduce]; *)
TensorReduce[in_] := Block[
{tensors,vars,pres,printCell,tenCount,tenOrder},
	tensors = Expand[in];
	vars = Complement[GetVars[tensors],GetVars[tensors,{lk[_,_],le[_,_],lE[_,_]}]];
	If[vars === {},(*No tensors*)
		Return[{{tensors,1}}];
	];
	tensors = CoefficientRules[tensors,vars];
	pres = tensors[[All,2]];
	tensors = tensors[[All,1]];
	tensors = Table[
		Flatten[Table[ConstantArray[vars[[j]],tensors[[i,j]]],{j,1,Length[tensors[[i]]]}]]
	,{i,1,Length[tensors]}];
	
	tenCount = 0;
	tenOrder = 0;
	printCell = PrintTemporary[
					"Reducing tensor ",Dynamic[tenCount]," of "<>ToString[Length[tensors]],"\n",
					"Tensor length: ",Dynamic[tenOrder]
				];
	tensors = Table[
				tenCount = i;
				tenOrder = Length[tensors[[i]]];
				TensorReduce0[tensors[[i]],pres[[i]]]
			,{i,1,Length[tensors]}];
	NotebookDelete[printCell];

	Join@@tensors			
];


(* This can only handle a list of {le[5,1], lk[6,2], lk[6,2],...} *)
ClearAll[TensorReduce0];
TensorReduce0[in_List,pre_] := Block[
{toRed = in,input,lReps,indexReps,sps,tens,soln,red},

	If[toRed === {},
		Return[{{pre,1}}];(* If no tensors...*)
	];
	toRed = GatherBy[toRed,#/.{le[i_,j_]:>i, lk[i_,j_]:>i, lE[i_,j_]:>i}&];
	toRed = Reverse[SortBy[toRed,Length]];

	input = Length/@toRed;
	If[OddQ[Plus@@input],
		Return[{{0,0}}]; (* If odd number of tensors, return 0 *)
	];

	lReps = Thread[ Range[Length[input]] -> Table[toRed[[i,1]]/.{lk[i_,j_]:>i,le[i_,j_]:>i,lE[i_,j_]:>i}, {i,1,Length[toRed]}] ];
	lReps = Dispatch[(# -> ll@@(#/.lReps)) &/@ Tuples[Range[Length[input]],2]];

	indexReps = Thread[ Range[Length[Flatten[toRed]]] -> (Flatten[toRed]/.{lk[i_,j_]:>p[j],le[i_,j_]:>\[Epsilon][j],lE[i_,j_]:>\[Epsilon]2[j]}) ];
	indexReps = Dispatch[(# -> ((AngleBracket@@(#/.indexReps))/. ZEROREPS)) &/@ Tuples[Range[Length[indexReps]],2]];

	sps = GetScalarProducts[input];
	sps = Times@@@(sps /. lReps);

	tens = GetTensorStructs[input];
	tens = Table[Plus@@(Times@@@(tens[[i]] /. indexReps)),{i,1,Length[tens]}];

	soln = GetSoln[input];
	soln = soln/.Rule[x_,y_]:>y;

	red = Table[{pre*sps[[i]], (soln.tens)[[i]]},{i,1,Length[sps]}];

	red
];


(* ::Section::Closed:: *)
(*Numerics*)


(* ::Subsection::Closed:: *)
(*Random Component Generatation*)


(*Random Component for a vector*)
(*Try larger values*)
ClearAll[RandComp]
RandComp[wp_:$MinPrecision]:=Sign[RandomReal[{-1,1}]]*RandomReal[{1,randSize},WorkingPrecision->wp];

(*Random complex component for a vector*)
ClearAll[RandCompC]
RandCompC[wp_:$MinPrecision]:=Sign[RandomReal[{-1,1}]]*RandomComplex[{-randSize-randSize I,randSize+ randSize I},WorkingPrecision->wp];


(* ::Subsection::Closed:: *)
(*Metric Matrix*)


(*Metric in +------... signature*)
g[i_,j_]:= (-1)^Boole[i>1] KroneckerDelta[i,j];

Mink[a_List, b_List]:= If[Length[a]===Length[b], Sum[a[[i]]*g[i,j]*b[[j]],{i,1,Length[a]},{j,1,Length[b]}]];


(* ::Subsection::Closed:: *)
(*4-D Spinor Technology*)


(*Follow Srednicki (except for signature of metric)*)
(*Britto: http://ipht.cea.fr/Docspht/articles/t11/138/public/Cours2011.pdf*)
(*S@M: http://arxiv.org/pdf/0710.5559v2.pdf.*)


(*\[Epsilon]^ab: upper indices*)
(*\[Epsilon]^12=1*)
epMetricUp=({
 {0, 1},
 {-1, 0}
});
(*Subscript[\[Epsilon], ab]: lower indices*)
epMetricDown=({
 {0, -1},
 {1, 0}
});

(*Contraction convection for dotted and undotted: Subscript[Overscript[\[Lambda], ~], adot] Overscript[\[Lambda], ~]^adot, \[Lambda]^a Subscript[\[Lambda], a]*)
(*Lorentz indices are up*)

(*Subscript[u, -](Subscript[p, i]) = |(Subscript[p, i]^-)> = |i] ~ Subscript[\[Lambda], a]    -    The last relationship disagrees with Britto and S@M*)
(*See Notes 2/3/2014 pg .2*)
(*Expression from Britto or convert Srednicki (50.8) to cartesian (p0,p1,p2,p3)*)
spinLambdaDown[vec_]:={-(vec[[2]]-I*vec[[3]])/Sqrt[vec[[1]]+vec[[4]]],(vec[[1]]+vec[[4]])/Sqrt[vec[[1]]+vec[[4]]]};
(*[i] ~ (\[Lambda]^a) *)
spinLambdaUp[vec_]:=Table[Sum[epMetricUp[[a,b]]spinLambdaDown[vec][[b]],{b,1,2}],{a,1,2}];
(*<i| ~ Subscript[Overscript[\[Lambda], ~], adot]*)
spinLambdaTildeDown[vec_]:={-(vec[[2]]+I*vec[[3]])/Sqrt[vec[[1]]+vec[[4]]],(vec[[1]]+vec[[4]])/Sqrt[vec[[1]]+vec[[4]]]};
(*Subscript[u, +](Subscript[p, i]) = |(Subscript[p, i]^+)> = |i> ~ (Overscript[\[Lambda], ~]^adot)    -    The last relationship disagrees with Britto and S@M*)
spinLambdaTildeUp[vec_]:=Table[Sum[epMetricUp[[adot,bdot]]spinLambdaTildeDown[vec][[bdot]],{bdot,1,2}],{adot,1,2}];


(*square bracket right*)sbr[vec_]:=spinLambdaDown[vec];
(*square bracket left*)sbl[vec_]:=spinLambdaUp[vec];
(*angle bracket right*)abr[vec_]:=spinLambdaTildeUp[vec];
(*angle bracket left*)abl[vec_]:=spinLambdaTildeDown[vec];

(*square bracket*) sb[i1_Integer,i2_Integer]:=sbl[VEC[i1]].sbr[VEC[i2]];
(*angle bracket*) ab[i1_Integer,i2_Integer]:=abl[VEC[i1]].abr[VEC[i2]];
(*square bracket*) sb[vec1_List,vec2_List]:=sbl[vec1].sbr[vec2];
(*angle bracket*) ab[vec1_List,vec2_List]:=abl[vec1].abr[vec2];


(*Subscript[(\[Sigma]^\[Mu]), a adot]*)
sigma={({
 {1, 0},
 {0, 1}
}),({
 {0, 1},
 {1, 0}
}),({
 {0, -I},
 {I, 0}
}),({
 {1, 0},
 {0, -1}
})};
(* \[Sigma]^(\[Mu] adot a) = (I, -\[Sigma]) *)
(*To get sigma bar, we should take the transpose too, as in Srednicki 35.19, i.e. switch a and adot*)
sigmaBar=Table[Table[Sum[epMetricUp[[a,b]]*epMetricUp[[adot,bdot]]*sigma[[i,b,bdot]],{b,1,2},{bdot,1,2}],{adot,1,2},{a,1,2}],{i,1,4}];


(*From Srednicki (60.7)*)
polPlus[k_,q_]:=Table[+(abl[q].sigmaBar[[\[Mu]]].sbr[k]/(Sqrt[2]abl[q].abr[k])),{\[Mu],1,4}]; (*Use + from Dixon instead of - from Srednicki. We get \[Epsilon]p.\[Epsilon]m=-1 then*)
polMinus[k_,q_]:=Table[-(sbl[q].sigma[[\[Mu]]].abr[k]/(Sqrt[2]sbl[q].sbr[k])),{\[Mu],1,4}];


(* ::Subsection::Closed:: *)
(*Generate and Set Vectors*)


SetVec[point_List]:=Module[{i},For[i=1,i<=Length[point],i++, VEC[point[[i,1]]]=point[[i,2]];];];


(*generate an on-shell d-dimensional vector, padded with n zeros*)
ClearAll[GenVectorN,GenVectorNC,GenVector];

GenVectorN[d_, n_, wp_:$MinPrecision]:=Module[{zeroPad,rand},
	zeroPad=ConstantArray[0,n];
	rand=Table[RandComp[wp],{i,1,d-n-1}];
	(*random whether incoming or outgoing*)
	Join[{Sign[RandomReal[{-1,1}]]Sqrt[Sum[rand[[i]]^2,{i,1,d-n-1}]]},Table[rand[[i]],{i,1,d-n-1}],zeroPad]
];

(*Complex*)
GenVectorNC[d_, n_, wp_:$MinPrecision]:=Module[{zeroPad,rand},
	zeroPad=ConstantArray[0,n];
	rand=Table[RandCompC[wp],{i,1,d-n-1}];
	(*random whether incoming or outgoing*)
	Join[{Sign[RandomReal[{-1,1}]]Sqrt[Sum[rand[[i]]^2,{i,1,d-n-1}]]},Table[rand[[i]],{i,1,d-n-1}],zeroPad]
];

(*Symbolic*)
GenVector[d_, n_, a_]:=Module[{zeroPad,aArray},
	zeroPad=ConstantArray[0,n];
	aArray=Table[a[i],{i,1,d-n-1}];
	(*random whether incoming or outgoing*)
	Join[{Sign[RandomReal[{-1,1}]]Sqrt[Sum[aArray[[i]]^2,{i,1,d-n-1}]]},Table[aArray[[i]],{i,1,d-n-1}],zeroPad]
];


(* ::Subsection::Closed:: *)
(*Generate Polarization Vectors for External States*)


(* The physical dimension is always 4, but d can be >= 4 *)
(* d = 6, dPhys = 4, means treat d-dPhys as internal -2\[Epsilon] dimensions *)
(*ClearAll[GenPol];*)
GenPol[helList_List] := Module[
{n, dPhys, d, pad, refMom, zeroPad, epPlus, epMinus, i},
	dPhys = 4;
	n = Length[helList];
	(* Check if the n momenta, vec[i], have already been generated *)
	If[ Length[Position[Head/@Table[VEC[i],{i,1,n}],VEC]] > 0,
		Print["Need to generate numerical momenta, VEC[i], before generating polarizations."];
		Return[];
	];
	d = Length[VEC[1]];
	pad = d - dPhys;
	refMom = Table[GenVectorNC[d,pad],{i,1,n}];
	zeroPad = ConstantArray[0,pad];

	(*pad 4-dim polarization vectors with zeros*)
	epPlus[p_,q_] := Join[polPlus[p,q],zeroPad];
	epMinus[p_,q_]:= Join[polMinus[p,q],zeroPad];	
	(*Generate the four external polarization vectors*)
	For[i = 1, i <= n, i++,
		Switch[helList[[i]],
			"p",
				POL[i] = epPlus[VEC[i],refMom[[i]]];,
			"m",
				POL[i] = epMinus[VEC[i],refMom[[i]]];,
			_,
				Print["Helicity list input needs to be ''p'' or ''m''"];
		];		
	];
];


(* The physical dimension is always 4, but d can be >= 4 *)
(* d = 6, dPhys = 4, means treat d-dPhys as internal -2\[Epsilon] dimensions *)
(*ClearAll[GenPol];*)
GenPol2[helList_List] := Module[
{n, d, dPhys, pad, refMom, zeroPad, epPlus, epMinus, i},
	dPhys = 4;
	n = Length[helList];
	(* Check if the n momenta, vec[i], have already been generated *)
	If[ Length[Position[Head/@Table[VEC[i],{i,1,n}],VEC]] > 0,
		Print["Need to generate numerical momenta, VEC[i], before generating polarizations."];
		Return[];
	];
	d = Length[VEC[1]];
	pad = d - dPhys;
	refMom = Table[GenVectorNC[d,pad],{i,1,n}];
	zeroPad = ConstantArray[0,pad];

	(*pad 4-dim polarization vectors with zeros*)
	epPlus[p_,q_] := Join[polPlus[p,q],zeroPad];
	epMinus[p_,q_]:= Join[polMinus[p,q],zeroPad];	
	(*Generate the four external polarization vectors*)
	For[i = 1, i <= n, i++,
		Switch[helList[[i]],
			"p",
				POL2[i] = epPlus[VEC[i],refMom[[i]]];,
			"m",
				POL2[i] = epMinus[VEC[i],refMom[[i]]];,
			_,
				Print["Helicity list input needs to be ''p'' or ''m''"];
		];		
	];
];


(* ::Subsection::Closed:: *)
(*Convert Symbolic to Numeric*)


(*Changed sign in front of NDot mud rules because we are in the mostly-minus convention*)


(*Convert dot products (ee, EE, ke, kE, kk, s) to numerics*)
(*For gravitons, we use the same helicity for e and E, i.e. the L and R copies*)
NDot[expr_]:=Module[{rep},
	rep={
		ee[i_,j_]:>Mink[POL[i],POL[j]],
		eE[i_,j_]:>Mink[POL[i],POL2[j]],
		EE[i_,j_]:>Mink[POL2[i],POL2[j]],
		ke[i_,j_]:>Mink[VEC[i],POL[j]],
		kE[i_,j_]:>Mink[VEC[i],POL2[j]],
		kk[i_,j_]:>Mink[VEC[i],VEC[j]],
		le[i_,j_]:>Mink[VEC[i],POL[j]],
		lE[i_,j_]:>Mink[VEC[i],POL2[j]],
		ll[i_,j_]:>Mink[VEC[i],VEC[j]],
		lk[i_,j_]:>Mink[VEC[i],VEC[j]],
		s[a__]:>Mink[Plus@@(VEC/@{a}),Plus@@(VEC/@{a})]
		(*
		(*Need and -1 in front since we are using mostly minus convention*)
		mud[i_Integer,j_Integer]:>-Mink[Join[ConstantArray[0,dPhys],VEC[i][[dPhys+1;;d]]],Join[ConstantArray[0,dPhys],VEC[j][[dPhys+1;;d]]]],
		mud[a_,b_]:>-Mink[Join[ConstantArray[0,dPhys],(a/.p->VEC)[[dPhys+1;;d]]],Join[ConstantArray[0,dPhys],(b/.p->VEC)[[dPhys+1;;d]]]],
		(*Dual vector*)
		lv[i_]:>-Det[Take[#,dPhys]&/@{Sign[i]VEC[Abs[i]],VEC[1],VEC[2],VEC[3]}],
		kv[i_]:>-Det[Take[#,dPhys]&/@{Sign[i]VEC[Abs[i]],VEC[1],VEC[2],VEC[3]}],
		lv[i_,j_,k_,l_]:>-Det[Take[#,dPhys]&/@{Sign[i]VEC[Abs[i]],Sign[j]VEC[Abs[j]],Sign[k]VEC[Abs[k]],Sign[l]VEC[Abs[l]]}],
		kv[i_,j_,k_,l_]:>-Det[Take[#,dPhys]&/@{Sign[i]VEC[Abs[i]],Sign[j]VEC[Abs[j]],Sign[k]VEC[Abs[k]],Sign[l]VEC[Abs[l]]}]
		*)
		};

	(expr/.R[a_]:>AngleBracket[a,a]) /. rep
];


(*Convert free indices to numerics*)
NFreeIndex[expr_]:=Module[{freeIndexRep},
	freeIndexRep=
		{\[Eta][i_,j_]:>g[i,j], 
		p[i_][j_]:>VEC[i][[j]],
		\[Epsilon][i_][j_]:>POL[i][[j]],
		\[Epsilon]2[i_][j_]:>POL2[i][[j]]};

	expr /. freeIndexRep
];


(* ::Subsection::Closed:: *)
(*Cut Numerics*)


(* VERY messy, but it works very well and very fast! *)


(* ::Subsubsection::Closed:: *)
(*Generate Cut Numerics*)
(*(Input: Contact vertex sturcture)*)


(*ClearAll[CutNumerics];*)
CutNumerics[verts_List,loopEdges_List,d_Integer,dPhys_Integer,wp_:$MinPrecision]:=Module[
{pad,extVerts,intVerts,contactVerts,edges,extEdges,intEdges,cutEdges,edgePos,indepExtEdges,indepEdges,depEdges,
momConsEqs,momConsRep,momList,mom,
onShellEdges,positionDuplicates,noOnShellCondEdges,dependence,
indepAppear,indepChoice,availIndepExtComp,availIndepLoopComp,indepWithMaxAvailIndepComp,
otherPositionsOfEdge,onShellIndepChoice,numIndepCompsNeeded,
a,vList,aExtReps,aLoopReps,aReps,
onShellEqs,eqsWithLoop,eqsWOLoop,eqsWOLoop2,eqsWithLoop2,solnWOLoop,solnWithLoop,soln,
solnWOLoopList,eqsWithLoopList,solnWithLoopList,solnList,numsList,solvevars,
ratPrec,analyticCutSoln,
vars,flag,counter,point,i},

	pad=d-dPhys;	

	extVerts = Cases[verts,a_ /; Length[a]==1];
	intVerts = Cases[verts,a_ /; Length[a]>=3];
	contactVerts = Cases[intVerts,a_ /; Length[a]>3];

	edges=Sort[DeleteDuplicates[Flatten[Abs[verts]]]];
	(*Since we don't expect every edge to appear in our contact graphs, this function will tell you the position of the edge number*)
	(*eg. edges={1,2,3,4,5,7,9,10,11}, edgePos[7]=6*)
	edgePos[edge_]:=Position[edges,edge][[1,1]];

	extEdges=Sort[Abs[Flatten[extVerts]]];
	intEdges=Sort[Delete[edges,List/@(edgePos/@extEdges)]];
	cutEdges=intEdges;
	
	(*the independent momenta are n-1 external momenta and the loop momenta*)
	indepExtEdges=Delete[extEdges,{-1}];
	indepEdges=Join[indepExtEdges,loopEdges];
	depEdges=Complement[edges,indepEdges];

	(*Momentum Conservation*)
	momConsEqs = #==0&/@Flatten[(Sum[(Sign[#[[i]]]mom[Abs[#[[i]]]]),{i,1,Length[#]}])&/@intVerts];
	momConsRep = Solve[momConsEqs,mom/@depEdges][[1]];
	momList=(mom/@edges)/.momConsRep;

	onShellEdges=Sort[Join[extEdges,cutEdges]];
	(*Find the positions of edges with the same momenta up to a sign - like cutting the intermediate props of internal bubbles*)
	(*We don't want to impose the same cut condition more than once*)
	(*This would allocate more free components than we need, leaving the system over-determined*)
	positionDuplicates[list_] := DeleteCases[GatherBy[Range@Length[list], Sort[Expand[{list[[#]],-list[[#]]}]] &],{_}];
	onShellEdges=Complement[onShellEdges,Flatten[Delete[#,{1}]&/@positionDuplicates[momList]]];
	noOnShellCondEdges=Complement[edges,onShellEdges];
	dependence=Cases[{#},mom[i_]:>i,\[Infinity]]&/@momList;

	(*The momenta that depend on the independent momenta - can depend on more than one*)
	(*Figure out how many components of each independent momenta to randomly generate*)
	(*Need one component of an indep mom to put a vector on-shell*)
	(*Unless we do max cuts, we won't put everything on-shell*)
	indepAppear=(mom[edges[[#]]]&/@DeleteDuplicates[Position[momList,mom[#],\[Infinity]][[All,1]]])&/@edges;
	(*For each edge, choose the independent momentum that we will take a component from to make the edge on-shell*)
	(*Need to make sure we don't choose an indep mom more than d (or dPhys) times*)
	(*Give preference to loop momenta so we can solve for ext mom separately*)
	(*If it depends on multiple indep moms, choose the indep mom that has the least number of edges dependent on it*)
	indepChoice={};
	For[i=1,i<=Length[edges],i++,
		(*Number of available components for the independent momenta after demanding all of its dependent momenta on-shell*)
		availIndepExtComp = (If[MemberQ[extEdges,#],dPhys-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		availIndepLoopComp = (If[MemberQ[loopEdges,#],d-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		(*Use loop momenta before external momenta*)
		indepWithMaxAvailIndepComp = If[Length[DeleteCases[availIndepLoopComp,Null]]>0,
										dependence[[i]][[Position[availIndepLoopComp,Max[DeleteCases[availIndepLoopComp,Null]]][[1,1]]]],
										dependence[[i]][[Position[availIndepExtComp,Max[DeleteCases[availIndepExtComp,Null]]][[1,1]]]]
									];
		indepChoice = AppendTo[indepChoice,indepWithMaxAvailIndepComp];
		otherPositionsOfEdge = DeleteCases[Position[indepAppear,mom[edges[[i]]],\[Infinity]],{edgePos[indepWithMaxAvailIndepComp],__},\[Infinity]];
		indepAppear = Delete[indepAppear,otherPositionsOfEdge];
	];
	(*only need indep mom component for ext legs and cut legs*)
	onShellIndepChoice=ReplacePart[indepChoice,(#->0)&/@noOnShellCondEdges];
	(*The number of independent components needed to put desired edges on-shell*)
	numIndepCompsNeeded=Count[onShellIndepChoice,#]&/@indepEdges;

	vList=Join[
			Join[Table[a[#,i],{i,1,dPhys}],ConstantArray[0,pad]]&/@extEdges,
			Table[a[#,i],{i,1,d}]&/@intEdges
			];
	vList=((mom/@edges)/.momConsRep)/.mom[i_]:>vList[[edgePos[i]]];

	(*indepEdge = Join[indepExtEdges,loopEdges] in that order*)
	aExtReps=Flatten[Table[a[indepEdges[[i]],j]->RandCompC[wp],{i,1,Length[indepExtEdges]},{j,1+numIndepCompsNeeded[[i]],dPhys}]];
	(*Make sure to have one undetermined coefficient in the d-dPhys dimensions to put cut loop edge on shell*)
	(*Otherwise, can lead to inconsistent equations. See max cuts of triangle diagram. See Simon Badger's paper. No problem with box or bubble.*)
	aLoopReps=Flatten[Table[a[indepEdges[[i]],j]->RandCompC[wp],{i,Length[indepExtEdges]+1,Length[indepEdges]},{j,(*1+*)numIndepCompsNeeded[[i]],d-1}]];
	aReps=Join[aExtReps,aLoopReps];

	vList=vList/.aReps;

	ratPrec=rationalizePrecision;(*Precision that we rationalize with. Doesn't work if too big for some reason*)
	(*On-Shell Condition*)
	onShellEqs=(AngleBracket[#,#]==0&/@(momList[[edgePos/@onShellEdges]]/.mom->p));
	eqsWithLoop=Sort[DeleteDuplicates[Flatten[Position[onShellEqs/.Equal->Plus,kk[#,_]][[All,1]]&/@loopEdges]]];
	eqsWOLoop=Complement[Range[Length[onShellEqs]],eqsWithLoop];

	(*Solve equation that only depend on ext mom first*)
	eqsWOLoop2=Solve[onShellEqs[[eqsWOLoop]]][[1]];
	eqsWOLoop=(eqsWOLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]];

(*
	solnWOLoopList=NSolve[Rationalize[eqsWOLoop,10^-ratPrec],WorkingPrecision->wp];//Timing
*)

	(*****New Solution Method*******)
	vars=Variables[eqsWOLoop/.Equal[x_,y_]:>x-y];
	flag=False;
	counter=0;
	While[And[flag===False,counter<100],
		solnWOLoop=FindRoot[eqsWOLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
		flag=And@@(#==0&/@Chop[(eqsWOLoop/.Equal[x_,y_]:>x-y)/.solnWOLoop]);
		counter=counter+1;
	];
	(********************************)
	(*Solve the other on-shell equation*)
	eqsWithLoop2=Solve[onShellEqs[[eqsWithLoop]]/.eqsWOLoop2][[1]];
	eqsWithLoop=(Expand[((eqsWithLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]])//.solnWOLoop]//Chop);

(*
	solnWithLoopList=NSolve[Rationalize[#,10^-ratPrec],WorkingPrecision->wp]&/@eqsWithLoopList;//Timing
*)

	If[eqsWithLoop==={},
		solnWithLoop={};,

		(*****New Solution Method*******)
		vars=Variables[eqsWithLoop/.Equal[x_,y_]:>x-y];
		flag=False;
		counter=0;
			While[And[flag===False,counter<100],
				solnWithLoop=FindRoot[eqsWithLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
				flag=And@@(#==0&/@Chop[(eqsWithLoop/.Equal[x_,y_]:>x-y)/.solnWithLoop]);
				counter=counter+1;
			];
		(********************************)
	];
	
	(*FindRoot only finds one solution, but keep in mind to following??*)
	(*NOTE: The solution choices become important for max cuts.*)
	(*For instance, a three vertex (1,2,9) needs either <12> = <19> = <29> = 0 or square brackets vanishing*)
	(*If square brackets vanish, then A(1+,2+,9-)=0*)
	(*solnList=Flatten[Table[Join[solnWOLoopList[[i]],solnWithLoopList[[i,j]]],{i,1,Length[solnWOLoopList]},{j,1,Length[solnWithLoopList[[i]]]}],1];*)
	soln=Join[solnWOLoop,solnWithLoop];

	point=(vList/.soln);
	point=Table[{edges[[i]],point[[i]]},{i,1,Length[edges]}];

	point
];


(* ::Subsubsection::Closed:: *)
(*Generate Cut Numerics (Soft)*)
(*(Always make 2st external leg soft)*)


(*ClearAll[CutNumericsSoft];*)
CutNumericsSoft[verts_List,loopEdges_List,d_Integer,dPhys_Integer,softScaling_,wp_:$MinPrecision]:=Module[

{pad,extVerts,intVerts,contactVerts,edges,extEdges,intEdges,cutEdges,edgePos,indepExtEdges,indepEdges,depEdges,
momConsEqs,momConsRep,momList,mom,
onShellEdges,positionDuplicates,noOnShellCondEdges,dependence,
indepAppear,indepChoice,availIndepExtComp,availIndepLoopComp,indepWithMaxAvailIndepComp,
otherPositionsOfEdge,onShellIndepChoice,numIndepCompsNeeded,
a,vList,aExtReps,aLoopReps,aReps,
onShellEqs,eqsWithLoop,eqsWOLoop,eqsWOLoop2,eqsWithLoop2,solnWOLoop,solnWithLoop,soln,
solnWOLoopList,eqsWithLoopList,solnWithLoopList,solnList,numsList,solvevars,
ratPrec,analyticCutSoln,
vars,flag,counter,point,i},

	pad=d-dPhys;	

	extVerts = Cases[verts,a_ /; Length[a]==1];
	intVerts = Cases[verts,a_ /; Length[a]>=3];
	contactVerts = Cases[intVerts,a_ /; Length[a]>3];

	edges=Sort[DeleteDuplicates[Flatten[Abs[verts]]]];
	(*Since we don't expect every edge to appear in our contact graphs, this function will tell you the position of the edge number*)
	(*eg. edges={1,2,3,4,5,7,9,10,11}, edgePos[7]=6*)
	edgePos[edge_]:=Position[edges,edge][[1,1]];

	extEdges=Sort[Abs[Flatten[extVerts]]];
	intEdges=Sort[Delete[edges,List/@(edgePos/@extEdges)]];
	cutEdges=intEdges;
	
	(*the independent momenta are n-1 external momenta and the loop momenta*)
	indepExtEdges=Delete[extEdges,{-1}];
	indepEdges=Join[indepExtEdges,loopEdges];
	depEdges=Complement[edges,indepEdges];

	(*Momentum Conservation*)
	momConsEqs = #==0&/@Flatten[(Sum[(Sign[#[[i]]]mom[Abs[#[[i]]]]),{i,1,Length[#]}])&/@intVerts];
	momConsRep = Solve[momConsEqs,mom/@depEdges][[1]];
	momList=(mom/@edges)/.momConsRep;

	onShellEdges=Sort[Join[extEdges,cutEdges]];
	(*Find the positions of edges with the same momenta up to a sign - like cutting the intermediate props of internal bubbles*)
	(*We don't want to impose the same cut condition more than once*)
	(*This would allocate more free components than we need, leaving the system over-determined*)
	positionDuplicates[list_] := DeleteCases[GatherBy[Range@Length[list], Sort[Expand[{list[[#]],-list[[#]]}]] &],{_}];
	onShellEdges=Complement[onShellEdges,Flatten[Delete[#,{1}]&/@positionDuplicates[momList]]];
	noOnShellCondEdges=Complement[edges,onShellEdges];
	dependence=Cases[{#},mom[i_]:>i,\[Infinity]]&/@momList;

	(*The momenta that depend on the independent momenta - can depend on more than one*)
	(*Figure out how many components of each independent momenta to randomly generate*)
	(*Need one component of an indep mom to put a vector on-shell*)
	(*Unless we do max cuts, we won't put everything on-shell*)
	indepAppear=(mom[edges[[#]]]&/@DeleteDuplicates[Position[momList,mom[#],\[Infinity]][[All,1]]])&/@edges;
	(*For each edge, choose the independent momentum that we will take a component from to make the edge on-shell*)
	(*Need to make sure we don't choose an indep mom more than d (or dPhys) times*)
	(*Give preference to loop momenta so we can solve for ext mom separately*)
	(*If it depends on multiple indep moms, choose the indep mom that has the least number of edges dependent on it*)
	indepChoice={};
	For[i=1,i<=Length[edges],i++,
		(*Number of available components for the independent momenta after demanding all of its dependent momenta on-shell*)
		availIndepExtComp = (If[MemberQ[extEdges,#],dPhys-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		availIndepLoopComp = (If[MemberQ[loopEdges,#],d-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		(*Use loop momenta before external momenta*)
		indepWithMaxAvailIndepComp = If[Length[DeleteCases[availIndepLoopComp,Null]]>0,
										dependence[[i]][[Position[availIndepLoopComp,Max[DeleteCases[availIndepLoopComp,Null]]][[1,1]]]],
										dependence[[i]][[Position[availIndepExtComp,Max[DeleteCases[availIndepExtComp,Null]]][[1,1]]]]
									];
		indepChoice = AppendTo[indepChoice,indepWithMaxAvailIndepComp];
		otherPositionsOfEdge = DeleteCases[Position[indepAppear,mom[edges[[i]]],\[Infinity]],{edgePos[indepWithMaxAvailIndepComp],__},\[Infinity]];
		indepAppear = Delete[indepAppear,otherPositionsOfEdge];
	];
	(*only need indep mom component for ext legs and cut legs*)
	onShellIndepChoice=ReplacePart[indepChoice,(#->0)&/@noOnShellCondEdges];
	(*The number of independent components needed to put desired edges on-shell*)
	numIndepCompsNeeded=Count[onShellIndepChoice,#]&/@indepEdges;

	vList=Join[
			Join[Table[a[#,i],{i,1,dPhys}],ConstantArray[0,pad]]&/@extEdges,
			Table[a[#,i],{i,1,d}]&/@intEdges
			];
	vList=((mom/@edges)/.momConsRep)/.mom[i_]:>vList[[edgePos[i]]];

	(*indepEdge = Join[indepExtEdges,loopEdges] in that order*)
	aExtReps=Flatten[Table[a[indepEdges[[i]],j]->(If[i===2,softScaling,1]*RandCompC[wp]),{i,1,Length[indepExtEdges]},{j,1+numIndepCompsNeeded[[i]],dPhys}]];
	(*Make sure to have one undetermined coefficient in the d-dPhys dimensions to put cut loop edge on shell*)
	(*Otherwise, can lead to inconsistent equations. See max cuts of triangle diagram. See Simon Badger's paper. No problem with box or bubble.*)
	aLoopReps=Flatten[Table[a[indepEdges[[i]],j]->RandCompC[wp],{i,Length[indepExtEdges]+1,Length[indepEdges]},{j,(*1+*)numIndepCompsNeeded[[i]],d-1}]];
	aReps=Join[aExtReps,aLoopReps];

	vList=vList/.aReps;

	ratPrec=rationalizePrecision;(*Precision that we rationalize with. Doesn't work if too big for some reason*)
	(*On-Shell Condition*)
	onShellEqs=(AngleBracket[#,#]==0&/@(momList[[edgePos/@onShellEdges]]/.mom->p));
	eqsWithLoop=Sort[DeleteDuplicates[Flatten[Position[onShellEqs/.Equal->Plus,kk[#,_]][[All,1]]&/@loopEdges]]];
	eqsWOLoop=Complement[Range[Length[onShellEqs]],eqsWithLoop];

	(*Solve equation that only depend on ext mom first*)
	eqsWOLoop2=Solve[onShellEqs[[eqsWOLoop]]][[1]];
	eqsWOLoop=(eqsWOLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]];

(*
	solnWOLoopList=NSolve[Rationalize[eqsWOLoop,10^-ratPrec],WorkingPrecision->wp];//Timing
*)

	(*****New Solution Method*******)
	vars=Variables[eqsWOLoop/.Equal[x_,y_]:>x-y];
	flag=False;
	counter=0;
	While[And[flag===False,counter<100],
		solnWOLoop=FindRoot[eqsWOLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
		flag=And@@(#==0&/@Chop[(eqsWOLoop/.Equal[x_,y_]:>x-y)/.solnWOLoop]);
		counter=counter+1;
	];
	(********************************)
	(*Solve the other on-shell equation*)
	eqsWithLoop2=Solve[onShellEqs[[eqsWithLoop]]/.eqsWOLoop2][[1]];
	eqsWithLoop=(Expand[((eqsWithLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]])//.solnWOLoop]//Chop);

(*
	solnWithLoopList=NSolve[Rationalize[#,10^-ratPrec],WorkingPrecision->wp]&/@eqsWithLoopList;//Timing
*)

	If[eqsWithLoop==={},
		solnWithLoop={};,

		(*****New Solution Method*******)
		vars=Variables[eqsWithLoop/.Equal[x_,y_]:>x-y];
		flag=False;
		counter=0;
			While[And[flag===False,counter<100],
				solnWithLoop=FindRoot[eqsWithLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
				flag=And@@(#==0&/@Chop[(eqsWithLoop/.Equal[x_,y_]:>x-y)/.solnWithLoop]);
				counter=counter+1;
			];
		(********************************)
	];
	
	(*FindRoot only finds one solution, but keep in mind to following??*)
	(*NOTE: The solution choices become important for max cuts.*)
	(*For instance, a three vertex (1,2,9) needs either <12> = <19> = <29> = 0 or square brackets vanishing*)
	(*If square brackets vanish, then A(1+,2+,9-)=0*)
	(*solnList=Flatten[Table[Join[solnWOLoopList[[i]],solnWithLoopList[[i,j]]],{i,1,Length[solnWOLoopList]},{j,1,Length[solnWithLoopList[[i]]]}],1];*)
	soln=Join[solnWOLoop,solnWithLoop];

	point=(vList/.soln);
	point=Table[{edges[[i]],point[[i]]},{i,1,Length[edges]}];

	point
];


(* ::Subsubsection::Closed:: *)
(*Generate Cut Numerics (Mu Scaling) [Needs to be updated to match cut numerics]*)
(*(Input: Contact vertex sturcture)*)


ClearAll[CutNumericsMuScaling];
CutNumericsMuScaling[verts_List,loopEdges_List,d_Integer,dPhys_Integer,muScaling_,wp_:$MinPrecision]:=Module[
{pad,extVerts,intVerts,contactVerts,edges,extEdges,intEdges,cutEdges,edgePos,indepExtEdges,indepEdges,depEdges,
momConsEqs,momConsRep,momList,mom,
onShellEdges,positionDuplicates,noOnShellCondEdges,dependence,
indepAppear,indepChoice,availIndepExtComp,availIndepLoopComp,indepWithMaxAvailIndepComp,
otherPositionsOfEdge,onShellIndepChoice,numIndepCompsNeeded,
a,vList,aExtReps,aLoopReps,aReps,
onShellEqs,eqsWithLoop,eqsWOLoop,eqsWOLoop2,eqsWithLoop2,solnWOLoop,solnWithLoop,soln,
solnWOLoopList,eqsWithLoopList,solnWithLoopList,solnList,numsList,solvevars,
ratPrec,analyticCutSoln,
vars,flag,counter,point,i},

	pad=d-dPhys;	

	extVerts = Cases[verts,a_ /; Length[a]==1];
	intVerts = Cases[verts,a_ /; Length[a]>=3];
	contactVerts = Cases[intVerts,a_ /; Length[a]>3];

	edges=Sort[DeleteDuplicates[Flatten[Abs[verts]]]];
	extEdges=Sort[Abs[Flatten[extVerts]]];
	intEdges=Sort[Delete[edges,List/@extEdges]];
	cutEdges=intEdges;
	
	(*Since we don't expect every edge to appear in our contact graphs, this function will tell you the position of the edge number*)
	(*eg. edges={1,2,3,4,5,7,9,10,11}, edgePos[7]=6*)
	edgePos[edge_]:=Position[edges,edge][[1,1]];

	(*the independent momenta are n-1 external momenta and the loop momenta*)
	indepExtEdges=Delete[extEdges,{-1}];
	indepEdges=Join[indepExtEdges,loopEdges];
	depEdges=Complement[edges,indepEdges];

	(*Momentum Conservation*)
	momConsEqs = #==0&/@Flatten[(Sum[(Sign[#[[i]]]mom[Abs[#[[i]]]]),{i,1,Length[#]}])&/@intVerts];
	momConsRep = Solve[momConsEqs,mom/@depEdges][[1]];
	momList=(mom/@edges)/.momConsRep;

	onShellEdges=Sort[Join[extEdges,cutEdges]];
	(*Find the positions of edges with the same momenta up to a sign - like cutting the intermediate props of internal bubbles*)
	(*We don't want to impose the same cut condition more than once*)
	(*This would allocate more free components than we need, leaving the system over-determined*)
	positionDuplicates[list_] := DeleteCases[GatherBy[Range@Length[list], Sort[Expand[{list[[#]],-list[[#]]}]] &],{_}];
	onShellEdges=Complement[onShellEdges,Flatten[Delete[#,{1}]&/@positionDuplicates[momList]]];
	noOnShellCondEdges=Complement[edges,onShellEdges];
	dependence=Cases[{#},mom[i_]:>i,\[Infinity]]&/@momList;

	(*The momenta that depend on the independent momenta - can depend on more than one*)
	(*Figure out how many components of each independent momenta to randomly generate*)
	(*Need one component of an indep mom to put a vector on-shell*)
	(*Unless we do max cuts, we won't put everything on-shell*)
	indepAppear=(mom[edges[[#]]]&/@DeleteDuplicates[Position[momList,mom[#],\[Infinity]][[All,1]]])&/@edges;
	(*For each edge, choose the independent momentum that we will take a component from to make the edge on-shell*)
	(*Need to make sure we don't choose an indep mom more than d (or dPhys) times*)
	(*Give preference to loop momenta so we can solve for ext mom separately*)
	(*If it depends on multiple indep moms, choose the indep mom that has the least number of edges dependent on it*)
	indepChoice={};
	For[i=1,i<=Length[edges],i++,
		(*Number of available components for the independent momenta after demanding all of its dependent momenta on-shell*)
		availIndepExtComp = (If[MemberQ[extEdges,#],dPhys-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		availIndepLoopComp = (If[MemberQ[loopEdges,#],d-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		(*Use loop momenta before external momenta*)
		indepWithMaxAvailIndepComp = If[Length[DeleteCases[availIndepLoopComp,Null]]>0,
										dependence[[i]][[Position[availIndepLoopComp,Max[DeleteCases[availIndepLoopComp,Null]]][[1,1]]]],
										dependence[[i]][[Position[availIndepExtComp,Max[DeleteCases[availIndepExtComp,Null]]][[1,1]]]]
									];
		indepChoice = AppendTo[indepChoice,indepWithMaxAvailIndepComp];
		otherPositionsOfEdge = DeleteCases[Position[indepAppear,mom[edges[[i]]],\[Infinity]],{edgePos[indepWithMaxAvailIndepComp],__},\[Infinity]];
		indepAppear = Delete[indepAppear,otherPositionsOfEdge];
	];
	(*only need indep mom component for ext legs and cut legs*)
	onShellIndepChoice=ReplacePart[indepChoice,(#->0)&/@noOnShellCondEdges];
	(*The number of independent components needed to put desired edges on-shell*)
	numIndepCompsNeeded=Count[onShellIndepChoice,#]&/@indepEdges;

	vList=Join[
			Join[Table[a[#,i],{i,1,dPhys}],ConstantArray[0,pad]]&/@extEdges,
			Table[a[#,i],{i,1,d}]&/@intEdges
			];
	vList=((mom/@edges)/.momConsRep)/.mom[i_]:>vList[[edgePos[i]]];

	(*indepEdge = Join[indepExtEdges,loopEdges] in that order*)
	aExtReps=Flatten[Table[a[indepEdges[[i]],j]->RandCompC[wp],{i,1,Length[indepExtEdges]},{j,1+numIndepCompsNeeded[[i]],dPhys}]];
	
	(*Make sure to have one undetermined coefficient in the d-dPhys dimensions to put cut loop edge on shell*)
	(*Otherwise, can lead to inconsistent equations. See max cuts of triangle diagram. See Simon Badger's paper. No problem with box or bubble.*)
	(*Actually, we want to scalar the mu's (extra dim pieces) here, so we'll try our luck*)
	(*Make random components small*)
	aLoopReps=Flatten[Table[a[indepEdges[[i]],j]->(If[j>dPhys,muScaling,1]*RandCompC[wp]),{i,Length[indepExtEdges]+1,Length[indepEdges]},{j,1+(numIndepCompsNeeded[[i]]),d}]];
	aReps=Join[aExtReps,aLoopReps];

	vList=vList/.aReps;

	(*ratPrec=rationalizePrecision;(*Precision that we rationalize with. Doesn't work if too big for some reason*)*)
	(*On-Shell Condition*)
	onShellEqs=(AngleBracket[#,#]==0&/@(momList[[edgePos/@onShellEdges]]/.mom->p));
	eqsWithLoop=Sort[DeleteDuplicates[Flatten[Position[onShellEqs/.Equal->Plus,kk[#,_]][[All,1]]&/@loopEdges]]];
	eqsWOLoop=Complement[Range[Length[onShellEqs]],eqsWithLoop];

	(*Solve equation that only depend on ext mom first*)
	eqsWOLoop2=Solve[onShellEqs[[eqsWOLoop]]][[1]];
	eqsWOLoop=(eqsWOLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]];

(*
	solnWOLoopList=NSolve[Rationalize[eqsWOLoop,10^-ratPrec],WorkingPrecision->wp];//Timing
*)

	(*****New Solution Method*******)
	vars=Variables[eqsWOLoop/.Equal[x_,y_]:>x-y];
	flag=False;
	counter=0;
	While[And[flag===False,counter<100],
		solnWOLoop=FindRoot[eqsWOLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
		flag=And@@(#==0&/@Chop[(eqsWOLoop/.Equal[x_,y_]:>x-y)/.solnWOLoop]);
		counter=counter+1;
	];
	(********************************)
	(*Solve the other on-shell equation*)
	eqsWithLoop2=Solve[onShellEqs[[eqsWithLoop]]/.eqsWOLoop2][[1]];
	eqsWithLoop=(Expand[((eqsWithLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]])//.solnWOLoop]//Chop);

(*
	solnWithLoopList=NSolve[Rationalize[#,10^-ratPrec],WorkingPrecision->wp]&/@eqsWithLoopList;//Timing
*)

	If[eqsWithLoop==={},
		solnWithLoop={};,

		(*****New Solution Method*******)
		vars=Variables[eqsWithLoop/.Equal[x_,y_]:>x-y];
		flag=False;
		counter=0;
			While[And[flag===False,counter<100],
				solnWithLoop=FindRoot[eqsWithLoop,{#,Random[Complex,{0,1000}]}&/@vars,WorkingPrecision->wp];
				flag=And@@(#==0&/@Chop[(eqsWithLoop/.Equal[x_,y_]:>x-y)/.solnWithLoop]);
				counter=counter+1;
			];
		(********************************)
	];
	
	(*FindRoot only finds one solution, but keep in mind to following??*)
	(*NOTE: The solution choices become important for max cuts.*)
	(*For instance, a three vertex (1,2,9) needs either <12> = <19> = <29> = 0 or square brackets vanishing*)
	(*If square brackets vanish, then A(1+,2+,9-)=0*)
	(*solnList=Flatten[Table[Join[solnWOLoopList[[i]],solnWithLoopList[[i,j]]],{i,1,Length[solnWOLoopList]},{j,1,Length[solnWithLoopList[[i]]]}],1];*)
	soln=Join[solnWOLoop,solnWithLoop];

	(*Hard-coded for 4-points*)
	analyticCutSoln=Solve[onShellEqs][[1]];
	solvevars=Flatten[Join[{kk[1,3]},Table[kk[i,i],{i,1,3}],Table[kk[loopEdges[[i]],j],{i,1,Length[loopEdges]},{j,1,3}],Table[kk[loopEdges[[i]],loopEdges[[j]]],{i,1,Length[loopEdges]},{j,i,Length[loopEdges]}]]];
	analyticCutSoln=ToRules[Reduce[analyticCutSoln/.Rule->Equal,solvevars]];

	point=(vList/.soln);
	point=Table[{edges[[i]],point[[i]]},{i,1,Length[edges]}];

	{point,momConsRep/.mom->p,analyticCutSoln}

];


(*
verts={{-1},{-2},{-3},{-4},{1,2,5,6},{3,4,-6,-5}};
d=8;
dPhys=4;
muScaling=10^-12;

{point,momcons,cutsoln}=CutNumericsNewMuScaling[verts,Quiet[GetLoopEdges[verts]],d,dPhys,muScaling];//Timing

((Plus@@((mom/@(#))/.mom[i_]:>Sign[i]*mom[Abs[i]]))/.mom[i_]:>point[[i,2]])&/@Select[verts,Length[#]>1&]
Sum[point[[i,2]],{i,1,4}]
Table[Mink[point[[i,2]],point[[i,2]]],{i,1,Length[point]}]
*)


Clear[RandCompR]
RandCompR[wp_:$MinPrecision]:=Sign[RandomReal[{-1,1}]]*RandomReal[{-randSize,randSize},WorkingPrecision->wp];

Clear[CutNumericsRealMuScaling];
CutNumericsRealMuScaling[verts_List,loopEdges_List,d_Integer,dPhys_Integer,muScaling_,wp_:$MinPrecision]:=Module[
{pad,extVerts,intVerts,contactVerts,edges,extEdges,intEdges,cutEdges,edgePos,indepExtEdges,indepEdges,depEdges,
momConsEqs,momConsRep,momList,mom,
onShellEdges,positionDuplicates,noOnShellCondEdges,dependence,
indepAppear,indepChoice,availIndepExtComp,availIndepLoopComp,indepWithMaxAvailIndepComp,
otherPositionsOfEdge,onShellIndepChoice,numIndepCompsNeeded,
a,vList,aExtReps,aLoopReps,aReps,
onShellEqs,eqsWithLoop,eqsWOLoop,eqsWOLoop2,eqsWithLoop2,solnWOLoop,solnWithLoop,soln,
solnWOLoopList,eqsWithLoopList,solnWithLoopList,solnList,numsList,solvevars,
ratPrec,analyticCutSoln,
vars,flag,counter,point},

	pad=d-dPhys;	

	extVerts = Cases[verts,a_ /; Length[a]==1];
	intVerts = Cases[verts,a_ /; Length[a]>=3];
	contactVerts = Cases[intVerts,a_ /; Length[a]>3];

	edges=Sort[DeleteDuplicates[Flatten[Abs[verts]]]];
	extEdges=Sort[Abs[Flatten[extVerts]]];
	intEdges=Sort[Delete[edges,List/@extEdges]];
	cutEdges=intEdges;
	
	(*Since we don't expect every edge to appear in our contact graphs, this function will tell you the position of the edge number*)
	(*eg. edges={1,2,3,4,5,7,9,10,11}, edgePos[7]=6*)
	edgePos[edge_]:=Position[edges,edge][[1,1]];

	(*the independent momenta are n-1 external momenta and the loop momenta*)
	indepExtEdges=Delete[extEdges,{-1}];
	indepEdges=Join[indepExtEdges,loopEdges];
	depEdges=Complement[edges,indepEdges];

	(*Momentum Conservation*)
	momConsEqs = #==0&/@Flatten[(Sum[(Sign[#[[i]]]mom[Abs[#[[i]]]]),{i,1,Length[#]}])&/@intVerts];
	momConsRep = Solve[momConsEqs,mom/@depEdges][[1]];
	momList=(mom/@edges)/.momConsRep;

	onShellEdges=Sort[Join[extEdges,cutEdges]];
	(*Find the positions of edges with the same momenta up to a sign - like cutting the intermediate props of internal bubbles*)
	(*We don't want to impose the same cut condition more than once*)
	(*This would allocate more free components than we need, leaving the system over-determined*)
	positionDuplicates[list_] := DeleteCases[GatherBy[Range@Length[list], Sort[Expand[{list[[#]],-list[[#]]}]] &],{_}];
	onShellEdges=Complement[onShellEdges,Flatten[Delete[#,{1}]&/@positionDuplicates[momList]]];
	noOnShellCondEdges=Complement[edges,onShellEdges];
	dependence=Cases[{#},mom[i_]:>i,\[Infinity]]&/@momList;

	(*The momenta that depend on the independent momenta - can depend on more than one*)
	(*Figure out how many components of each independent momenta to randomly generate*)
	(*Need one component of an indep mom to put a vector on-shell*)
	(*Unless we do max cuts, we won't put everything on-shell*)
	indepAppear=(mom[edges[[#]]]&/@DeleteDuplicates[Position[momList,mom[#],\[Infinity]][[All,1]]])&/@edges;
	(*For each edge, choose the independent momentum that we will take a component from to make the edge on-shell*)
	(*Need to make sure we don't choose an indep mom more than d (or dPhys) times*)
	(*Give preference to loop momenta so we can solve for ext mom separately*)
	(*If it depends on multiple indep moms, choose the indep mom that has the least number of edges dependent on it*)
	indepChoice={};
	For[i=1,i<=Length[edges],i++,
		(*Number of available components for the independent momenta after demanding all of its dependent momenta on-shell*)
		availIndepExtComp = (If[MemberQ[extEdges,#],dPhys-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		availIndepLoopComp = (If[MemberQ[loopEdges,#],d-Length[indepAppear[[edgePos[#]]]]]
						)&/@(dependence[[i]]);
		(*Use loop momenta before external momenta*)
		indepWithMaxAvailIndepComp = If[Length[DeleteCases[availIndepLoopComp,Null]]>0,
										dependence[[i]][[Position[availIndepLoopComp,Max[DeleteCases[availIndepLoopComp,Null]]][[1,1]]]],
										dependence[[i]][[Position[availIndepExtComp,Max[DeleteCases[availIndepExtComp,Null]]][[1,1]]]]
									];
		indepChoice = AppendTo[indepChoice,indepWithMaxAvailIndepComp];
		otherPositionsOfEdge = DeleteCases[Position[indepAppear,mom[edges[[i]]],\[Infinity]],{edgePos[indepWithMaxAvailIndepComp],__},\[Infinity]];
		indepAppear = Delete[indepAppear,otherPositionsOfEdge];
	];
	(*only need indep mom component for ext legs and cut legs*)
	onShellIndepChoice=ReplacePart[indepChoice,(#->0)&/@noOnShellCondEdges];
	(*The number of independent components needed to put desired edges on-shell*)
	numIndepCompsNeeded=Count[onShellIndepChoice,#]&/@indepEdges;

	vList=Join[
			Join[Table[a[#,i],{i,1,dPhys}],ConstantArray[0,pad]]&/@extEdges,
			Table[a[#,i],{i,1,d}]&/@intEdges
			];
	vList=((mom/@edges)/.momConsRep)/.mom[i_]:>vList[[edgePos[i]]];

	(*indepEdge = Join[indepExtEdges,loopEdges] in that order*)
	aExtReps=Flatten[Table[a[indepEdges[[i]],j]->RandCompR[wp],{i,1,Length[indepExtEdges]},{j,1+numIndepCompsNeeded[[i]],dPhys}]];
	
	(*Make sure to have one undetermined coefficient in the d-dPhys dimensions to put cut loop edge on shell*)
	(*Otherwise, can lead to inconsistent equations. See max cuts of triangle diagram. See Simon Badger's paper. No problem with box or bubble.*)
	(*Actually, we want to scalar the mu's (extra dim pieces) here, so we'll try our luck*)
	(*Make random components small*)
	aLoopReps=Flatten[Table[a[indepEdges[[i]],j]->(If[j>dPhys,muScaling,1]*RandCompR[wp]),{i,Length[indepExtEdges]+1,Length[indepEdges]},{j,1+(numIndepCompsNeeded[[i]]),d}]];
	aReps=Join[aExtReps,aLoopReps];

	vList=vList/.aReps;

	(*ratPrec=rationalizePrecision;(*Precision that we rationalize with. Doesn't work if too big for some reason*)*)
	(*On-Shell Condition*)
	onShellEqs=(AngleBracket[#,#]==0&/@(momList[[edgePos/@onShellEdges]]/.mom->p));
	eqsWithLoop=Sort[DeleteDuplicates[Flatten[Position[onShellEqs/.Equal->Plus,kk[#,_]][[All,1]]&/@loopEdges]]];
	eqsWOLoop=Complement[Range[Length[onShellEqs]],eqsWithLoop];

	(*Solve equation that only depend on ext mom first*)
	eqsWOLoop2=Solve[onShellEqs[[eqsWOLoop]]][[1]];
	eqsWOLoop=(eqsWOLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]];

(*
	solnWOLoopList=NSolve[Rationalize[eqsWOLoop,10^-ratPrec],WorkingPrecision->wp];//Timing
*)

	(*****New Solution Method*******)
	vars=Variables[eqsWOLoop/.Equal[x_,y_]:>x-y];
	flag=False;
	counter=0;
	While[And[flag===False,counter<100],
		solnWOLoop=FindRoot[eqsWOLoop,{#,Random[Real,{-1000,1000}]}&/@vars,WorkingPrecision->wp];
		flag=And@@(#==0&/@Chop[(eqsWOLoop/.Equal[x_,y_]:>x-y)/.solnWOLoop]);
		counter=counter+1;
	];
	(********************************)
	(*Solve the other on-shell equation*)
	eqsWithLoop2=Solve[onShellEqs[[eqsWithLoop]]/.eqsWOLoop2][[1]];
	eqsWithLoop=(Expand[((eqsWithLoop2/.Rule->Equal)/.kk[i_,j_]:>Mink[vList[[edgePos[i]]],vList[[edgePos[j]]]])//.solnWOLoop]//Chop);

(*
	solnWithLoopList=NSolve[Rationalize[#,10^-ratPrec],WorkingPrecision->wp]&/@eqsWithLoopList;//Timing
*)

	If[eqsWithLoop==={},
		solnWithLoop={};,

		(*****New Solution Method*******)
		vars=Variables[eqsWithLoop/.Equal[x_,y_]:>x-y];
		flag=False;
		counter=0;
			While[And[flag===False,counter<100],
				solnWithLoop=FindRoot[eqsWithLoop,{#,Random[Real,{-1000,1000}]}&/@vars,WorkingPrecision->wp];
				flag=And@@(#==0&/@Chop[(eqsWithLoop/.Equal[x_,y_]:>x-y)/.solnWithLoop]);
				counter=counter+1;
			];
		(********************************)
	];
	
	(*FindRoot only finds one solution, but keep in mind to following??*)
	(*NOTE: The solution choices become important for max cuts.*)
	(*For instance, a three vertex (1,2,9) needs either <12> = <19> = <29> = 0 or square brackets vanishing*)
	(*If square brackets vanish, then A(1+,2+,9-)=0*)
	(*solnList=Flatten[Table[Join[solnWOLoopList[[i]],solnWithLoopList[[i,j]]],{i,1,Length[solnWOLoopList]},{j,1,Length[solnWithLoopList[[i]]]}],1];*)
	soln=Join[solnWOLoop,solnWithLoop];

	(*Hard-coded for 4-points*)
	analyticCutSoln=Solve[onShellEqs][[1]];
	solvevars=Flatten[Join[{kk[1,3]},Table[kk[i,i],{i,1,3}],Table[kk[loopEdges[[i]],j],{i,1,Length[loopEdges]},{j,1,3}],Table[kk[loopEdges[[i]],loopEdges[[j]]],{i,1,Length[loopEdges]},{j,i,Length[loopEdges]}]]];
	analyticCutSoln=ToRules[Reduce[analyticCutSoln/.Rule->Equal,solvevars]];

	point=(vList/.soln);
	point=Table[{edges[[i]],point[[i]]},{i,1,Length[edges]}];

	{point,momConsRep/.mom->p,analyticCutSoln}

];


(* ::Section::Closed:: *)
(*Sewing Technology*)


(* ::Subsection::Closed:: *)
(*Physical State Projector*)


(*eg. H[\[Mu][5,1],\[Mu][5,3]]:   Used to sew p[i][\[Mu][leg,1]]*\[Epsilon][j][\[Mu][leg,3]] -> -ke[i,j]+... *)
(*projChoice=0 keeps unphysical states*)
H[index1_,index2_,projChoice_Integer:1]:=
Module[{leg,il1,il2,refMom,compRel,proj},

	(*il = "index label"*)
	{leg,il1}=List@@index1;
	il2=(List@@index2)[[2]];
	refMom=ToExpression["ref"<>ToString[leg]];

	(*Completeness Relation in See Dixon Eq. 30*)
		compRel=(-\[Eta][index1,index2]+(p[leg][index1]p[refMom][index2]+p[leg][index2]p[refMom][index1])/AngleBracket[p[leg],p[refMom]]);
		proj=Switch[projChoice,
			(* No physical state projector *)
				0,
					-\[Eta][index1,index2],
			(* physical states *)
				1,
					compRel
			];
	proj
];


(*eg. H[\[Mu][5,1],\[Mu][5,3]]:   Used to sew p[i][\[Mu][leg,1]]*\[Epsilon][j][\[Mu][leg,3]] -> -ke[i,j]+... *)
(*projChoice=0 keeps unphysical states*)
HnoProp[index1_,index2_,projChoice_Integer:1]:=
Module[{leg,il1,il2,refMom,compRel,proj},

	(*il = "index label"*)
	{leg,il1}=List@@index1;
	il2=(List@@index2)[[2]];
	refMom=ToExpression["ref"<>ToString[leg]];

	(*Completeness Relation in See Dixon Eq. 30*)
		compRel=(-2*AngleBracket[p[leg],p[refMom]]*\[Eta][index1,index2]+2p[leg][index1]p[refMom][index2]+2p[leg][index2]p[refMom][index1]);
		proj=Switch[projChoice,
			(* No physical state projector *)
				0,
					{-\[Eta][index1,index2],1},
			(* physical states *)
				1,
					{compRel,R[p[leg]+p[refMom]]}
			];

	proj

];


(* ::Subsection::Closed:: *)
(*Gravity Projectors*)


(*Example*)
(*
2------7->-(1,3)    (2,4)-7->-------3
   |                       |
   |                       |
   ^ 6                     v 8
   |                       |
1------5-<-(1,3)    (2,4)-5-<-------4
*)
(*
             (1,2)|(3,4) (3,4)|(5,6)         
     2 ---<-------->-------->-------->--- 3
              |   7    |   6     |
              |        |         |               
			  |        |         |
              |   8    |   5     |
     1 ----<------->-------->-------->--- 4
             (1,2)|(3,4) (3,4)|(5,6)        
*)

(*type:
	0 = Gravitons, Anti-sym, Dilatons, No physical state projector
	1 = Gravitons, Anti-sym, Dilatons, physical states
	2 = Gravitons, Dilatons, physical states
	3 = Gravitons, physical states
*)
MakeGravProj[legLabel_Integer,{leftIndex1_Integer,rightIndex1_Integer},{leftIndex2_Integer,rightIndex2_Integer},type_Integer]:=
	Module[{proj},
	
	proj=Switch[type,
		(* Gravitons, Anti-sym, Dilatons, No physical state projector *)
			0,
				h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2],0]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2],0],
		(* Gravitons, Anti-sym, Dilatons, physical states *)
			1,
				h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2]],
		(* Gravitons, Dilatons, physical states *)
			2,
				1/2(h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2]]+
					h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,rightIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,leftIndex2]]),
		(* Gravitons, physical states *)
			3,
				1/2(h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2]]+
					h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,rightIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,leftIndex2]])-
				1/(d-2)(h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,rightIndex1]]h0[\[Mu][legLabel,leftIndex2],\[Mu][legLabel,rightIndex2]]),
		(* Gravitons, Anti-sym, Dilatons, physical states - separate left and right copies *)
			4,
				times[h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2]],h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2]]],
		(* Gravitons, Anti-sym, Dilatons, no physical state projector - separate left and right copies *)
			5,
				times[h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,leftIndex2],0],h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,rightIndex2],0]],
		(* The symmeterizing part of the GD projector*)
			6,
				h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,rightIndex2]]h0[\[Mu][legLabel,rightIndex1],\[Mu][legLabel,leftIndex2]],
		(* Dilaton *)
			7,
				-1/(d-2)(h0[\[Mu][legLabel,leftIndex1],\[Mu][legLabel,rightIndex1]]h0[\[Mu][legLabel,leftIndex2],\[Mu][legLabel,rightIndex2]])
		];
		
	proj
];

NMakeGravProj[legLabel_Integer,{leftIndex1_Integer,rightIndex1_Integer},{leftIndex2_Integer,rightIndex2_Integer},type_Integer]:=
	Module[{proj},
		
	proj=MakeGravProj[legLabel,{leftIndex1,rightIndex1},{leftIndex2,rightIndex2},type];
	
	proj /. {h0[index1_, index2_, x___] :> NMatH[index1,index2,x]}
];


(* ::Subsection::Closed:: *)
(*Numeric Matrices*)


(*Makes numeric matrix from output of something like StripPol4[___,1]*)
(*Error message is harmless*)
NMatNoLR[expr_]:=Module[{d, matDotN, indices},
	(* Get the the dimension *)
	d = Length[VEC[1]];

	matDotN=NDot[expr];
	indices=Sort[DeleteDuplicates[Cases[matDotN,\[Mu][__],\[Infinity]]]];
	(*Plugs in all possible values for \[Mu][i,j]'s*)
	(*Error message is harmless*)
	Quiet[(Table@@Prepend[({#,1,d}&/@indices),NFreeIndex[matDotN]])[Sequence@@indices],{Part::pspec}]
];


(*Makes numeric matrix from output of something like StripPol4LR[Mtree4[{1,2,7,8}],1,2]*)
(*Error message is harmless*)
NMat[expr_]:=Module[{d, matDotN, indicesL, indicesR, timesRep},
	(* Get the the dimension *)
	d = Length[VEC[1]];

	matDotN=NDot[expr];
	{indicesL,indicesR}=Sort[DeleteDuplicates[Cases[#,\[Mu][__],\[Infinity]]]]&/@(matDotN /. times->List);
	(*Plugs in all possible values for \[Mu][i,j]'s*)
	timesRep={times[x__,y__] :> Times[
								(Table@@Prepend[({#,1,d}&/@indicesL),NFreeIndex[x]])[Sequence@@indicesL],
								(Table@@Prepend[({#,1,d}&/@indicesR),NFreeIndex[y]])[Sequence@@indicesR]
							]};
	(*Error message is harmless*)
	Quiet[matDotN /.timesRep,{Part::pspec}]
];


(*Numerical physical state projector matrix*)
NMatH[index1_,index2_,x___]:=Module[{d, freeLowerIndexRep,hTemp},
	(* Get the the dimension *)
	d = Length[VEC[1]];

	(*Need to use lower indices for these!!*)
	freeLowerIndexRep=
		{\[Eta][i_,j_]:>g[i,j], 
		p[i_][j_]:>Sum[g[j,k]*VEC[i][[k]],{k,1,d}]};
	
	hTemp = NDot[H[index1,index2,x]];
	Table[hTemp/.freeLowerIndexRep,{index1,1,d},{index2,1,d}][index1,index2]
];


(* ::Subsection::Closed:: *)
(*Matrix Multiplication Method*)


(*Clear[ContractFull]*)
ContractFull[matProduct_]:=Module[{expr,twoContrRep},
	expr=Expand[matProduct];
	twoContrRep={a_[\[Mu][leg_,i_],x_]*b_[y___,\[Mu][leg_,i_],z___]:> ContractMat[a[\[Mu][leg,i],x], b[y,\[Mu][leg,i],z]]};
	expr=expr //. twoContrRep;
	(*If our matrix expression is a sum, make it a list*)
	expr=If[Head[expr]===Plus,List@@expr,{expr}];
	(*Chop[*)Plus@@(LongTensorsContract/@expr)(*]*)
];


ClearAll[ContractMat];
ContractMat[mat1_,mat2_]:=Module[{\[Mu]List1,\[Mu]List2,contIndexList,pos1,pos2,flatList1,flatList2,unconIndices,newMatrixHead},
	(*strip off the indices of the matrices*)
	\[Mu]List1=List@@mat1;
	\[Mu]List2=List@@mat2;
	
	(*The list of identical indices to contract*)
	contIndexList=Intersection[\[Mu]List1,\[Mu]List2];

	(*Find the position of the index in each list*)
	pos1=Flatten[Position[\[Mu]List1,#]&/@contIndexList];
	pos2=Flatten[Position[\[Mu]List2,#]&/@contIndexList];

	(*Make lists to rearrange matrices so that they can be multiplied*)
	flatList1=Join[List/@Delete[Range[Length[\[Mu]List1]],List/@pos1],{pos1}];
	flatList2=Join[{pos2},List/@Delete[Range[Length[\[Mu]List2]],List/@pos2]];

	(*Keep track of uncontracted indices*)
	unconIndices=Join[Delete[\[Mu]List1,List/@pos1],Delete[\[Mu]List2,List/@pos2]];
	
	(*Return the multiplied matrices with the uncontracted indices*)
	newMatrixHead=Flatten[Head[mat1],flatList1].Flatten[Head[mat2],flatList2];
	
	If[NumberQ[newMatrixHead],
		newMatrixHead,
		newMatrixHead@@unconIndices
	]
];


(*This uses the best route (tensors with the most common indices) for contracting a product of tensors*)
(*Aims to control intermediate tensor length*)
(*This step follows all two-index contractions, this will always include H contractions*)
ClearAll[LongTensorsContract];
LongTensorsContract[prodOfTensors_]:=Module[{list,listOfIndices,numCommonIndices,max,maxPos},
	
	list=List@@prodOfTensors;
	
	While[Not[And@@(NumberQ/@list)],
		
		(*Pull all of the indices (if there are any, i.e. not a number) off of each tensor*)
		listOfIndices=Cases[#,\[Mu][__],\[Infinity]]&/@list;

		(*Make a matrix of the number of common indices between terms*)
		numCommonIndices=
			Table[
				If[i<j,
					Length[Intersection[listOfIndices[[i]],listOfIndices[[j]]]],
					0],
			{i,1,Length[listOfIndices]},{j,1,Length[listOfIndices]}];
		
		(*Find the tensor that share the maximum number of indices*)
		max=Max[numCommonIndices];
		maxPos=Position[numCommonIndices,max][[1]];
		
		(*Contract the tensor that have the most indices in common*)
		(*Make a list of that contaction and the other elements left in the list*)
		list=Append[
			Delete[list,List/@maxPos],
			ContractMat[list[[maxPos[[1]]]],list[[maxPos[[2]]]]]];
	];
	(*Return the number*)
	Times@@list
];


(* ::Section::Closed:: *)
(*Numerical Sewing to Get Cut*)


(* Gravity only supports 3-, 4-, and 5-point component trees *)
(* Modified to handle more than 4 external legs, but not tested *)
(* For n-points, external legs must be labeled 1,2,...,n in input *)


(* ::Subsection::Closed:: *)
(*YM*)


(*Clear[SewYM];*)
SewYM[verts_List,projLegChoices_:{}]:=Module[
{d,tree,extEdges,intEdges,cutVerts,cutEdges,n,treesYM,matTreesYM,nProj,ymProjTypes,edgeMatch,ymProjArgs,matHYM,i},
	
	If[Head[VEC[1]] === VEC,
		Print["Numerical vectors, VEC[i], are not set"];
		Return;
	];
	(* Get the the dimension *)
	d = Length[VEC[1]];

	(*All internal edges are cut*)
	{extEdges,intEdges,cutVerts} = GraphInfo[verts];
	cutVerts = tree/@cutVerts;
	cutEdges = intEdges;
	n = Length[extEdges];

	(*YM Component Trees*)
	treesYM = Table[
		cutVerts[[i]]
			/.{tree[a_] :> If[Length[a] === 3, 
								Tree3[a,i], 
								StripPol[n,Atree[a],i]]}
			,{i,1,Length[cutVerts]}];
	(*Make them numeric matrices*)
	matTreesYM = Chop[Expand[NMatNoLR/@treesYM]];

	(*Cut Projectors*)
	nProj=1;(*For physical state projectors*)
	(*projLegChoices tell what cut legs to use the projectors on. Others are assigned -eta_munu*)
	If[projLegChoices==={}, 
		ymProjTypes = ConstantArray[nProj,Length[cutEdges]];,
		ymProjTypes = ReplacePart[ConstantArray[0,Length[cutEdges]],Rule[#,nProj]&/@Flatten[Position[cutEdges,#]&/@projLegChoices]];
	];

	edgeMatch = {#,Position[cutVerts,#][[All,1]][[1]],Position[cutVerts,-#][[All,1]][[1]]}&/@cutEdges;
	ymProjArgs = {\[Mu][#[[1]],#[[2]]],\[Mu][#[[1]],#[[3]]]}&/@edgeMatch;
	ymProjArgs = Table[AppendTo[ymProjArgs[[i]],ymProjTypes[[i]]],{i,1,Length[ymProjArgs]}];
	(*Generate ref momenta for cut leg projectors*)
	For[i=1,i<=Length[cutEdges],i++,VEC[ToExpression["ref"<>ToString[cutEdges[[i]]]]] = GenVectorNC[d,0]];
	(*Physical state projectors for the cut legs in YM*)
	matHYM = NMatH@@#&/@ymProjArgs;

	(*YM*)
	ContractFull[(Times@@matHYM)*(Times@@matTreesYM)]
];


(* ::Subsection::Closed:: *)
(*Gravity*)


SewGDA[verts_List,projLegChoices_:{}]:= SewGrav[verts,1,projLegChoices];
SewGD[verts_List,projLegChoices_:{}]:= SewGrav[verts,2,projLegChoices];
SewG[verts_List,projLegChoices_:{}]:= SewGrav[verts,3,projLegChoices];


Clear[SewGrav];
SewGrav[verts_List,nProj_Integer,projLegChoices_:{}]:=Module[
{d,tree,extEdges,intEdges,cutVerts,cutEdges,n,treesGrav,matTreesGrav,gravProjTypes,edgeMatch,gravProjArgs,matHGrav,i},
	If[Head[VEC[1]] === VEC,
		Print["Numerical vectors, VEC[i], are not set"];
		Return;
	];
	(* Get the the dimension *)
	d = Length[VEC[1]];

	(*All internal edges are cut*)
	{extEdges,intEdges,cutVerts} = GraphInfo[verts];
	cutVerts = tree/@cutVerts;
	cutEdges = intEdges;
	n = Length[extEdges];

	(*Gravity Component Trees*)
	treesGrav=Table[cutVerts[[i]]
					/.{tree[a_]:>Switch[Length[a],
									3,Mtree3[a,2*i-1,2*i],
									4,StripPolLR[n,Mtree4[a],2*i-1,2*i],
									5,StripPolLR[n,Mtree5[a],2*i-1,2*i]]}
			,{i,1,Length[cutVerts]}];
	(*Make them numeric matrices*)
	matTreesGrav=Chop[NMat/@treesGrav];

	(*Cut Projectors*)
	(*projLegChoices tell what cut legs to use the projectors on. Others are assigned -eta_munu*)
	If[projLegChoices==={}, 
		gravProjTypes=ConstantArray[nProj,Length[cutEdges]];,
		gravProjTypes=ReplacePart[ConstantArray[0,Length[cutEdges]],Rule[#,nProj]&/@Flatten[Position[cutEdges,#]&/@projLegChoices]];
	];

	edgeMatch={#,Position[cutVerts,#][[All,1]][[1]],Position[cutVerts,-#][[All,1]][[1]]}&/@cutEdges;
	gravProjArgs={#[[1]],{2*#[[2]]-1,2*#[[2]]},{2*#[[3]]-1,2*#[[3]]}}&/@edgeMatch;
	gravProjArgs=Table[AppendTo[gravProjArgs[[i]],gravProjTypes[[i]]],{i,1,Length[gravProjArgs]}];
	(*Generate ref momenta for cut leg projectors*)
	For[i=1,i<=Length[cutEdges],i++,VEC[ToExpression["ref"<>ToString[cutEdges[[i]]]]]=GenVectorNC[d,0]];
	(*The physical state projector for gravitons, dilaton, and antisym is H[\[Mu],\[Rho]]H[\[Nu],\[Sigma]]*)
	matHGrav=NMakeGravProj@@#&/@gravProjArgs;

	ContractFull[(Times@@matHGrav)*(Times@@matTreesGrav)]
];


(* ::Section::Closed:: *)
(*Find Graphs That Contribute to Cut*)


(* ::Subsection::Closed:: *)
(*Cut Tree Amplitude Structure *)


(*Edge data structure*)
(*Canonicalizes the sign of momentum combination inside*)
(*This could mess up idea of incoming/outgoing, especially for external legs*)
(*especially if automorphisms switch + and - ed's*)
ed /: ed[a_] := -ed[-a] /; Or[Count[{a[[1]]},Times[-1,_]] > 0, Count[{a},Times[-1,_]]>0];


(* Takes the contact graph vertex structures and builds the 
contibuting diagrams using color info *)
Clear[TreeAmpStruc]
TreeAmpStruc[legLabels_List,ColorOrderedQ_]:=Module[
{n,uniqueGraphs,ib,momCons,verts,basis,rep,relab},
	
	n = Length[legLabels];
	
	If[n === 3,
		uniqueGraphs = {Vs[{-ed[p[1]]},{-ed[p[2]]},{ed[p[1]+p[2]]},{ed[p[1]],ed[p[2]],ed[-p[1]-p[2]]}]};,

		ib = GetGraphBasis[n,0];
		momCons = (MomRules/@ib)/.l[x_]:>p[x];
		verts = Table[Vs@@((ib[[i]]/.{x_Integer:>Sign[x]*ed[p[Abs[x]]]})/.(momCons[[i]])),{i,1,Length[ib]}];
		If[ColorOrderedQ === False, 
			basis = GetSymBasis[n];, 
			basis = GetCOBasis[n];
		];
		rep = {fColor[j_,list_]:>((verts[[j]]/.Table[p[k]->p[list[[k]]],{k,1,Length[list]}])/.(momCons[[j]]))};
		uniqueGraphs = basis /. rep;
	];
	
	relab[list_] := Join[
					Table[p[i]->Sign[list[[i]]]p[Abs[list[[i]]]],{i,1,Length[list]-1}],
					{ed[Sum[p[i],{i,1,Length[list]-1}]]:>ed[-Sign[list[[-1]]]p[Abs[list[[-1]]]]]}
				];

	uniqueGraphs/.relab[legLabels]
];


(* ::Subsection::Closed:: *)
(*Join Trees*)


Clear[JoinVsTrees]
(*Joins many Vs trees together*)
JoinVsTrees[vsList_List]:=Module[{extEds,extEdsFinal,intEds,extVars,intVars,intReps,intEdsFinal,preFactor},
	
	extEds=Cases[#,{_},\[Infinity]]&/@vsList;
	extEdsFinal=(Join@@extEds)//.{{a___,{b_},c___,{-b_},d___}:>{a,c,d},{a___,{-b_},c___,{b_},d___}:>{a,c,d}};

	intEds=Cases[#,{_,_,_},\[Infinity]]&/@vsList;
	extVars=Variables[#]&/@extEds;
	intVars=Table[Complement[Variables[intEds[[i]]],extVars[[i]]],{i,1,Length[vsList]}];

	(*Label each internal edge by the tree where it came from*)
	(*Need to worry about, eg., ed[p[1]+p[2]] showing up multiple times
	in say
	\                          /
	 \___bubble___bubble ___ _/
	 /                        \
	/                          \
	*)
	intReps=Table[intVars[[i]]/.ed[a_,n___]:>(ed[a,n]->ed[a,n+i]),{i,1,Length[intVars]}];
	intEdsFinal=Join@@Table[intEds[[i]]/.intReps[[i]],{i,1,Length[intEds]}];

	(*If the vsList comes with coefficients, make sure to multiply them*)
	preFactor=Times@@(vsList/.Vs[__]:>1);

	preFactor*(Vs@@Join[extEdsFinal,intEdsFinal])
];


(* ::Subsection::Closed:: *)
(*Isomorphisms*)
(*(Input: Vs vertex list structure, vertex basis structures)*)


(*Find the basis element that is isomorphic to our input vertex list*)
FindIBIso[vertList_Vs, vertBasis_List]:=Module[
{ib,vList,extEds,intEds,vars,repN,vertListN,
isoIndex,isoSig,iso,
invRepN,isopRep,iso\[Epsilon]Rep,isoRep,isoTemp,Temp},

	(* If vertBasis is in BG form, change from BG to vert structure *)
	If[Head[vertBasis[[1]]] === BG,
		ib = vertBasis /. BG[x__] :> x[[1]];,
		ib = vertBasis;
	];

	vList = vertList/.Vs->List;
	(*replace the ed[p[i]]'s with integers so we can build graphs*)
	(*Need external edges to be 1-n*)
	extEds = Variables[Cases[vList,{_},\[Infinity]]];
	intEds = Complement[Variables[Cases[vList,{_,_,_},\[Infinity]]],extEds];
	(*Make sure external edges are first in your variable list so they get the lowest numbers*)
	vars = Join[extEds,intEds];
	repN = Table[vars[[i]]->i,{i,1,Length[vars]}];
	
	(*use reps to convert ed[]'s to numbers*)
	vertListN = vList/.repN;

	{isoIndex, isoSig, iso} = FindIso[vertListN, ib];
	
	(*Change back to original labeling*)
	invRepN=(repN/.Rule->Temp)/.Temp[x_,y_] :> Rule[p[y],x];
	isopRep = invRepN/.ed[a_,m___]:>a;
	iso\[Epsilon]Rep = (invRepN/.Rule[a_,b_]:>Rule[a/.p->\[Epsilon],b/.ed[x_,m___]:>\[Epsilon][x]])/.{\[Epsilon][p[i_]]:>\[Epsilon][i], \[Epsilon][x_Plus]:>\[Epsilon][Length[List@@(x)]+1]};
	isoRep = Join[isopRep,iso\[Epsilon]Rep];
	isoTemp = iso/.{Rule[a_,b_]:>Temp[a,b/.l->p]};
	iso = isoTemp/.Temp[a_,b_]:>Rule[a,b/.isoRep];
	
	{isoIndex,isoSig,iso}
];


(* ::Section::Closed:: *)
(*Analytical Cutting*)


(* ::Subsection::Closed:: *)
(*Get Cut Isomorphisms*)


(* Lines up with GetCutIsos as long as we consistently find intVerts (and Tuples works the same) *)
Clear[GetJacOnCut];
GetJacOnCut[verts_List] := Module[{intVerts,colorFactors,colorEqs,colorSoln,jacOnTheCuts},

	intVerts = GraphInfo[verts][[3]];

	colorFactors=
		Expand[
			Times@@#&/@(
				Tuples[
					If[Length[#]===3,
						{fColor[3]},
						(GetSymBasis[Length[#]]/.GetDDMReps[Length[#]]) /. fColor[i_,list_]:>fColor[{i},{#[[Ordering[list]]]}]
					]&/@intVerts]
			)/.{fColor[3]:>1}
		];
	colorFactors = colorFactors//.{fColor[i1_List,list1_List]*fColor[i2_List,list2_List]:>fColor[Join[i1,i2],Join[list1,list2]]};
	colorEqs = Table[COLOR[i] == colorFactors[[i]],{i,1,Length[colorFactors]}];
	colorSoln = ToRules[Reduce[colorEqs, Variables[colorFactors]]];
	
	jacOnTheCuts = Table[{COLOR[i]/.colorSoln},{i,1,Length[colorFactors]}]
];


(* Input: Contact diagram structure, Color-ordered?, Loop edges choice, graph basis (vertex or BG structure) *)
(*Clear[GetCutIsos];*)
GetCutIsos[verts_List, loopEdges_List, ib_List, colorOrderedQ_]:=Module[
{intVerts,edges,extEdges,momConsRep,ampGraphs,ampGraphsSigns,len,isos,jac,i},

	intVerts = GraphInfo[verts][[3]];

	(*Momentum Conservation*)
	momConsRep = MomRules[verts, loopEdges];
	
	ampGraphs = JoinVsTrees/@Tuples[ TreeAmpStruc[#,colorOrderedQ]&/@intVerts ];
	(*signs or prefactors that could come from color-ordered basis *)
	ampGraphsSigns = ampGraphs/.Vs[__]:>1;
	ampGraphs = ampGraphs*ampGraphsSigns;
	
	len = Length[ampGraphs];
	isos = Range[len];
	For[i=1,i<=Length[ampGraphs],i++,
		(*Print["Finding isomorphism "<>ToString[i]<>" of "<>ToString[len]];*)
		isos[[i]] = FindIBIso[ampGraphs[[i]],ib]/.{l->p};
	];
	isos = isos/. Rule[x_,y_]:>Rule[x,y/.momConsRep];
	(* Include ampGraphsSigns in isos signature if color-ordered; otherwise, make all signatures 1*)
	If[colorOrderedQ === True,
		isos[[All,2]] = ampGraphsSigns*(isos[[All,2]]);,

		(*isos[[All,2]]=Abs[isos[[All,2]]];*)
		isos[[All,2]] = ampGraphsSigns*(isos[[All,2]]);
		(*Append color information*)
		jac = GetJacOnCut[verts];
		isos = ArrayFlatten[{{isos,jac}}];
	];

	isos
];


(* ::Subsection::Closed:: *)
(*Analytical Cutting*)


(* Color-ordered *)
CutYM[verts_List, loopEdges_List, exprList0_List, ib_List, listQ_:False] := CutAnalytic[verts,loopEdges,exprList0,ib,True,listQ,False];
(* Color-dressed *)
CutYMColor[verts_List, loopEdges_List, exprList0_List, ib_List, listQ_:False] := CutAnalytic[verts,loopEdges,exprList0,ib,False,listQ,False];
(* Gravity *)
CutGravity[verts_List, loopEdges_List, exprList0_List, ib_List, listQ_:False] := CutAnalytic[verts,loopEdges,exprList0,ib,False,listQ,True];


(*
Note: Final factors of (-I)^# of cut legs is good for YM and Gravity.
Reason: 
	YMprop = (I/R)*[-\[Eta]]  ~> [-H] (physical state projector) for sewing, 
	GravProp = (I/R)*[(1/2)(\[Eta]\[Eta] + \[Eta]\[Eta]) - (1/(Ds-2))\[Eta]\[Eta]] ~> [(1/2)(HH + HH) - (1/(Ds-2))HH] for sewing
*)


(* verts is the graph where all external edges are cut *)
(* exprList0 are the Feynman or ansatze expressions *)
(*The list should be in a form with numerators and propagators separated*)
(* = {{num_Plus,props_List},{num,props},{num,props},...}*)
Clear[CutAnalytic];
CutAnalytic[verts_List, loopEdges_List, exprList0_List, ib_List, colorOrderedQ_:False, listQ_:False, gravityQ_:False]:=Module[
{isos,exprList,vars,cutReps,dispatchIsos,appliedIsosList,
remainingProps,propVars,propPowers,maxPowers,preFactorPropPowers,preFactors,cutEdges},

	(******Isomorphism data*********)
	(*****************************************)
	isos = GetCutIsos[verts,loopEdges,ib,colorOrderedQ];
	If[gravityQ,
		isos = isos/.(Rule[\[Epsilon][i_],\[Epsilon][j_]]:>Sequence[Rule[\[Epsilon][i],\[Epsilon][j]],Rule[\[Epsilon]2[i],\[Epsilon]2[j]]]);
	];
	(*Delete any unfound isos - like triangles in N=4*)
	isos = DeleteCases[isos,{0,__}];
	(*****************************************)
	(*****************************************)

	(******Ansatz or Feyman rule data*********)
	(*****************************************)
	(*The list should be in a form with numerators and propagators separated*)
	(* = {{num_Plus,props_List},{num,props},...}*)
	exprList = exprList0/.{lk->kk,ll->kk,le->ke,lE->kE};
	(*Add all graphs that contribute to cut*)
	vars = Table[GetVars[exprList[[isos[[i,1]]]],SYMBOLSLIST],{i,1,Length[isos]}];
	(*****************************************)
	(*****************************************)

	(* Cut information *)
	cutReps = GetCutReps[verts,loopEdges];
	
	(*To speed up isomorphism application make dispatch list*)
	dispatchIsos = Table[Dispatch[(#->(Relab[isos[[i,3]],#]/.cutReps))&/@vars[[i]]],{i,1,Length[isos]}];
	(*Apply isomorphisms*)
	appliedIsosList = Table[exprList[[isos[[i,1]]]]/.dispatchIsos[[i]],{i,1,Length[isos]}]/.R[x_]:>Expand[AngleBracket[x,x]/.cutReps](*Check!!*);
	(*Make sure to use signature in color-ordered YM! (Not sure about color-dressed)*)
	(*Be careful to only apply to numerator not props*)
	(*Don't apply for gravity*)
	appliedIsosList = Table[
						{If[(*colorOrderedQ &&*)!gravityQ,isos[[i,2]],1] * appliedIsosList[[i,1]], appliedIsosList[[i,2]]}
					,{i,1,Length[appliedIsosList]}];

(*If we just want the relabelings return appliedIsosList*)
If[listQ===True,
	
	appliedIsosList,	

	(*No prefactors or color decomposition necessary*)
	(*Essentially multiplies by the common denominator*) 
	remainingProps = Table[DeleteCases[appliedIsosList[[i,2]],1],{i,1,Length[appliedIsosList]}];
	propVars = DeleteDuplicates[Flatten[remainingProps]];
	propPowers = Table[Count[#,propVars[[i]]],{i,1,Length[propVars]}]&/@remainingProps;
	maxPowers = Max[#]&/@Transpose[propPowers];
	preFactorPropPowers = (maxPowers-#)&/@propPowers;
	preFactors = ((Times@@(propVars^#))&/@preFactorPropPowers)//Simplify;

	cutEdges = GraphInfo[verts][[2]];
	(*{num, denom}*)
	{
		(-I)^Length[cutEdges]*
			Sum[
				(*color*) If[(!colorOrderedQ && !gravityQ),isos[[i,4]],1]*
				preFactors[[i]]*appliedIsosList[[i,1]],
			{i,1,Length[preFactors]}],
		Times@@(propVars^maxPowers)
	}
]

];


(* ::Subtitle::Closed:: *)
(****************************************)
(*End*)
(****************************************)


End[]; (* End `Private` *)
EndPackage[];
