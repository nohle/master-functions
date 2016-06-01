(* ::Package:: *)

(* Don't forget to set MASTERFUNCTIONSDIRECTORY in the *)
(* Public Initializations of MasterFunctions.m ! *)
masterFunctionsDirectory = $HomeDirectory<>"/Desktop/MasterFunctions/";
Quiet[Get[masterFunctionsDirectory<>"MasterFunctions.m"];];


(* ::Subtitle::Closed:: *)
(**********************************************************)
(*Basic Examples*)
(**********************************************************)


(* ::Subsection::Closed:: *)
(*4-Point Tree Using Basic Functions*)


(* Let's look at a 4-point color-ordered tree amplitude *)
n = 4;
L = 0;


(* We get the graph topologies from... *)
graphs = GetGraphBasis[n, L]
(* We can also read off the fabc's from the 3-vertices *)


(* To visualize the graphs use... *)
DisplayGraph/@graphs


(* To look at the relabeling symmetries of the graph we use *)
(* (We transpose the result so it doesn't look so weird) *)
LabelwiseAutomorphisms[ graphs[[1]] ] //Transpose
(* 
These are:
	- identity
	- antisym. in 3<->4
	- antisym. in 1<->2
	- sym. in 1<->2 and 3<->4
	- sym. in reflection
	- 1<->2 + reflection
	- 3<->4 + reflection
	- 1<->2 + 3<->4 + reflection
*)


(* Let's look at color *)
(* Here, fColor[i, list] gives the fabc contractions of graphs[[i]] with external ordering "list" *)

(* What color factors are independent under (anti)symmetries of the graph? *)
GetSymBasis[n]
(* What are the color factors of the DDM basis (arXiv:hep-ph/9910563) ? *)
GetDDMBasis[n]
(* For fColor[1,{1,2,3,4}], what graphs contribute to the color-ordered amplitude?*)
GetCOBasis[n]


(* If we want to see how other color factors fit into the bases, we can use, for instance... *)
symReps = GetSymReps[n]
(* So, *)
fColor[1,{2,1,3,4}] /. symReps
(* Be careful though! For higher points, there can be A LOT of color factors! *)


(* What Jacobi identities do the graph(s) of GetGraphBasis[n] satisfy? *)
test = GetGraphJacobis[n, L];
test /. RG[x_, y_, z_]:>RG[x,y[[1]],z[[1]]]
(* At tree level, we can get this information from color information as well, specifically, GetDDMReps[n] *)


(* These two graphs of GetCOBasis[n]---the s and t channels---are computed using... *)
amp = Atree[{1,2,3,4}]
(* Note: Atree[n] == Atree[Range[n]] *)


(* How can we strip off the polarization vectors? *)
(* We will use the structures given in "Structure Definitions" *)
ampStripped = amp /. {ee[i_,j_] :> \[Eta][\[Mu][i],\[Mu][j]], ke[i_,j_] :> p[i][\[Mu][j]]}


(* We can recontract the polarization vectors using "Propagate" *)
amp2 = Factor[
	Propagate[\[Epsilon][1][\[Mu][1]]\[Epsilon][2][\[Mu][2]]\[Epsilon][3][\[Mu][3]]\[Epsilon][4][\[Mu][4]] * ampStripped]
]

(* Are they the same? *)
Expand[amp2 - amp] === 0
(* Yes *)


(* Now, how can we get the other color-ordered amplitude, 1324, in the DDM basis? *)
GetDDMBasis[n][[2]]

(* We can just use *)
amp1324 = Atree[{1,3,2,4}]

(* Or we can relabel "amp", the 1234 ordering *)
reps = {
			p[2] -> p[3], \[Epsilon][2] -> \[Epsilon][3],
			p[3] -> p[2], \[Epsilon][3] -> \[Epsilon][2]
};
amp1324Other = Relab[reps, amp]
(* Also, on-shell conditions*)
amp1324Other = amp1324 /. kk[i_,i_]:>0

(* Are they the same?? *)
Expand[amp1324 - amp1324Other] === 0
(* Yes *)


(* Notice that in "amp" and "amp1324" we see kk[1,2], kk[2,3], kk[1,3], etc. *)
(* We can see the variables in a function using GetVars[] (or Variables[]) *)
GetVars[amp]
GetVars[amp1324]
(* The nice thing about GetVars is that we can exclude some variables *)
GetVars[amp, {ee[1,2]}]
(* or if we just want to look at kk's or s's *)
GetVars[amp, {ee[__], ke[__]}]
GetVars[amp1324, {ee[__], ke[__]}]


(* We see s, t, and u, but s+t+u=0 *)
(* So, we should really apply all mom cons / on-shell / transversality conditions *)
amp1324 = MomConsN[n, amp1324]

(* Ugly, so Factor it *)
amp1324 = Factor[amp1324]


(* Sometimes, using s's can make things prettier *)
(* To do that, use a "True" argument in MomConsN *)
Factor[MomConsN[n, amp1324, True]]


(* ::Subsection::Closed:: *)
(*6-Point Tree*)


(* Let's look at 6-point trees now *)
n = 6;
L = 0;


graphs = GetGraphBasis[n, L]
DisplayGraph/@graphs


(* The color-oredered amplitude (in terms of s's) is...*)
amp = Atree[6];
(* Atree[n] should automatically put things in terms of s's not kk's, if not,
we can use the "True" argument in MomConsN *)
(*amp = ExpandAll[MomConsN[n, Atree[6], True]];*)


(* This is built out of the Feynman rules for the graphs *)
GetFeynYM[n,L]
(* Note that the numerator and propagators are separated *)
(* The quantities R[x] represent x^2 *)
R[p[1]+p[2]+p[3]]
MomConsN[n, R[p[1]+p[2]+p[3]], True]
(* or *)
R[p[1]+p[2]+p[3]] /. R[x__]:>AngleBracket[x,x]
(* or *)
MomConsN[n, R[p[1]+p[2]+p[3]] /. R[x__]:>AngleBracket[x,x] ]


(* We can build Feynman rules one graph at a time (though these are not saved) *)
feyn2 = FeynmanYM[graphs[[2]]];
(* We need to do a little extra work though *)
(* It doesn't replace internal legs using momentum conservation (useful for speed with numerics) *)
(* We can use MomRules[] to find out the momentum conservation rules for the graph *)
momRules = MomRules[graphs[[2]]]
feyn2MomCons = Expand[MomConsN[n, Relab[momRules, feyn2], True]]


(* Say we are given a 6-point graph... *)
unknownGraph = {{-10},{-4},{-2},{-1},{-5},{-12},{10,4,16},{2,1,17},{5,12,18},{-16,-17,-18}};
DisplayGraph[unknownGraph]
(* Let's compute it's Feynman rules *)
feynUnknown = FeynmanYM[unknownGraph];


(* This graph looks familiar though *)
(* Can we find a graph in the basis that this is isomorphic to? *)
FindIso[unknownGraph, graphs]
(* Yes! The output, {i, sign, replacements} tell us that "unknownGraph" is isomorphic to
	graphs[[i]], with the replacement rules "replacements" and sign "sign". *)

(* If we want to map the other way, i.e. from unknownGraph to graphs[[2]], then... *)
iso = FindIso[unknownGraph, graphs, True]
(* In both cases, there are many possible isomorphisms (due to automorphisms), but this function
	just chooses one. *)


(* Now, we can check to see that relabeling feynUnknown with iso should give us feyn2 *)
(* Note: Use {iso[[2]],1} so we don't change sign of propagators too *)
feyn2Check = {iso[[2]],1} * Relab[iso[[3]], feynUnknown];
(* Are they the same ??? *)
feyn2Check === feyn2
(* No! *)
(* Doesn't look good! *)
(* What about mom cons? *)
feyn2mom = Relab[momRules, feyn2];
feyn2Check = Relab[momRules, feyn2Check];
(* Are they the same? *)
feyn2Check === feyn2mom
(* No! *)
(* BUT! Remember there are a lot of different isomorphisms from one graph to the other,
	so we shouldn't expect immediate agreement. If we apply MomConsN[] and expand,
	they should match.*)
(* The propagators match *)
feyn2Check[[2]] === feyn2mom[[2]]
(* Yes! *)
(* The numerators... *)
Expand[MomConsN[n, feyn2Check[[1]] - feyn2mom[[1]]]] === 0
(* ... match!*)


(* Let's pretend that it's expensive to apply MomConsN and Expand[] *)
(* Generating numerics can help a lot!*)

(* The "contact diagram" input. Was initially intended for loops *)
verts = {{-1},{-2},{-3},{-4},{-5},{-6},{1,2,3,4,5,6}};
(* Cut loop momenta are in 4-2\[Epsilon] dimensions. Say \[Epsilon] = -3. Then, dim = 10 *)
(* For trees, it doesn't matter *)
dim = 10;
dimPhys = 4;
(* Generates numerical vectors for our momenta *)
(* The {} are the loop momentum labels... there are none here *)
Print["Numerical Vectors:"];
numPoint = CutNumerics[verts,{},dim,dimPhys]
(* Can do some momentum conservation and on-shell checks here if you're skeptical*)
(* .....................*)

(* Sets the global variables *)
SetVec[numPoint]
(* Sets the global variables POL[1],...,POL[6] for an MHV amplitude *)
GenPol[{"m","m","p","p","p","p"}]

(* Converts symbolic to numeric *)
(* Numerators: *)
Print["Check numerators numerically (should be the same):"];
num = NDot[ feyn2Check[[1]] ]
num2 = NDot[ feyn2mom[[1]] ]
Round[num/num2]===1
(* Propagator lists: *)
Print["Check propagators numerically (should be the same):"];
prop = NDot[ Times@@(feyn2Check[[2]]) ]
prop2 = NDot[ Times@@(feyn2mom[[2]]) ]
Round[prop/prop2]===1
(* They agree! *)

(* Remember the easier way of getting the Feynman results for all diagrams *)
feynFull = GetFeynYM[n,L];
{num2, props2} = {feynFull[[2,1]], Times@@(feynFull[[2,2]])};
Print["Numerator:"];
NDot[ num2 ]
Print["Propagator:"];
NDot[ props2 ]
(* Checks out! *)


(* ::Subsection::Closed:: *)
(*6-Point Tree Numerics (Compare to Parke-Taylor Formula)*)


n = 6;
L = 0;
amp = Atree[Range[n]];


(* Generate numerics... *)
verts = {{-1},{-2},{-3},{-4},{-5},{-6},{1,2,3,4,5,6}};
dim = 5;
dimPhys = 4;
numPoint = CutNumerics[verts,{},dim,dimPhys];
SetVec[numPoint];
GenPol[{"m","m","p","p","p","p"}];

(* Compare "amp" to the known MHV expression *)
ampN = NDot[amp]
(* Compare to MHV amplitude in spinor helicity *)
(* "ab" stands for "angle bracket". "sb" is "square bracket" *)
(* See "4-D Spinor Technology" under "Functions" -> "General Numerics" *)
ampMHV = (I * ab[VEC[1],VEC[2]]^4)/(ab[VEC[1],VEC[2]]ab[VEC[2],VEC[3]]ab[VEC[3],VEC[4]]ab[VEC[4],VEC[5]]ab[VEC[5],VEC[6]]ab[VEC[6],VEC[1]])
(* Or for the lazy *)
ampMHV = (I * ab[1,2]^4)/(ab[1,2]ab[2,3]ab[3,4]ab[4,5]ab[5,6]ab[6,1])
Round[ ampN/ampMHV ] === 1


(* ::Subtitle::Closed:: *)
(**********************************************************)
(*Two-Loop BCJ-on-the-Cuts: Unitarity Cuts Check*)
(**********************************************************)


(* We check the results of http://arxiv.org/abs/1510.03448 *)
(* Download the "numerators.m" ancillary file and rename it "numeratorsBCJ2Loop.m" *)


(*
By "check," we mean verify the spanning set of unitarity cuts through 
a variety of methods. Thus, we are also checking the correctness of our 
technology for cutting, sewing, numerics, isomorphisms, Feynman rules, etc. 
*)


(* ::Subsection::Closed:: *)
(*General Setup (Run me)*)


(*********Spanning set of cuts**********)
(* Non-planar are only relevant for color-ordered YM *)
(* Non-planar cuts should follow automatically *)
cutVVverts={{-1},{-2},{-3},{-4},{1,2,-7,-5},{3,4,-6,-8},{5,7,8,6}};
cutHVverts={{-1},{-2},{-3},{-4},{1,9,10,6},{2,8,-10,-9},{3,4,-6,-8}};
cutHVNPverts={{-1},{-2},{-3},{-4},{1,9,10,6},{2,8,-9,-10},{3,4,-6,-8}};
cut3verts={{-1},{-2},{-3},{-4},{4,1,9,10,11},{2,3,-11,-10,-9}};
cut3NPverts={{-1},{-2},{-3},{-4},{4,1,9,10,11},{2,3,-11,-9,-10}};


(* Expression from http://arxiv.org/abs/1510.03448 *)
SetDirectory[NotebookDirectory[]];
Get["numeratorsBCJ2Loop.m"];//Timing
(*put info in old form*)
bcj =
	Table[
		{
			numerators[[i]],
			Sort[Flatten[List@@(propagators[[i]])/.{x_^2:>R[x], x_^4:>{R[x],R[x]}, x_^6:>{R[x],R[x],R[x]}}]/.k->p]
		}
	,{i,1,Length[numerators]}];
ibBCJ = vertices;


secondCopyReps = {ee->EE,ke->kE};


(* ::Subsection::Closed:: *)
(*Check Yang-Mills*)


(* Get two-loop Feynman expression *)
(* Note that the graph basis here, ibFeyn, is different than ib (or vertices) *)
(* The cutting technology doesn't care. A cut is a unique object, so it shouldn't matter *)
feyn = Expand[GetFeynYM[4,2]];//Timing
ibFeyn = GetGraphBasis[4,2];//Timing


(******************************************)
(*Choose cut*)
(******************************************)
verts = cutVVverts;
Print["We check the cut:"];
verts
loopEdges = ChooseLoopMom[verts][[1]];
Print["Contributions to cut:"];
Print["{Diagram #, Signature, Relabeling from Diagram Basis 'ibFeyn'}"];
GetCutIsos[verts,loopEdges,ibFeyn,True]
(******************************************)

(********Generate cut numerics*******)
Print["Generating cut numerics..."];
d=6;
dPhys=4;
point = CutNumerics[verts,loopEdges,d,dPhys];//Timing
SetVec[point];
helList={"m","m","p","p"};
GenPol[helList];

(********Get cut by sewing trees******)
Print["Getting numerical cut by sewing trees..."];
cutSew = SewYM[verts,loopEdges];//Timing

(********Get cut by cutting Feynman rules*********)
Print["Getting numerical cut by cutting Feynman rules..."];
{cutFeynNum,cutFeynDenom} = CutYM[verts,loopEdges,feyn,ibFeyn];//Timing
cutFeynN = cutFeynNum/cutFeynDenom;
vars = GetVars[cutFeynN,{Ds}];
NReps = Dispatch[#->NDot[#]&/@vars];
cutFeynN = cutFeynN/.NReps;//Timing

(********Get cut by cutting ansatz *********)
Print["Getting numerical cut by cutting BCJ-on-the-cuts..."];
{cutBCJNum,cutBCJDenom} = CutYM[verts,loopEdges,bcj,ibBCJ];//Timing
cutBCJN = cutBCJNum/cutBCJDenom;
vars = GetVars[cutBCJN,{a[__]}];
NReps = Dispatch[#->NDot[#]&/@vars];
cutBCJN = cutBCJN/.NReps;//Timing

(******** Three different methods **********)
cs = cutSew;
cf = Chop[Expand[cutFeynN /. {Ds -> d}]];
(*Put I^L back in*)
cb = Chop[Expand[I^2 * cutBCJN /. Dispatch[{Ds->d}]]];
Print["Numerical sewing matches Feynman matches BCJ-on-the-cuts?"];
(Round[Expand[cs/cf]] === 1) && (Round[Expand[cs/cb]] === 1)

(*********** Analytical Check ***********)
Print["Checking analytical cut of BCJ-on-the-cuts against Feynman rules..."];
ratio = (cutFeynDenom/cutBCJDenom);
bool = (Expand[cutFeynNum - ratio*(I^2*cutBCJNum)]===0);//Timing
Print["Same analytical cut as Feynman rules?"];
bool


(* ::Subsection::Closed:: *)
(*Check Gravity ( Double Copy: (BCJ-on-the-Cuts)^2 )*)


(******************************************)
(*Choose cut*)
(******************************************)
verts = cutHVverts;
loopEdges = ChooseLoopMom[verts][[1]];
GetCutIsos[verts,loopEdges,ibBCJ,True];
(******************************************)

(********Generate cut numerics*******)
Print["Generating cut numerics..."];
d = 6;
dPhys = 4;
point = CutNumerics[verts,loopEdges,d,dPhys];//Timing
SetVec[point];
helList = {"m","m","p","p"};
GenPol[helList];
helList2 = {"m","m","p","p"};
GenPol2[helList2];

(********Get cut by sewing trees******)
Print["Getting numerical cut by sewing trees..."];
cutSew = SewGDA[verts,loopEdges];//Timing

(********Get cut by cutting (BCJ-on-the-cuts)^2 *********)
Print["Getting numerical cut by cutting (BCJ-on-the-cuts)^2..."];
grav = Table[{(bcj[[i,1]]*(bcj[[i,1]]/.secondCopyReps))/.Ds->d, bcj[[i,2]]}, {i,1,Length[bcj]}];
{cutBCJNum,cutBCJDenom} = CutGravity[verts,loopEdges,grav,ibBCJ];//Timing
cutBCJN = cutBCJNum/cutBCJDenom;
Clear[cutBCJNum,cutBCJDenom];
vars = GetVars[cutBCJN,{Ds}];
NReps = Dispatch[#->NDot[#]&/@vars];
cutBCJN = cutBCJN/.NReps;//Timing

cs = cutSew;
(*
The factor of I^3 comes from:
	1) We factored out an I^L from the YM numerators (See Eq. 2.1 of http://arxiv.org/pdf/1303.6605.pdf)
	2) M ~ I^(L+1) and YM ~ I^L (See Eq. 2.5 of http://arxiv.org/pdf/1303.6605.pdf)
*)
cb = I^3*Chop[Expand[cutBCJN]];

Print["Numerical sewing matches cutting (BCJ-on-the-cuts)^2 ?"];
Round[cb/cs] === 1


(* ::Subsection::Closed:: *)
(*Intermission: Generate new Feynman expressions (with fermions too for N = 4) in the same graph basis as BCJ-on-the-Cuts*)


feyn = FeynmanYM[#,1]&/@vertices;//Timing
loopMom = ChooseLoopMom[#][[1]]&/@vertices;
feyn = Table[{
				Expand[MomConsN[4, Relab[MomRules[vertices[[i]], loopMom[[i]]], feyn[[i,1]]]]],
				ExpandAll[Relab[MomRules[vertices[[i]]], feyn[[i,2]]]]
			},
		{i,1,Length[feyn]}];//Timing
ibFeyn = ibBCJ;


(* ::Subsection::Closed:: *)
(*Check Gravity ( Double Copy: BCJ-on-the-Cuts * Feynman )*)


(******************************************)
(*Choose cut*)
(******************************************)
verts = cutHVverts;
loopEdges = ChooseLoopMom[verts][[1]];
GetCutIsos[verts,loopEdges,ibBCJ,True];
(******************************************)

(********Generate cut numerics*******)
Print["Generating cut numerics..."];
d = 6;
dPhys = 4;
point = CutNumerics[verts,loopEdges,d,dPhys];//Timing
SetVec[point];
helList = {"m","m","p","p"};
GenPol[helList];
helList2 = {"m","m","p","p"};
GenPol2[helList2];

(********Get cut by sewing trees******)
Print["Getting numerical cut by sewing trees..."];
cutSew = SewGDA[verts,loopEdges];//Timing

(********Get cut by cutting BCJ-on-the-Cuts * Feynman *********)
Print["Getting numerical cut by cutting BCJ-on-the-Cuts * Feynman..."];
(* Get rid of fermion contribution for nonSUSY test *)
grav = Table[{(bcj[[i,1]]*(((1/I^2)*feyn[[i,1]]/.Df->0)/.secondCopyReps))/.Ds->d, bcj[[i,2]]}, {i,1,Length[bcj]}];
{cutBCJNum,cutBCJDenom} = CutGravity[verts,loopEdges,grav,ibBCJ];//Timing
cutBCJN = cutBCJNum/cutBCJDenom;
Clear[cutBCJNum,cutBCJDenom];
vars = GetVars[cutBCJN,{Ds}];
NReps = Dispatch[#->NDot[#]&/@vars];
cutBCJN = cutBCJN/.NReps;//Timing

cs = cutSew;
(*
The factor of I^3 comes from:
	1) We factored out an I^L from the YM numerators (See Eq. 2.1 of http://arxiv.org/pdf/1303.6605.pdf)
	2) M ~ I^(L+1) and YM ~ I^L (See Eq. 2.5 of http://arxiv.org/pdf/1303.6605.pdf)
*)
cb = I^3*Chop[Expand[cutBCJN]];

Print["Numerical sewing matches cutting BCJ-on-the-Cuts * Feynman ?"];
Round[cb/cs] === 1


(* ::Subsection::Closed:: *)
(*Check N = 4 SUGRA (Double Copy: BCJ-on-the-Cuts * N = 4 Feynman  &  N = 4 BCJ * N = 0 Feynman )*)


(* N = 4 BCJ representation *)
bcjN4 = Table[{0,bcj[[i,2]]},{i,1,Length[bcj]}];
n4 = Expand[MomConsN[4,s[1,2]^2 s[2,3] Atree[4]]]; (* s^2 t Atree *)
bcjN4[[1,1]] = n4; (* planar double box *)
bcjN4[[3,1]] = n4; (* nonplanar double box *)

(* BCJ-on-the-Cuts * N = 4 Feynman *)
N4Grav1 = Table[{
			(bcj[[i,1]]/.Ds->4) * ((feyn[[i,1]]/.secondCopyReps)/.{Ds->10,Df->8}),
			bcj[[i,2]]},
		{i,1,Length[bcj]}];
(* N = 4 BCJ * N = 0 Feynman *)
N4Grav2 = Table[{
			(bcjN4[[i,1]]/.secondCopyReps) * (feyn[[i,1]]/.{Ds->4,Df->0}),
			bcjN4[[i,2]]},
		{i,1,Length[bcjN4]}];


(******************************************)
(*Choose cut*)
(******************************************)
verts = cutHVverts;
loopEdges = ChooseLoopMom[verts][[1]];
GetCutIsos[verts,loopEdges,ibBCJ,True];
(******************************************)

(********Cutting N=0 BCJ * N=4 Feynman*********)
Print["Cutting N=0 BCJ * N=4 Feynman..."];
{cutNum1,cutDenom1} = CutGravity[verts,loopEdges,N4Grav1,ibBCJ];//Timing
cut1 = cutNum1/cutDenom1;

(********Cutting N=4 BCJ * N=0 Feynman*********)
Print["Cutting N=4 BCJ * N=0 Feynman..."];
{cutNum2,cutDenom2} = CutGravity[verts,loopEdges,N4Grav2,ibBCJ];//Timing
cut2 = cutNum2/cutDenom2;

vars = Union[GetVars[cut1,{Ds,Df}],GetVars[cut2,{Ds,Df}]];//Timing
NReps = Dispatch[
			Table[
				vars[[i]]->((-1)^RandomInteger[])*RandomPrime[1000000](*RandomComplex[WorkingPrecision->100]*)
			,{i,1,Length[vars]}]
		];
Ncut1 = Expand[cut1/.NReps];//Timing
Ncut2 = Expand[cut2/.NReps];//Timing

Print["Different prescriptions give same N = 4 SUGRA cut?"];
Ncut1 === Ncut2

(*
(*Analytic Check*)
Print["Checking equality..."];
ratio = (cutDenom1/cutDenom2)
bool = (Expand[cutNum1 - ratio*(cutNum2)]===0);//Timing
Print["Same analytic cut?"];
bool
*)


(* ::Subtitle::Closed:: *)
(**********************************************************)
(*One-Loop BCJ: Integration Check*)
(**********************************************************)


(* We checked our cutting technology at two loops above, but let's integrate one loop as a solid check. *)
(* We take the expressions from the ancillary file of http://arxiv.org/abs/1309.7416 *)
(* Download "numerators.m" and rename it "numeratorsBCJ1Loop.m" *)


(* ::Section::Closed:: *)
(*1-Loop: YM^2 Integration Benchmark *)
(*(Four Dilaton UV Divergence and Vanishing Four Graviton Divergence)*)


(* ::Subsection::Closed:: *)
(*Projectors & Propagators (Run me)*)


(****************************
Note:
Need to be careful using the traceless part of the projector.
Need to use physical state projectors to get the proper Ds-2 counting.
****************************)

(*Graviton state projector... sort of*)
Clear[projFull];
projFull[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(\[Eta][mu1,mu3]\[Eta][mu2,mu4]+\[Eta][mu1,mu4]\[Eta][mu2,mu3])
			-(1/(Ds-2))\[Eta][mu1,mu2]\[Eta][mu3,mu4]
			);

(*symmetrizing*)
Clear[projSym];
projSym[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(\[Eta][mu1,mu3]\[Eta][mu2,mu4]+\[Eta][mu1,mu4]\[Eta][mu2,mu3])
			);
(*anti-symmetrizing*)
Clear[projAntiSym];
projAntiSym[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(\[Eta][mu1,mu3]\[Eta][mu2,mu4]-\[Eta][mu1,mu4]\[Eta][mu2,mu3])
			);

(*trace ?? Ds-2??*)
Clear[projTrace];
projTrace[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/(Ds-2))\[Eta][mu1,mu2]\[Eta][mu3,mu4]
			);

(*Anti-sym. three-form*)
Clear[proj3];
proj3[{mu1a_,mu1b_,mu1c_},{mu2a_,mu2b_,mu2c_}]:= Module[{side1,side2,perms,sigs},
	side1 = {mu1a,mu1b,mu1c};
	side2 = {mu2a,mu2b,mu2c};
	perms = Permutations[side2];
	sigs = Signature[side2]*(Signature/@perms);
	
	1/Length[perms]*
	Plus@@
		(
			sigs
			*Times@@@(
				Table[
					\[Eta][side1[[i]],#[[i]]]
					,{i,1,Length[side1]}
				]&/@perms
			)
		)
];


(*Checked factors of -1 and I*)
(*Graviton propagator*)
Clear[hProp];
hProp[i1_,i2_]:=I*projFull[{\[Mu][i1,1],\[Mu][i1,2]},{\[Mu][i2,1],\[Mu][i2,2]}]/R[p[i1]];
(*Scalar propagator*)
Clear[\[Phi]Prop];
\[Phi]Prop[i1_,i2_]:=I/R[p[i1]];
(*Antisymmetric propagator*)
Clear[BProp];
BProp[i1_,i2_]:=I*projAntiSym[{\[Mu][i1,1],\[Mu][i1,2]},{\[Mu][i2,1],\[Mu][i2,2]}]/R[p[i1]];
(*Guess for anti-symmetric 3-form propagator*)
Clear[threeFormProp];
threeFormProp[i1_,i2_]:=I*proj3[{\[Mu][i1,1],\[Mu][i1,2],\[Mu][i1,3]},{\[Mu][i2,1],\[Mu][i2,2],\[Mu][i2,3]}]/R[p[i1]];


Clear[projGDATemp];
Clear[projGDA];
projGDATemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			H[mu1,mu3]H[mu2,mu4]
			);
projGDA[i_,{L1_,R1_},{L2_,R2_}]:=projGDATemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

Clear[projGTemp];
Clear[projG];
projGTemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(H[mu1,mu3]H[mu2,mu4]+H[mu1,mu4]H[mu2,mu3])
			-(1/(Ds-2))H[mu1,mu2]H[mu3,mu4]
			);
projG[i_,{L1_,R1_},{L2_,R2_}]:=projGTemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

Clear[projGDTemp];
Clear[projGD];
projGDTemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(H[mu1,mu3]H[mu2,mu4]+H[mu1,mu4]H[mu2,mu3])
			);
projGD[i_,{L1_,R1_},{L2_,R2_}]:=projGDTemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

Clear[projGATemp];
Clear[projGA];
projGATemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			H[mu1,mu3]H[mu2,mu4]
			-(1/(Ds-2))H[mu1,mu2]H[mu3,mu4]
			);
projGA[i_,{L1_,R1_},{L2_,R2_}]:=projGATemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

Clear[projATemp];
Clear[projA];
projATemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/2)(H[mu1,mu3]H[mu2,mu4]-H[mu1,mu4]H[mu2,mu3])
			);
projA[i_,{L1_,R1_},{L2_,R2_}]:=projATemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

Clear[projDTemp];
Clear[projD];
projDTemp[{mu1_,mu2_},{mu3_,mu4_}]:=(
			(1/(Ds-2))H[mu1,mu2]H[mu3,mu4]
			);
projD[i_,{L1_,R1_},{L2_,R2_}]:=projDTemp[{\[Mu][i,L1],\[Mu][i,R1]},{\[Mu][i,L2],\[Mu][i,R2]}];

(*Guess for 3-form*)
Clear[proj3FormTemp];
Clear[proj3Form];
proj3FormTemp[{mu1a_,mu1b_,mu1c_},{mu2a_,mu2b_,mu2c_}]:=(
			proj3[{mu1a,mu1b,mu1c},{mu2a,mu2b,mu2c}]/.\[Eta]->H
			);
proj3Form[i_,{side1a_,side1b_,side1c_},{side2a_,side2b_,side2c_}]:=proj3FormTemp[{\[Mu][i,side1a],\[Mu][i,side1b],\[Mu][i,side1c]},{\[Mu][i,side2a],\[Mu][i,side2b],\[Mu][i,side2c]}];

(*Guess for 3-form*)
Clear[proj3FormTemp0];
Clear[proj3Form0];
proj3FormTemp0[{mu1a_,mu1b_,mu1c_},{mu2a_,mu2b_,mu2c_}]:=(
			proj3[{mu1a,mu1b,mu1c},{mu2a,mu2b,mu2c}]/.\[Eta]->h0
			);
proj3Form0[i_,{side1a_,side1b_,side1c_},{side2a_,side2b_,side2c_}]:=proj3FormTemp0[{\[Mu][i,side1a],\[Mu][i,side1b],\[Mu][i,side1c]},{\[Mu][i,side2a],\[Mu][i,side2b],\[Mu][i,side2c]}];


(********************************)
(*Apply projectors*)
(********************************)
Clear[ApplyProj];
ApplyProj[expr0_,projList_List]:=Module[{expr,stripRep,stripDispatch,indexLabelL,indexLabelR,epProj,vars,reps},

If[Length[projList]!=4,"Need four projector choices",

	expr = expr0;

	indexLabelL = 3;
	indexLabelR = 4;
	
	stripRep={
		ke[i_,j_]:>p[i][\[Mu][j,indexLabelL]],
		ee[i_,j_]:>\[Eta][\[Mu][i,indexLabelL],\[Mu][j,indexLabelL]],
		kE[i_,j_]:>p[i][\[Mu][j,indexLabelR]],
		EE[i_,j_]:>\[Eta][\[Mu][i,indexLabelR],\[Mu][j,indexLabelR]],
		eE[i_,j_]:>\[Eta][\[Mu][i,indexLabelL],\[Mu][j,indexLabelR]]
	};

	epProj[iLeg_,particleName_]:=
		Simplify[
			Propagate[
				\[Epsilon][iLeg][\[Mu][iLeg,1]]
				*\[Epsilon]2[iLeg][\[Mu][iLeg,2]]
				*ToExpression["proj"<>particleName][iLeg,{1,2},{indexLabelL,indexLabelR}]
			]/.{ke[i_,i_]:>0,kE[i_,i_]:>0,kk[i_,i_]:>0(*,p[iLeg][x__]:>0*)}(*CANNOT kill p^mu's even if amp is gauge inv. See soft paper*)
		];
	
	vars = GetVars[expr,{Ds,Dsm2,kk[__]}];
	stripDispatch = Dispatch[Thread[vars -> (vars/.stripRep)]];

	Print["Strip ext polarization"];
	Print[expr = expr /. stripDispatch;//Timing];

	reps =
		{
			ke[ref1,i_]:>0,kE[ref1,i_]:>0,
			ke[ref2,i_]:>0,kE[ref2,i_]:>0,
			ke[ref3,i_]:>0,kE[ref3,i_]:>0,
			ke[ref4,i_]:>0,kE[ref4,i_]:>0,
			kk[ref1,2]->0,kk[ref1,3]->0,
			kk[ref2,1]->0,kk[ref2,3]->0,
			kk[ref3,1]->0,kk[ref3,2]->0,
			kk[ref1,ref2]->0,kk[ref1,ref3]->0,
			kk[ref1,ref4]->0,kk[ref2,ref3]->0,
			kk[ref2,ref4]->0,kk[ref3,ref4]->0
		};
	
	Print["Begin applying projectors"];
	Print[expr = (Propagate[(epProj[1,projList[[1]]]*expr)/.{Ds->Dsm2+2}]/.{Ds->Dsm2+2})/.reps;//Timing];
	Print[expr = Expand[MomConsN[4,expr]]/.reps;//Timing];
	Print[expr = (Propagate[(epProj[2,projList[[2]]]*expr)/.{Ds->Dsm2+2}]/.{Ds->Dsm2+2})/.reps;//Timing];
	Print[expr = Expand[MomConsN[4,expr]]/.reps;//Timing];
	Print[expr = (Propagate[(epProj[3,projList[[3]]]*expr)/.{Ds->Dsm2+2}]/.{Ds->Dsm2+2})/.reps;//Timing];
	Print[expr = Expand[MomConsN[4,expr]]/.reps;//Timing];
	Print[expr = (Propagate[(epProj[4,projList[[4]]]*expr)/.{Ds->Dsm2+2}]/.{Ds->Dsm2+2})/.reps;//Timing];
	Print[expr = Expand[MomConsN[4,expr]]/.reps;//Timing];

	Print["Togethering"];
	Print[expr = Together[expr];//Timing];

	Print["Gauge invariance check:"];
	Print[Table[Expand[MomConsN[4,Relab[{\[Epsilon][i]->p[i]},expr]]],{i,1,4}]];
	Print[Table[Expand[MomConsN[4,Relab[{\[Epsilon]2[i]->p[i]},expr]]],{i,1,4}]];

	expr
]
];
(*******************************)
(*******************************)

ApplyProj[expr0_,theoryID_String]:=ApplyProj[expr0,{theoryID,theoryID,theoryID,theoryID}]


processIDList = {"hhhh","hhhD","hhDD","hDDD","DDDD","hhBB","BBBB","hBBD","BBDD","hBBB","hhhB","hhBD","hBDD","BDDD","BBBD"};

applyProjArg[processID_] := Switch[processID,
				"hhhh",
					{"G","G","G","G"},
				"hhhD",
					{"G","G","G","D"},
				"hhDD",
					{"G","G","D","D"},
				"hDDD",
					{"G","D","D","D"},
				"DDDD",
					{"D","D","D","D"},
				"hhBB",
					{"G","G","A","A"},
				"BBBB",
					{"A","A","A","A"},
				"hBBD",
					{"G","A","A","D"},
				"BBDD",
					{"A","A","D","D"},
				"hBBB",
					{"G","A","A","A"},
				"hhhB",
					{"G","G","G","A"},
				"hhBD",
					{"G","G","A","D"},
				"hBDD",
					{"G","A","D","D"},
				"BDDD",
					{"A","D","D","D"},
				"BBBD",
					{"A","A","A","D"}
			]


(* ::Subsection::Closed:: *)
(*Load Numerators (Run me)*)


vertices1loop = { 
	{{-1},{-2},{-3},{-4},{1,6,-5},{2,7,-6},{3,8,-7},{4,5,-8}},
	{{-1},{-2},{-3},{-4},{1,2,6},{-5,-6,7},{3,8,-7},{4,5,-8}},
	{{-1},{-2},{-3},{-4},{1,2,6},{-5,-6,7},{5,-7,8},{3,4,-8}}
};

SetDirectory[NotebookDirectory[]];
Get["numeratorsBCJ1Loop.m"];
nums = Expand[MomConsN[4,numerators[[{1,4,6}]]/.{Df->0,Dg->Ds-2,Nf->0,Ng->1,Ns->0,R[x_]:>AngleBracket[x/.k->p,x/.k->p]}]];
nums = Table[I^(1+1) nums[[i]]*(nums[[i]]/.{ke->kE,le->lE,ee->EE}),{i,1,Length[nums]}];
props = Factor[Expand[MomConsN[4,List@@@(propagators[[{1,4,6}]])/.{Df->0,Dg->Ds-2,Nf->0,Ng->1,Ns->0,R[x_]:>AngleBracket[x/.k->p,x/.k->p]}]]];


(* ::Subsection::Closed:: *)
(*UV*)


(* Expand to vacuum integrals, extract UV, tensor reduce, sum over permutations *)
(* (I'm sure there a cleaner method. This is quick and dirty) *)


integrands = {0,0,0}; (*box, triangle, bubble*)
d = 4;

For[i=1, i<=3, i++,
	Print["Working on diagram "<>ToString[i]];
	powReps = {
		lk[i_,j_]:>\[Alpha][i]*lk[i,j],
		ll[i_,i_]:>\[Alpha][i]^2*ll[i,i],
		R[l[i_]]:>\[Alpha][i]^2R[l[i]],
		InvR[i_]:>(1/\[Alpha][i])^2 InvR[i],
		le[i_,j_]:>\[Alpha][i]*le[i,j],
		lE[i_,j_]:>\[Alpha][i]*lE[i,j]};
	invPowReps = {
		lk[i_,j_]:>(1/\[Beta][i])*lk[i,j],
		ll[i_,i_]:>(1/\[Beta][i])^2*ll[i,i],
		R[l[i_]]:>(1/\[Beta][i])^2 R[l[i]],
		InvR[i_]:>(\[Beta][i])^2 InvR[i],
		le[i_,j_]:>(1/\[Beta][i])*le[i,j],
		lE[i_,j_]:>(1/\[Beta][i])*lE[i,j]};
	pr = props[[i]] /. {ll[5,5]:>1/InvR[5]};
	Print[nu = Expand[nums[[i]]]/.powReps;//Timing];
	\[Alpha]s = DeleteDuplicates[Cases[nu,\[Alpha][_],\[Infinity]]];

	maxPow = 10; (*Should be 8 but set to 10 to be safe*)
	Print[numPows = Expand[Table[1/i! D[nu,{\[Alpha][5],i}]/.\[Alpha][5]->0,{i,0,maxPow}]];//Timing];

	degOfDiv[d_, pow_] := (d + pow) - 2*Length[Position[pr,InvR[5]]];

	maxDegOfDiv = degOfDiv[d,maxPow];
	denomEx =(InvR[5]^Length[Position[pr,InvR[5]]])*Series[1/(InvR[5]^Length[Position[pr,InvR[5]]] * Times@@pr) /. invPowReps,{\[Beta][5],0,maxDegOfDiv}];
	numEx = Table[
		dod = degOfDiv[d,pow];
		If[dod < 0, 0,
			Expand[SeriesCoefficient[denomEx,dod]]/.{\[Beta][5]->1}
		]
	,{pow,0,maxPow}];

	(* +1 because indexing starts at one but pow starts at 0 *)
	integrand = Sum[numEx[[pow+1]]*numPows[[pow+1]],{pow,0,maxPow}];
	Print[integrand = TensorReduce[integrand];//Timing];
	integrand = DeleteCases[integrand,{_,0}];
	Print[integrand = Expand[integrand /. {ll[5,5] -> 1/InvR[5] (*+ m^2*)}];//Timing];
	integrand = Plus@@(Times@@@integrand);
	Print[integrand = Expand[MomConsN[4,integrand/.Ds->d]];//Timing];

	graph = vertices1loop[[i]];
	permReps = Join[Thread[p/@{1,2,3,4}->p/@#],Thread[\[Epsilon]/@{1,2,3,4}->\[Epsilon]/@#],Thread[\[Epsilon]2/@{1,2,3,4}->\[Epsilon]2/@#]]&/@Permutations[{1,2,3,4}];
	symFactor = Length[Transpose[EdgewiseAutomorphisms[graph]]];
	Print[integrand = (1/symFactor) * Sum[MomConsN[4,Relab[permReps[[i]],integrand]],{i,1,Length[permReps]}];//Timing];
	Print[integrand = Expand[integrand];//Timing];

	integrands[[i]] = integrand;
];

uv = ExpandAll[Plus@@integrands];


(********Generate cut numerics*******)
Print["Generating cut numerics..."];
d = 4;
dPhys = 4;
point = CutNumerics[{{-1},{-2},{-3},{-4},{1,2,3,4}},{}, d, dPhys];//Timing
SetVec[point];
helList={"m","m","p","p"};
GenPol[helList];
GenPol2[helList]; (*Gravitons have same helicity on both sides of double copy*)

NDot[uv] (*4 gravitons should vanish*)


uvDDDD = ApplyProj[uv,applyProjArg["DDDD"]];//Timing


uv1 = Factor[uvDDDD/.Dsm2->d-2]/.eE[_,_]:>Sqrt[2]
(*From Eq. (4.3) of http://arxiv.org/pdf/1303.6605.pdf*)
uv2 = (3d^2 - 92d + 1132)/120 * Factor[MomConsN[4,s[1,2]^2+s[2,3]^2+s[1,3]^2]] * InvR[5]^2

Print["Agrees with literature?"];
uv1 === uv2


(* ::Subtitle::Closed:: *)
(**********************************************************)
(*NLSM Examples*)
(**********************************************************)


(* Check our NLSM Feynman rules *)
(* We still need to fix factors of -1 and I ! *)


(* Be careful with 8 points. The color information takes a long time to generate if the result isn't already (automatically) saved *)
(* Check the "aux" -> "colorTree" folders to see if the 8-point information is already there *)
(* I didn't upload 8-point information to GitHub because of file size constraints *)

(* There are also better ways of finding the diagrams that contribute to a color-ordered amplitude; *)
(* however, that hasn't been implemented in MasterFunctions yet. *)
(* See comments in section "fabc Color Basis" in MasterFunctions *)


(* ::Section::Closed:: *)
(*Tree-Level Check*)


(* Note -1 and I fudge factors *)
(* We compare with the results of arXiv: 1304.3048 *)


(* ::Subsection::Closed:: *)
(*4-Point Color-Ordered Amplitude*)
(*( Compare with arXiv:1304.3048, Eq. (2.33) )*)


n = 4;
fudge = I;
me = Expand[fudge * NLSMtree[n]/.F->1];
correct = Expand[ExpandAll[Plus@@(MomConsN[n,Relab[#,
	(
		(1/2)s[1,2]
	)
],True]&/@(Table[p[j]->p[#[[j]]],{j,1,n}]&/@Table[RotateLeft[Range[n],i],{i,0,n-1}]))]];

me == correct


(* ::Subsection::Closed:: *)
(*6-Point Color-Ordered Amplitude*)
(*( Compare with arXiv: 1304.3048, Eq. (2.35) )*)


n = 6;
fudge = -I;
me = Expand[fudge * NLSMtree[n]/.F->1];
correct = Expand[ExpandAll[Plus@@(MomConsN[n,Relab[#,
	(
		-(1/2) * (s[1,2] + s[2,3]) * (s[1,2,3,4] + s[4,5]) / s[1,2,3] + s[1,2]
	)
],True]&/@(Table[p[j]->p[#[[j]]],{j,1,n}]&/@Table[RotateLeft[Range[n],i],{i,0,n-1}]))]];

me == correct


(* ::Subsection::Closed:: *)
(*8-Point Color-Ordered Amplitude*)
(*( Compare with arXiv: 1304.3048, Eq. (B.1) )*)


n = 8;
fudge = I;
me = Expand[fudge * NLSMtree[n]/.F->1];
correct = Expand[ExpandAll[Plus@@(MomConsN[8,Relab[#,
	(
		(1/2)*((s[1,2]+s[2,3])(s[1,2,3,4]+s[4,5,6,7])(s[5,6]+s[6,7]))/(s[1,2,3]s[5,6,7])
		+((s[1,2]+s[2,3])(s[1,2,3,4]+s[4,5])(s[6,7]+s[7,8]))/(s[1,2,3]s[6,7,8])
		-((s[1,2]+s[2,3])(s[4,5]+s[4,5,6,7]+s[5,6]+s[5,6,7,8]+s[6,7]+s[7,8]))/s[1,2,3]
		+2*s[1,2] + (1/2)s[1,2,3,4]
	)
],True]&/@(Table[p[j]->p[#[[j]]],{j,1,n}]&/@Table[RotateLeft[Range[n],i],{i,0,n-1}]))]];

me == correct


(* ::Section::Closed:: *)
(*Tree Amplitude Generation from Yang-Mills*)


(* 
Inspired by http://arxiv.org/abs/1412.3479, we give a prescription for
obtaining NLSM amplitudes from Yang-Mills amplitudes.

Prescription:
	1. Atree = YM amplitude
	2. Set ke[i,j]'s in Atree to zero
	3. Choose the coefficient of ANY ee[i,j]
	4. Set the remaining ee[i,j]'s to s[i,j]'s

This is like computing Yang-Mills coupled to external scalars of different
flavors, where ee[i,j] should be thought of as a Kronecker delta in flavor.
We privilege two external flavors (by picking the coefficient of any ee[i,j]).
Then, we exploit a sort of flavor-kinematics duality by setting the remaining
ee[i,j]'s to s[i,j]'s
*)


n = 6;
fudge = (-1)^(n/2);
nlsm = Expand[fudge * NLSMtree[n]/.F->1];//Timing
(* This can take a long time to load (and very long time to compute) for 8 points *)
ym = Atree[n];//Timing 

ym = Expand[ym /. ke[__]:>0];//Timing
ees = GetVars[ym,{s[__]}];
nlsmList = Expand[MomConsN[n, Coefficient[ym, #] /. ee[i_,j_]:>s[i,j], True]] &/@ ees;//Timing

# === nlsm &/@ nlsmList


(* ::Section::Closed:: *)
(*1-Loop UV Integration: DBI ~ NLSM x N = 4 *)


vertices1loop = { 
	{{-1},{-2},{-3},{-4},{1,6,-5},{2,7,-6},{3,8,-7},{4,5,-8}},
	{{-1},{-2},{-3},{-4},{1,2,6},{-5,-6,7},{3,8,-7},{4,5,-8}},
	{{-1},{-2},{-3},{-4},{1,2,6},{-5,-6,7},{5,-7,8},{3,4,-8}}
};

n = 4;
loops = 1;

graph = vertices1loop[[1]];
{num,props} = FeynNLSM[graph];//Timing
num = (I^loops)*(-1)^(n/2)(-I)num; (*Fudge to match Henrik*)

numNLSMBox = Factor[num/.{s[i_,j_]:>2kk[i,j],kk[5,5]->ll[5,5],kk[i_,5]:>lk[5,i]}];
propsNLSMBox = Expand[MomConsN[4,Times@@(props/.R[x__]:>AngleBracket[x,x])]];
propsNLSMBox = Factor[propsNLSMBox/.{s[i_,j_]:>2kk[i,j],kk[5,5]->ll[5,5],kk[i_,5]:>lk[5,i]}];

tree4 = Factor[MomConsN[4,Atree[4],True]];
numN4Box = Factor[I*s[1,2]s[2,3]*tree4];

dc = Factor[numNLSMBox * numN4Box];


powReps = {
	lk[i_,j_]:>\[Alpha][i]*lk[i,j],
	ll[i_,i_]:>\[Alpha][i]^2*ll[i,i],
	R[l[i_]]:>\[Alpha][i]^2R[l[i]],
	InvR[i_]:>(1/\[Alpha][i])^2 InvR[i],
	le[i_,j_]:>\[Alpha][i]*le[i,j]};
invPowReps = {
	lk[i_,j_]:>(1/\[Beta][i])*lk[i,j],
	ll[i_,i_]:>(1/\[Beta][i])^2*ll[i,i],
	R[l[i_]]:>(1/\[Beta][i])^2 R[l[i]],
	InvR[i_]:>(\[Beta][i])^2 InvR[i],
	le[i_,j_]:>(1/\[Beta][i])*le[i,j]};
pr = propsNLSMBox /. {ll[5,5]:>1/InvR[5]};
pre = 1/(45*F^4);
nu = Expand[numNLSMBox / pre]/.powReps;

maxPow = 8;
numPows = Expand[Table[1/i! D[nu,{\[Alpha][5],i}]/.\[Alpha][5]->0,{i,0,maxPow}]];

degOfDiv[d_, pow_] := (d + pow) - 2*Length[Position[pr,InvR[5]]];


uvList = {0,0,0,0,0,0};

For[d = 2, d <= 12, d+=2,
	maxDegOfDiv = degOfDiv[d,maxPow];
	denomEx =(InvR[5]^Length[Position[pr,InvR[5]]])*Series[1/(InvR[5]^Length[Position[pr,InvR[5]]] * Times@@pr) /. invPowReps,{\[Beta][5],0,maxDegOfDiv}];
	numEx = Table[
		dod = degOfDiv[d,pow];
		If[dod < 0, 0,
			Expand[SeriesCoefficient[denomEx,dod]]/.{\[Beta][5]->1}
		]
	,{pow,0,maxPow}];

	(* +1 because indexing starts at one but pow starts at 0 *)
	integrand = Sum[numEx[[pow+1]]*numPows[[pow+1]],{pow,0,maxPow}];
	Print[integrand = TensorReduce[integrand];//Timing];
	integrand = DeleteCases[integrand,{_,0}];
	Print[integrand = Expand[integrand /. {ll[5,5] -> 1/InvR[5] (*+ m^2*)}];//Timing];
	integrand = Plus@@(Times@@@integrand);
	Print[integrand = Expand[MomConsN[4,integrand/.Ds->d]];//Timing];

	permReps = Thread[p/@{1,2,3,4}->p/@#]&/@Permutations[{1,2,3,4}];
	symFactor = Length[Transpose[EdgewiseAutomorphisms[graph]]];
	uv = Factor[ pre * (1/symFactor) * Sum[MomConsN[4,Relab[permReps[[i]],integrand]],{i,1,Length[permReps]}]];

	uvList[[d/2]] = uv;
];


Print["Cleaned up expressions for UV divergence in D = 2,4,6,8,10,12 dimensions:"];
uvList2 =
{
	0,
	+(1/2) * (s^2 + t^2 + u^2) * InvR[5]^2,
	-(1/2) * (s * t * u) * InvR[5]^3,
	+(1/40) * (s^2 + t^2 + u^2)^2 * InvR[5]^4,
	-(1/28) * (s * t * u)*(s^2 + t^2 + u^2) * InvR[5]^5,
	+(1/1008) * ((s^2 + t^2 + u^2)^3 + 12*(s * t * u)^2) * InvR[5]^6
}*(1/F^4)

Print["Matches our result?"];
Factor[MomConsN[4,uvList2/.{s->s[1,2],t->s[2,3],u->s[1,3]}]] ===  uvList


(* ::Subtitle::Closed:: *)
(**********************************************************)
(*Soft Gravity Dilaton Example*)
(**********************************************************)


(* Show that leading and subleading orders vanish for the *)
(* soft gravity dilaton as show in arXiv: 1512.03316 *)

(* Note that CutNumericsSoft always sets leg 2 to be soft. Quick and dirty adjustment... but it works! *)


tree4 = Mtree[4];
tree5 = Mtree[5];

scaling = 
Table[

verts = {{-1},{-2},{-3},{-4},{-5},{1,2,3,4,5}};
dim = 4;
dimPhys = 4;
softScaling = 10^(-i);
numPoint = CutNumericsSoft[verts, {}, dim, dimPhys, softScaling];
SetVec[numPoint];

GenPol[{"m","m","p","p","m"}];
GenPol2[{"p","p","p","p","m"}];
treeN1 = NDot[tree5]/.times->Times;

GenPol[{"p","m","p","p","m"}];
GenPol2[{"m","p","p","p","m"}];
treeN2 = NDot[tree5]/.times->Times;

GenPol[{"m","p","p","p","m"}];
GenPol2[{"p","m","p","p","m"}];
treeN3 = NDot[tree5]/.times->Times;

GenPol[{"p","p","p","p","m"}];
GenPol2[{"m","m","p","p","m"}];
treeN4 = NDot[tree5]/.times->Times;

(1/2)(treeN1 + treeN2 + treeN3 + treeN4),

{i,0,25,5}
];

(* See that scaling goes like q^1 (i.e. softScaling^1), i.e. subsubleading *)
scaling//MatrixForm
Print["At each step, we scale the soft momenta by an extra factor of 10^-5."];
Print["The log_10 of the ratio at each step is approximately: "];
Table[ Round[Re[Log[10,scaling[[i+1]]/scaling[[i]]]]], {i,1,5}]
Print["Because we generate different numerics each time we won't see the exactly correct scaling."];
Print["However, it approximately hovers around 10^-5, which ~(soft momentum)^1, i.e. subleading order."];



