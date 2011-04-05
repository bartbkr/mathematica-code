Print["This is Perturbation AIM Simple Version 2.4.1ets"] ;
(*                                                                          *)
(* Credits: original algorithm designed by Gary Anderson, 2002-3;           *)
(*   corrected, optimized, and recoded by Eric Swanson, 2004.               *)
(* Code is currently maintained by Eric Swanson, see                        *)
(*   http://www.ericswanson.pro for most recent version and for             *)
(*   background, instructions, capabilities, examples, and tips.            *)
(* This code will only work with Mathematica version 5.  While it is not    *)
(*   too hard to modify the code to work in earlier versions of             *)
(*   Mathematica, those earlier versions are *much* slower at solving       *)
(*   linear systems of equations, and a few of the function calls are also  *)
(*   a bit more awkward.                                                    *)
(* This is the "simple" version of the code.  A more complicated but        *)
(*   faster and more numerically stable version of the code is available    *)
(*   for download from the web site above.  This "simple" version should    *)
(*   be completely adequate for most purposes, however, and will be much    *)
(*   easier for the user to understand and even tailor to his/her           *)
(*   individual needs.                                                      *)
(*                                                                          *)

SetOptions["stdout",PageWidth->79] ;             (* adjust this as you like *)

(* Load two standard Mathematica packages: *)

<<DiscreteMath`Combinatorica`        (* we'll use the Compositions function *)
<<LinearAlgebra`MatrixManipulation`   (* we'll use the BlockMatrix function *)


(* The AIMZeroTol parameter should be thought of as an "economic zero."     *)
(*   In other words, AIMZeroTol should be set to the smallest possible      *)
(*   value that a coefficient in the model (either in inputs or outputs)    *)
(*   could take on and still have economic meaning (as opposed to being     *)
(*   simply an artifact of finite-precision machine arithmetic).            *)
(* The default value of AIMZeroTol is 10^-10, but you can change it and     *)
(*   in fact should play around with it for any given model to test the     *)
(*   numerical stability of the solution to the model.  Note that when      *)
(*   working with arbitrary-precision numbers in Mathematica, you have      *)
(*   more freedom to set this parameter to smaller values, but you should   *)
(*   still think of it as being an "economic zero" rather than a "machine   *)
(*   numerical zero" (the latter is left to the internals of Mathematica    *)
(*   to handle appropriately).  Thus, values like 10^-50 are probably       *)
(*   unnecessarily small and will only require you to input and compute     *)
(*   more digits of numerical precision than is really worthwhile.          *)

AIMZeroTol=10^-10 ;


(* Useful utilities that we will repeatedly call.  The syntax               *)
(*  FunctionName[vars]:= FunctionName[vars] = ...                           *)
(*  tells Mathematica to remember the answer, so that it does not have      *)
(*  to be recomputed every time the function is called.                     *)

AIMGenericArgs[eqns_]:= AIMGenericArgs[eqns] =
Table[Unique["arg"],{AIMNArgs[eqns]}] ;

AIMLagVars[eqns_]:= AIMLagVars[eqns] =
Flatten[Map[Table[#[t+i],{i,Min[Cases[eqns,#[t+j_.]->j,Infinity]],-1}]&,
							AIMVarNames[eqns]]] ;
AIMMaxLag[eqns_]:= AIMMaxLag[eqns] =
-Min[Cases[eqns,x_Symbol[t+i_.]->i,Infinity]] ;

AIMMaxLead[eqns_]:= AIMMaxLead[eqns] =
Max[Cases[eqns,x_Symbol[t+i_.]->i,Infinity]] ;

AIMNArgs[eqns_]:= AIMNArgs[eqns] = Length[AIMStateVars[eqns]] ;

AIMNEqns[eqns_]:= AIMNEqns[eqns] = Length[Flatten[{eqns}]] ;

AIMNVars[eqns_]:= AIMNVars[eqns] = Length[AIMVarNames[eqns]] ;

AIMShocks[eqns_]:= AIMShocks[eqns] = Union[Cases[eqns,eps[x_][t],Infinity]] ;

AIMSSSubs={Sigma->0, eps[_][_]->0, x_[t+_.]:>Symbol[SymbolName[x]<>"AIMSS"]}

AIMSSVars[eqns_]:= AIMSSVars[eqns] = 
Map[Symbol[SymbolName[#]<>"AIMSS"]&, AIMVarNames[eqns]] ;

AIMStateVars[eqns_]:= AIMStateVars[eqns] =
Join[AIMLagVars[eqns], AIMShocks[eqns], {Sigma}] ;

AIMVarNames[eqns_]:= AIMVarNames[eqns] =
Union[Cases[eqns,x_Symbol[t+i_.]->x,Infinity]] ;


(* Calculate the steady state *)

AIMSS[eqns_,___,opts___Rule]:= AIMSS[eqns,AIMZeroTol] =
Module[{ssguess,precision,sseqns,symbols,ss},
 ssguess = AIMSSGuess /.{opts} /.AIMSSGuess->Table[0,{AIMNVars[eqns]}] ;
 precision = AIMPrecision /.{opts} /.AIMPrecision->Precision[eqns]-2 ;
 sseqns = Thread[(eqns /.Equal->Subtract /.AIMSSSubs)==0] ;
 AIMModelDiagnostics[eqns] ;
 Print["Finding steady state, AIMSSGuess->",InputForm[ssguess]] ;
 symbols = Complement[Union[Cases[sseqns,_Symbol,Infinity]],
						Join[AIMSSVars[eqns],{E}]] ;
 If[Length[symbols]>0,Print["Warning: found symbols ",symbols," in equations"]];
 ss = Chop[FindRoot[sseqns, Transpose[{AIMSSVars[eqns],ssguess}],
		MaxIterations->1000, WorkingPrecision->precision], AIMZeroTol] ;
 If[Head[ss]===FindRoot, $Failed, ss]
] ;


(* Front end that does error-checking and formats the output *)

AIMSeries[eqns_,deg_Integer]:=
Module[{soln,const,argsubs},
 If[AIMSS[eqns,AIMZeroTol] === $Failed ||
	(soln=AIMSoln[eqns,deg,AIMZeroTol]) === $Failed, $Failed,
 Print["Formatting output, time is ", Date[]] ;
 const = AIMSSVars[eqns] /.AIMSS[eqns,AIMZeroTol] ;
 argsubs = Thread[AIMGenericArgs[eqns]->AIMStateVars[eqns]] ;
 Thread[Through[AIMVarNames[eqns][t]] == const + (soln /.argsubs)]
]] ;


(* Front end for Linear AIM *)

AIMSoln[eqns_,1,zerotol_]:= AIMSoln[eqns,1,zerotol] =
Module[{eqnseq0,allvars,hmat,epsmat,soln,cofb,s0inv,alllagvars},
 eqnseq0 = Flatten[{eqns}] /.Equal->Subtract ;
 allvars = Flatten[Map[Through[AIMVarNames[eqns] [t+#]]&,
			 Range[-Max[AIMMaxLag[eqns],1], AIMMaxLead[eqns]]]] ;
 hmat = Outer[D,eqnseq0,allvars] /.AIMSSSubs /.AIMSS[eqns,zerotol] ;
 epsmat = Outer[D,eqnseq0,AIMShocks[eqns]] /.AIMSSSubs /.AIMSS[eqns,zerotol] ;
 If[(soln=AIMLinearSoln[hmat,AIMMaxLead[eqns]]) === $Failed, $Failed,
 {cofb,s0inv} = soln ;
 alllagvars = allvars[[Range[Max[AIMMaxLag[eqns],1]*AIMNVars[eqns]]]] ;
 Chop[cofb .alllagvars - s0inv .epsmat .AIMShocks[eqns] /.
		Thread[AIMStateVars[eqns]->AIMGenericArgs[eqns]], zerotol]
]] ;


(* This is the heart of the program.  Derivatives are evaluated at steady   *)
(*  state, the expectation is taken, and coefficients are solved using      *)
(*  the method of undetermined coefficients.                                *)
(* The following two tricks are not strictly necessary and bloat the code   *)
(*  a bit, but increase speed and, more importantly, seem to reduce the     *)
(*  likelihood of numerical inaccuracies:                                   *)
(* 1. Solve a certainty-equivalent version of the problem first (Sigma=0).  *)
(* 2. Impose constraint that all terms linear in Sigma (even if nonlinear   *)
(*     in other variables) have zero coefficients, which is fairly easy to  *)
(*     show mathematically.                                                 *)

AIMSoln[eqns_,deg_Integer,zerotol_]:= AIMSoln[eqns,deg,zerotol] =
Module[{args,drvindxs,cedrvindxs,stdrvindxs,cecoeffs,stcoeffs,nextceTerms,
	nextstTerms,bsubs,ssderivs,cesystem,cesoln,dum,dropce,stsystem,stsoln},
 args = AIMGenericArgs[eqns] ;
 drvindxs = Compositions[deg,AIMNArgs[eqns]] ;
  cedrvindxs = Select[drvindxs, #[[-1]]==0 &] ;
  stdrvindxs = Select[drvindxs, #[[-1]] >1 &] ;
 cecoeffs = Table[Unique[],{AIMNVars[eqns]},{Length[cedrvindxs]}] ;
  stcoeffs = Table[Unique[],{AIMNVars[eqns]},{Length[stdrvindxs]}] ;
 nextceTerms = cecoeffs .Map[Apply[Times,Power[args,#]]&, cedrvindxs] ;
  nextstTerms = stcoeffs .Map[Apply[Times,Power[args,#]]&, stdrvindxs] ;
 bsubs = Thread[Map[bFunc,AIMVarNames[eqns]] -> Map[Apply[Function,
	 {args,#}]&, AIMSoln[eqns,deg-1,zerotol] + nextceTerms + nextstTerms]] ;
 Print["Differentiating eqns, time is ", Date[]] ;
 ssderivs = Chop[AIMDerivatives[eqns,deg][0] /.bsubs, zerotol] ;
 Print["Undetermined coefficients to solve: ",
		 AIMNVars[eqns] *(Length[cedrvindxs]+Length[stdrvindxs])] ;
 Print["Calculating CE solution, time is ", Date[]] ;
 cesoln = If[AIMNArgs[eqns]===1, PrependTo[args,dum]; {},
   cesystem = Flatten[Chop[Take[CoefficientArrays[ssderivs /.args[[-1]]->0,
						Drop[args,-1]], -1], zerotol]] ;
   Chop[Flatten[NSolve[cesystem, Flatten[cecoeffs]]], zerotol]] ;
 Print["Calculating Stoch solution, time is ", Date[]] ;
 dropce = Flatten[Drop[CoefficientArrays[ssderivs, args[[-1]]], 2]] ;
 stsystem = Chop[Expand[Flatten[CoefficientArrays[dropce, Drop[args,-1]]]] /.
		 eps[x_][_]^n_->mom[x,n] /.eps[_][_]->0 /.cesoln, zerotol] ;
 stsoln = Chop[Flatten[NSolve[stsystem, Flatten[stcoeffs]]], zerotol] ;
 AIMSoln[eqns,deg-1,zerotol] + (nextceTerms /.cesoln) + (nextstTerms /.stsoln)
] ;

(* That's essentially it.  The following routine calculates derivatives	    *)
(*   of the equations composed with the (unknown) solution functions b.	    *)
(* The trick of using univariate differentiation to calculate all the	    *)
(*   multivariate derivatives of Fob speeds up the code considerably; in    *)
(*   particular, variations on the obvious Outer[D,eqns,vars] are *much*    *)
(*   slower.								    *)

AIMDerivatives[eqns_,0]:=
Function[tAIM, Evaluate[AIMSubBFuncsIntoEqns[eqns] /.
		Thread[AIMStateVars[eqns]->tAIM*AIMGenericArgs[eqns]]]] ;

AIMDerivatives[eqns_,deg_Integer]:= AIMDerivatives[eqns,deg] =
AIMDerivatives[eqns,deg-1]' ;


AIMSubBFuncsIntoEqns[eqns_]:= AIMSubBFuncsIntoEqns[eqns] =
With[{deveqns = eqns /.x_Symbol[t+i_.]:>Symbol[SymbolName[x]<>"AIMSS"] +x[t+i] 
				/.Equal->Subtract /.AIMSS[eqns,AIMZeroTol]},
 AIMSubOutLeadVars[eqns,deveqns,AIMMaxLead[eqns]]
] ;


AIMSubOutLeadVars[origeqns_,eqnssofar_,0]:= 
eqnssofar /.x_Symbol[t]->Apply[bFunc[x],AIMStateVars[origeqns]] /.
						eps[x_][t+i_]->Sigma*eps[x][t+i]

AIMSubOutLeadVars[origeqns_,eqnssofar_,lead_Integer]:=
With[{fwdsv=AIMStateVars[origeqns] /.t->t+lead},
With[{reducedeqns=eqnssofar /.x_Symbol[t+lead]->Apply[bFunc[x],fwdsv]},
 AIMSubOutLeadVars[origeqns, reducedeqns, lead-1]
]] ;  (* note how the use of recursion makes handling multiple leads simple *)


(* Print out model diagnostics (obviously) *)

AIMModelDiagnostics[eqns_] := (
Print["\nModel Diagnostics:"] ;
Print["Number of equations: ",AIMNEqns[eqns]] ;
Print["Number of variables: ",AIMNVars[eqns]] ;
Print["Number of shocks:    ",Length[AIMShocks[eqns]]] ;
Print["Maximum lag:         ",AIMMaxLag[eqns]] ;
Print["Maximum lead:        ",AIMMaxLead[eqns]] ;
Print["Lagged variables: ",AIMLagVars[eqns]] ;
Print["Shocks: ",AIMShocks[eqns]] ;
Print[" together with Sigma, these yield ",AIMNArgs[eqns]," state variables\n"];
Print["List of all variables: ",AIMVarNames[eqns]] ;
Print["Treating steady state and final coeff values < ",
				N[AIMZeroTol], " as zero (AIMZeroTol)\n"] ;
) ;


(* Everything that follows is Linear AIM.  See Anderson and Moore (1985)    *)
(*   for a description of the algorithm.  There are ways to improve the     *)
(*   speed and numerical stability of this implementation for larger        *)
(*   models, at the cost of complicating the code.  In particular, one can: *)
(* 1. Implement "exact" ShiftRights to reduce the number of calls to the    *)
(*     singular value decomposition.                                        *)
(* 2. Use a Schur decomposition instead of eigenvectors to compute the      *)
(*     left invariant subspace corresponding to the large eigenvalues of    *)
(*     AR1Form[shiftedh].  (This requires a user-written routine, because   *)
(*     Mathematica's Schur implementation does not sort eigenvalues into    *)
(*     any order.)                                                          *)
(* These improvements are implemented in the more advanced version of the   *)
(*   code.                                                                  *)

AIMLinearSoln[hmat_,nleads_Integer] :=
Module[{hrows,hcols,shiftedh,qmat,ar1eigvals,stabconds,bmat,smat},
 {hrows,hcols} = Dimensions[hmat] ;
 {shiftedh,qmat} = NestWhile[AIMShiftRight, {hmat, {}},
 				MatrixRank[AIMLastBlock[#[[1]]]] <hrows &,
 				1, nleads*hrows] ;
 Print["Lead matrix is full rank; computing stability conditions"] ;
 ar1eigvals = Eigenvalues[AIMAR1Form[shiftedh]] ;
 stabconds = Chop[Eigenvectors[Transpose[AIMAR1Form[shiftedh]],
		Length[Cases[Abs[ar1eigvals],i_/;i>1+AIMZeroTol]]], AIMZeroTol];
 Print["Model has ",Length[stabconds], " unstable roots"] ;
 qmat = Join[qmat,stabconds] ;
 If[Length[qmat]<hrows*nleads, Print["Multiple Linear Solutions"]; $Failed,
 If[Length[qmat]>hrows*nleads || MatrixRank[AIMLastBlock[qmat]]<hrows*nleads,
 					 Print["No Linear Solutions"]; $Failed,
 bmat = LinearSolve[AIMLastBlock[qmat],-AIMFirstBlocks[qmat]] ;
 smat = hmat .BlockMatrix[{{IdentityMatrix[hcols-hrows*nleads]},
				{ZeroMatrix[hrows*nleads,hrows],bmat}}] ;
 Chop[{Take[bmat,hrows], Inverse[AIMLastBlock[smat]]}, AIMZeroTol]
]]] ;


AIMLinearSoln[hmat_,0] :=
With[{cofb=LinearSolve[AIMLastBlock[hmat],-AIMFirstBlocks[hmat]]},
 Chop[{cofb, Inverse[AIMLastBlock[hmat]]}, AIMZeroTol]
] ;


AIMShiftRight[{hmatold_,qmatsofar_}]:=
Module[{hrows=Length[hmatold],svdu,svdsig,hmatnew,zerorows,firstblocks},
 {svdu,svdsig} = Take[SingularValueDecomposition[AIMLastBlock[hmatold]],2] ;
 hmatnew = Transpose[svdu] .hmatold ;
 zerorows = Map[List,Range[MatrixRank[svdsig]+1,hrows]] ;
Print["Shifting ",Length[zerorows]," linear combinations of equations forward"];
 firstblocks = AIMFirstBlocks[hmatnew] ;
 Chop[{ReplacePart[hmatnew,BlockMatrix[{{ZeroMatrix[hrows,hrows],firstblocks}}],
 							 zerorows,zerorows],
   Join[qmatsofar,firstblocks[[Flatten[zerorows]]]]}, AIMZeroTol]
] ;


AIMAR1Form[shiftedh_]:= AIMAR1Form[shiftedh] =
With[{hrows=Length[shiftedh],hcols=Length[shiftedh[[1]]]},
 Chop[BlockMatrix[{{ZeroMatrix[hcols-2*hrows,hrows],
						IdentityMatrix[hcols-2*hrows]},
	{LinearSolve[AIMLastBlock[shiftedh],-AIMFirstBlocks[shiftedh]]}}],
	AIMZeroTol]
] ;


AIMLastBlock[matrix_]:=
With[{dims=Dimensions[matrix]},
 SubMatrix[matrix,{1,dims[[2]]-dims[[1]]+1},{dims[[1]],dims[[1]]}]
] ;


AIMFirstBlocks[matrix_]:=
With[{dims=Dimensions[matrix]},
 SubMatrix[matrix,{1,1},{dims[[1]],dims[[2]]-dims[[1]]}]
] ;
