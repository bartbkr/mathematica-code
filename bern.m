(* Recursive identification code for Bernanke VAR term premium model*)
(* Barton Baker*)
(* Created 11/30/10*)

Needs["MathOptimizer`Optimize`"];

delta1=Import["G:\\MacroProg\\Data\\berndelta1.txt","Table"]
delta0=Import["G:\\MacroProg\\Data\\berndelta0.txt","Table"]
phi=Import["G:\\MacroProg\\Data\\bernphi.txt","Table"]
sigma=Import["G:\\MacroProg\\Data\\bernsigma.txt","Table"]
constants=Import["G:\\MacroProg\\Data\\bernconstants.txt","Table"]

B[1]=-delta1

A[1]=-delta0

B[n_] := DotProduct[Transpose[(phi - sigma*lambda1)], B[n - 1]] - delta1

A[n_]:= A[n-1]+DotProduct[Transpose[B[n-1]],(mu-sigma*lambda0)]+0.5*DotProduct[DotProduct[DotProduct[Transpose[B[n-1]],sigma],Transpose[sigma]],B[n-1]]-delta0



(**Import Data zero-coupon bond data**)

zeroco=Import["G:\\MacroProg\\Data\\nocoup.txt","Table"]

sixmth=zeroco[[All,1]]
oneyr=zeroco[[All,2]]
twoyr=zeroco[[All,3]]
threeyr=zeroco[[All,4]]
fiveyr=zeroco[[All,5]]
sevenyr=zeroco[[All,6]]
tenyr=zeroco[[All,7]]

eqns=Sum[



