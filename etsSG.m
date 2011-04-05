(* ::Package:: *)

(* Equations of the model, in a comma-delimited list (standard Mathematica
  format).  "eps" denotes a stochastic shock.  Note that your shocks
  *must* have a label eps[label][t] on them.  The label can be a number,
  letter, or word.  Your model can also contain shocks dated t+1 or
  later, but *not* t-1 or earlier (to handle the latter situation, you
  must define an auxiliary variable x = eps[label][t] and then consider
  x[t-1], which will work fine.)
 Each equation may have an equals sign ("==") or not.  If there is no
  equals sign, the code sets that particular equation equal to 0.
*)

eqns={
  Y[t] == A[t] *K[t-1]^alpha,
  Log[A[t]] == rho *Log[A[t-1]] + eps[a][t],
  K[t] == (1-delta) *K[t-1] + Inv[t], (* Mathematica treats "I" as Sqrt[-1] *)
  Y[t] == C[t] + Inv[t],
  C[t]^-gamma == beta *(1+r[t+1]) *C[t+1]^-gamma,
  r[t] == alpha *A[t] *K[t-1]^(alpha-1) - delta,
  Welf[t] == C[t]^(1-gamma) /(1-gamma) + beta *Welf[t+1]
}

(* variables listed here are transformed to logs and then approximated in
  logs rather than in levels (analogous to log-linearizing a linear model):
*)
logvars = {A, C, K, Y}
logrules = Map[#[x_]->E^(#[x])&, logvars]

(* Note: the above doesn't change the name of the variable, even though you
  transformed it, so you have to remember that you did the log transformation.
  If you want the name of the variable to change, use the following:
logrules = Map[#[x_]->E^(Symbol["log"<>SymbolName[#]][x])&, logvars]
*)

(* parameter values *)
parametervals={
  alpha->0.3,
  beta->0.99,
  gamma->1.1,
  delta->0.1,
  rho->0.8
}

(* substitute log transformation rules and parameter values into equations: *)
sgmodel = eqns /.logrules //.parametervals


(* complete variable list (Mathematica places in alphabetical order):
{A, C, Inv, K, r, Welf, Y}
*)

(* find the steady state *)
AIMSS[sgmodel]

(* if Mathematica doesn't find the steady state starting from 0, you have
  to give it an initial guess:*)
ssguess={0,0,0,1,0.1,-100,0}
AIMSS[sgmodel,AIMSSGuess->ssguess]

(*find the first-order approximation to solution for each variable at
  date t:
AIMSeries[sgmodel,1]
 This produces a comma-separated list of equations as answers.  The left-hand
  side variables are in levels.  The right-hand side variables are dated
  t-1 and earlier (shocks dated t) and are in *deviations* from steady state.
 Also note that the Mathematica TableForm command:
AIMSeries[sgmodel,1] //TableForm
  (*will produce more readable output.*)
*)

 (*
(* Find the second-order approximation:*)
AIMSeries[sgmodel,2]
 (*Note that the variable "Sigma" in the output denotes the "scale factor"
  of the stochastic shocks eps[t+1].  In particular, "Sigma" (with a
  capital "S") is a reserved symbol that you should probably not use in
  your equations unless you really know what you are doing (lowercase
  "sigma" is not used by the code and thus you may use it freely).
*)*)


(* Find the third-order approximation:*)
AIMSeries[sgmodel,3] //TableForm




(* Alternative formulation of the model:
eqns={
  Y[t] == A[t] *K[t]^alpha,
  Log[A[t]] == rho *Log[A[t-1]] + eps[a][t],
  K[t] == (1-delta) *K[t-1] + Inv[t-1], (* Mathematica treats "I" as Sqrt[-1] *)
  Y[t] == C[t] + Inv[t],
  C[t]^-gamma == beta *(1+r[t+1]) *C[t+1]^-gamma,
  r[t] == alpha *A[t] *K[t]^(alpha-1) - delta,
  Welf[t] == C[t]^(1-gamma) /(1-gamma) + beta *Welf[t+1]
}
 Note that the timing on K is different, although economically the model
  is the same.  By convention, AIM expresses the solution to any model in
  terms of variables dated t-1 and earlier (and shocks dated t), so this
  formulation of the model reports the solution in terms of A[t-1],
  Inv[t-1], and K[t-1] (and eps[a][t] and Sigma).
 For "predetermined" variables in general, you can vary how AIM reports
  the solution by varying the timing convention on the predetermined
  variable.  (Both solutions are economically identical.)
*)
