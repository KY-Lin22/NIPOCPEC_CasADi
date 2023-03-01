function Var_SOC = SecondOrderCorrection(self, Var, Jac, KKT_Residual, KKT_Matrix, Fun_full)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
nu_G_Repmat = repmat(Option.RegularParam.nu_G, 1, nStages);

% SOC KKT residual
KKT_Residual_SOC = struct('G_Fsb', [], 'C_Fsb', [], 'F_Fsb', [], 'PHI_Fsb', [],...
   'HxTlambdaNext', KKT_Residual.HxTlambdaNext,  'HuT', KKT_Residual.HuT, 'HpT', KKT_Residual.HpT, 'HwT', KKT_Residual.HwT);

G_Fsb_Correct = FunObj.G_Fsb(Fun_full.PSIg, Jac.PSIgG_diagVec, nu_G_Repmat);
PHI_Fsb_Correct = FunObj.PHI_Fsb(Fun_full.PSIphi, Jac.PSIphiPHI_diagVec, nu_G_Repmat);

KKT_Residual_SOC.G_Fsb = KKT_Residual.G_Fsb + full(G_Fsb_Correct);
KKT_Residual_SOC.C_Fsb = KKT_Residual.C_Fsb + Fun_full.C;
KKT_Residual_SOC.F_Fsb = KKT_Residual.F_Fsb + Fun_full.F;
KKT_Residual_SOC.PHI_Fsb = KKT_Residual.PHI_Fsb + full(PHI_Fsb_Correct);

% compute second order correction direction and var
[dY_SOC, ~] = self.SearchDirection_Riccati(KKT_Residual_SOC, KKT_Matrix);

Var_SOC.sigma  = Var.sigma  + dY_SOC(              1 : Dim.Node(1), :);
Var_SOC.eta    = Var.eta    + dY_SOC(Dim.Node(1) + 1 : Dim.Node(2), :);
Var_SOC.lambda = Var.lambda + dY_SOC(Dim.Node(2) + 1 : Dim.Node(3), :);
Var_SOC.gamma  = Var.gamma  + dY_SOC(Dim.Node(3) + 1 : Dim.Node(4), :);
Var_SOC.x      = Var.x      + dY_SOC(Dim.Node(4) + 1 : Dim.Node(5), :);
Var_SOC.u      = Var.u      + dY_SOC(Dim.Node(5) + 1 : Dim.Node(6), :);
Var_SOC.p      = Var.p      + dY_SOC(Dim.Node(6) + 1 : Dim.Node(7), :);
Var_SOC.w      = Var.w      + dY_SOC(Dim.Node(7) + 1 : Dim.Node(8), :);
end

