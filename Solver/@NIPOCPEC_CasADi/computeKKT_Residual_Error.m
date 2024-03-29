function [KKT_Residual, KKT_Error] = computeKKT_Residual_Error(self, Var, Fun, Jac)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
nu_G = Option.RegularParam.nu_G;
lambdaNext = [Var.lambda(:, 2 : end), zeros(Dim.lambda, 1)];

%% KKT Residual
% primal feasibility
G_Fsb = FunObj.G_Fsb(Fun.PSIg, Jac.PSIgG_diagVec, nu_G);
C_Fsb = Fun.C;
F_Fsb = Fun.F;
PHI_Fsb = FunObj.PHI_Fsb(Fun.PSIphi, Jac.PSIphiPHI_diagVec, nu_G);
% dual feasibility
[HxT, HuT, HpT, HwT] = FunObj.HAM_grad(Var.sigma, Var.eta, Var.lambda, Var.gamma,...
    Jac.Lx, Jac.Lu, Jac.Lp, Jac.Lw,...
    Jac.Gx, Jac.Gu, Jac.Gp, ...
    Jac.Cx, Jac.Cu, Jac.Cp, Jac.Cw,...
    Jac.Fx, Jac.Fu, Jac.Fp, ...
    Jac.PHIp, Jac.PHIw);
%
KKT_Residual = struct('G_Fsb', full(G_Fsb), 'C_Fsb', C_Fsb, 'F_Fsb', F_Fsb, 'PHI_Fsb', full(PHI_Fsb),...
    'HxTlambdaNext', full(HxT) + lambdaNext, 'HuT', full(HuT), 'HpT', full(HpT), 'HwT', full(HwT));

%% KKT Error
% compute scaling parameter for KKT stationarity
scaling_max = 100;
LAMBDA_norm = norm(reshape(Var.sigma, [], 1), 1) + norm(reshape(Var.eta, [], 1), 1) ...
             + norm(reshape(Var.lambda, [], 1), 1) + norm(reshape(Var.gamma, [], 1), 1);          
LAMBDA_trial = (LAMBDA_norm)/(nStages * Dim.LAMBDA);
stationarityScaling = max([scaling_max, LAMBDA_trial])/scaling_max;

% compute KKT_Error
KKT_Error.Feasibility = max([norm(reshape(Fun.PSIg, [], 1), Inf),...
    norm(reshape(Fun.C, [], 1), Inf),...
    norm(reshape(Fun.F, [], 1), Inf),...
    norm(reshape(Fun.PSIphi, [], 1), Inf)]);

KKT_Error.Stationarity = max([norm(reshape(KKT_Residual.HxTlambdaNext, [], 1), Inf),...
    norm(reshape(KKT_Residual.HuT, [], 1), Inf),...
    norm(reshape(KKT_Residual.HpT, [], 1), Inf),...
    norm(reshape(KKT_Residual.HwT, [], 1), Inf)])/stationarityScaling;

KKT_Error.Total = max([KKT_Error.Feasibility, KKT_Error.Stationarity]);
end

