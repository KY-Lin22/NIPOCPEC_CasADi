function KKT_Matrix = computeKKT_Matrix(self, Fun, Jac, Hessian)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;

% regular parameter
nu_G_Repmat = repmat(Option.RegularParam.nu_G, 1, nStages);
nu_J_Repmat = repmat(Option.RegularParam.nu_J, 1, nStages);
nu_H_Repmat = repmat(Option.RegularParam.nu_H, 1, nStages);

% init
KKT_Matrix = struct('J', [], 'BL', [], 'BU', []);

%% KKT diagnal matrix: J
J = FunObj.J(nu_G_Repmat, nu_J_Repmat, nu_H_Repmat,...
    Fun.PSIgSigma_diagVec, Fun.PSIgG_diagVec, Fun.PSIphiGamma_diagVec, Fun.PSIphiPHI_diagVec,...
    Jac.Gx, Jac.Gu, Jac.Gp, Jac.Gw, Jac.Cx, Jac.Cu, Jac.Cp, Jac.Cw,...
    Jac.Fx, Jac.Fu, Jac.Fp, Jac.Fw, Jac.PHIx, Jac.PHIu, Jac.PHIp, Jac.PHIw,...
    Hessian);
KKT_Matrix.J = full(J);

%% KKT off-diagnal matrix: BL and BU
BL = [zeros(Dim.sigma + Dim.eta, Dim.Y);...
    zeros(Dim.lambda, Dim.LAMBDA), eye(Dim.x), zeros(Dim.lambda, Dim.u + Dim.p + Dim.w);...
    zeros(Dim.gamma + Dim.Z, Dim.Y)];
KKT_Matrix.BL = BL;  
KKT_Matrix.BU = BL';

end

