function KKT_Matrix = computeKKT_Matrix(self, Jac, Hessian)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;
Dim = OCPEC.Dim;

% regular parameter
nu_G = Option.RegularParam.nu_G;
nu_J = Option.RegularParam.nu_J;
nu_H = Option.RegularParam.nu_H;

% init
KKT_Matrix = struct('J', [], 'BL', [], 'BU', []);

%% KKT diagnal matrix: J
J = FunObj.J(nu_G, nu_J, nu_H,...
    Jac.PSIgSigma_diagVec, Jac.PSIgG_diagVec, Jac.PSIphiGamma_diagVec, Jac.PSIphiPHI_diagVec,...
    Jac.Gx, Jac.Gu, Jac.Gp,...
    Jac.Cx, Jac.Cu, Jac.Cp, Jac.Cw,...
    Jac.Fx, Jac.Fu, Jac.Fp,...
    Jac.PHIp, Jac.PHIw,...
    Hessian);
KKT_Matrix.J = full(J);

%% KKT off-diagnal matrix: BL and BU
BL = [zeros(Dim.sigma + Dim.eta, Dim.Y);...
    zeros(Dim.lambda, Dim.LAMBDA), eye(Dim.x), zeros(Dim.lambda, Dim.u + Dim.p + Dim.w);...
    zeros(Dim.gamma + Dim.Z, Dim.Y)];
KKT_Matrix.BL = BL;  
KKT_Matrix.BU = BL';

end

