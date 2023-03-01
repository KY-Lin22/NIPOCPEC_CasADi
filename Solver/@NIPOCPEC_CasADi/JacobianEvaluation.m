function Jac = JacobianEvaluation(self, Var, Fun, s, z, mode, FRP)
%JacobianEvaluation
%   return struct Jac with fileds including the following values
% Lx Lu Lp Lw
% Gx Gu Gp Gw
% Cx Cu Cp Cw
% Fx Fu Fp Fw
% PHIx PHIu PHIp PHIw
% PSIgSigma_diagVec
% PSIgG_diagVec
% PSIphiGamma_diagVec 
% PSIphiPHI_diagVec

OCPEC = self.OCPEC;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
s_Repmat = repmat(s, 1, nStages);
z_Repmat = repmat(z, 1, OCPEC.nStages);
% cost function Jacobian
switch mode
    case 'Regular'
        [Lx_T,Lu_T, Lp_T, Lw_T] = FunObj.L_T_grad(Var.x(:, end), Var.u(:, end), Var.p(:, end), Var.w(:, end));
        [Lx_S,Lu_S, Lp_S, Lw_S] = FunObj.L_S_grad(Var.x, Var.u, Var.p, Var.w);
        Lx = full(Lx_S);
        Lu = full(Lu_S);
        Lp = full(Lp_S);
        Lw = full(Lw_S);
        Lx(:, 1 + (nStages - 1) * Dim.x : end) = Lx(:, 1 + (nStages - 1) * Dim.x : end) + full(Lx_T);
        Lu(:, 1 + (nStages - 1) * Dim.u : end) = Lu(:, 1 + (nStages - 1) * Dim.u : end) + full(Lu_T);
        Lp(:, 1 + (nStages - 1) * Dim.p : end) = Lp(:, 1 + (nStages - 1) * Dim.p : end) + full(Lp_T);
        Lw(:, 1 + (nStages - 1) * Dim.w : end) = Lw(:, 1 + (nStages - 1) * Dim.w : end) + full(Lw_T);
    case 'FRP'
        [Lx, Lu, Lp, Lw] = FunObj.FRP_L_grad(Var.x, Var.u, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end

% constraint Jacobian
[Gx, Gu, Gp, Gw] = FunObj.G_grad(Var.x, Var.u, Var.p, Var.w);
[Cx, Cu, Cp, Cw] = FunObj.C_grad(Var.x, Var.u, Var.p, Var.w);
[Fx, Fu, Fp, Fw] = FunObj.F_grad(Var.x, Var.u, Var.p, Var.w);
[PHIx, PHIu, PHIp, PHIw] = FunObj.PHI_grad(Var.x, Var.u, Var.p, Var.w, s_Repmat);

% FB Jacobian for G and PHI
[PSIgSigma_diagVec, PSIgG_diagVec] = FunObj.FB_G_grad(Var.sigma, Fun.G, z_Repmat);
[PSIphiGamma_diagVec, PSIphiPHI_diagVec] = FunObj.FB_PHI_grad(Var.gamma, Fun.PHI, z_Repmat);
%
Jac = struct('Lx', full(Lx), 'Lu', full(Lu), 'Lp', full(Lp), 'Lw', full(Lw),...
    'Gx', full(Gx), 'Gu', full(Gu), 'Gp', full(Gp), 'Gw', full(Gw),...
    'Cx', full(Cx), 'Cu', full(Cu), 'Cp', full(Cp), 'Cw', full(Cw),...
    'Fx', full(Fx), 'Fu', full(Fu), 'Fp', full(Fp), 'Fw', full(Fw),...
    'PHIx', full(PHIx), 'PHIu', full(PHIu), 'PHIp', full(PHIp), 'PHIw', full(PHIw),...
    'PSIgSigma_diagVec', full(PSIgSigma_diagVec),...
    'PSIgG_diagVec', full(PSIgG_diagVec),...
    'PSIphiGamma_diagVec', full(PSIphiGamma_diagVec),...
    'PSIphiPHI_diagVec', full(PSIphiPHI_diagVec));

end

