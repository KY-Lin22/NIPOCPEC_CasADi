function Fun = FunctionEvaluation(self, Var, s, z, mode, FRP)
%FunctionEvaluation
%   return struct Fun with fileds including the following values
% L G C F PHI
% PSIg 
% PSIgSigma_diagVec
% PSIgG_diagVec
% PSIphi
% PSIphiGamma_diagVec 
% PSIphiPHI_diagVec

OCPEC = self.OCPEC;
FunObj = self.FunObj;
s_Repmat = repmat(s, 1, OCPEC.nStages);
z_Repmat = repmat(z, 1, OCPEC.nStages);
xPrev = [OCPEC.InitState, Var.x(:, 1 : end - 1)];

% cost function
switch mode
    case 'Regular'
        L_T = FunObj.L_T(Var.x(:, end), Var.u(:, end), Var.p(:, end), Var.w(:, end));
        L_S = FunObj.L_S(Var.x, Var.u, Var.p, Var.w);        
        L = full(L_S);
        L(:, end) = L(:, end) + full(L_T);
    case 'FRP'
        L = FunObj.FRP_L(Var.x, Var.u, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end
% constraint
G = FunObj.G(Var.x, Var.u, Var.p, Var.w);
C = FunObj.C(Var.x, Var.u, Var.p, Var.w);
F = FunObj.F(xPrev, Var.x, Var.u, Var.p, Var.w);
PHI = FunObj.PHI(Var.x, Var.u, Var.p, Var.w, s_Repmat);
% FB function and Jacobian for G
PSIg = FunObj.FB_G(Var.sigma, full(G), z_Repmat);
[PSIgSigma_diagVec, PSIgG_diagVec] = FunObj.FB_G_grad(Var.sigma, full(G), z_Repmat);
% FB function and Jacobian for PHI
PSIphi = FunObj.FB_PHI(Var.gamma, full(PHI), z_Repmat);
[PSIphiGamma_diagVec, PSIphiPHI_diagVec] = FunObj.FB_PHI_grad(Var.gamma, full(PHI), z_Repmat);
%
Fun = struct('L', full(L),...
    'G', full(G), 'C', full(C), 'F', full(F), 'PHI', full(PHI),...
    'PSIg', full(PSIg),...
    'PSIgSigma_diagVec', full(PSIgSigma_diagVec),...
    'PSIgG_diagVec', full(PSIgG_diagVec),...
    'PSIphi', full(PSIphi),...
    'PSIphiGamma_diagVec', full(PSIphiGamma_diagVec),...
    'PSIphiPHI_diagVec', full(PSIphiPHI_diagVec));

end

