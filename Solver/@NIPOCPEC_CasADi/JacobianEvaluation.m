function Jac = JacobianEvaluation(self, Var, s, mode, FRP)
%JacobianEvaluation
%   return struct Jac with fileds including the following values
% Lx Lu Lp Lw(regular recorded as sparse, FRP record as DM type)
% Gx Gu Gp(record as DM type)
% Cx Cu Cp Cw(record as DM type)
% Fx Fu Fp(record as DM type)
% PHIp PHIw(record as DM type)

OCPEC = self.OCPEC;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;

% cost function Jacobian
switch mode
    case 'Regular'
        [Lx_T, Lp_T] = FunObj.L_T_grad(Var.x(:, end), Var.p(:, end));
        [Lx_S,Lu_S, Lp_S, Lw_S] = FunObj.L_S_grad(Var.x, Var.u, Var.p, Var.w);
        Lx = sparse(Lx_S);
        Lu = sparse(Lu_S);
        Lp = sparse(Lp_S);
        Lw = sparse(Lw_S);
        Lx(:, 1 + (nStages - 1) * Dim.x : end) = Lx(:, 1 + (nStages - 1) * Dim.x : end) + sparse(Lx_T);
        Lp(:, 1 + (nStages - 1) * Dim.p : end) = Lp(:, 1 + (nStages - 1) * Dim.p : end) + sparse(Lp_T);
    case 'FRP'
        [Lx, Lu, Lp, Lw] = FunObj.FRP_L_grad(Var.x, Var.u, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
end

% constraint Jacobian
[Gx, Gu, Gp] = FunObj.G_grad(Var.x, Var.u, Var.p);
[Cx, Cu, Cp, Cw] = FunObj.C_grad(Var.x, Var.u, Var.p, Var.w);
[Fx, Fu, Fp] = FunObj.F_grad(Var.x, Var.u, Var.p);
[PHIp, PHIw] = FunObj.PHI_grad(Var.p, Var.w, s);

%
Jac = struct('Lx', Lx, 'Lu', Lu, 'Lp', Lp, 'Lw', Lw,...
    'Gx', Gx, 'Gu', Gu, 'Gp', Gp,...
    'Cx', Cx, 'Cu', Cu, 'Cp', Cp, 'Cw', Cw,...
    'Fx', Fx, 'Fu', Fu, 'Fp', Fp,...
    'PHIp', PHIp, 'PHIw', PHIw);

end

