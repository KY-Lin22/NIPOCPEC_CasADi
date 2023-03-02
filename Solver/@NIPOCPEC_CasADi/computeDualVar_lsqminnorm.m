function [eta, lambda] = computeDualVar_lsqminnorm(self, Var, s)
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;

Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
LAMBDA_threshold = Option.FRP.LAMBDA_threshold;

% Jacobian evaluation
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

[Gx, Gu, Gp, ~] = FunObj.G_grad(Var.x, Var.u, Var.p, Var.w);
[Cx, Cu, Cp, Cw] = FunObj.C_grad(Var.x, Var.u, Var.p, Var.w);
[Fx, Fu, Fp, Fw] = FunObj.F_grad(Var.x, Var.u, Var.p, Var.w);
[~, ~, PHIp, PHIw] = FunObj.PHI_grad(Var.x, Var.u, Var.p, Var.w, s);
Gx = full(Gx);
Gu = full(Gu);
Gp = full(Gp);
Cx = full(Cx);
Cu = full(Cu);
Cp = full(Cp);
Cw = full(Cw);
Fx = full(Fx);
Fu = full(Fu);
Fp = full(Fp);
Fw = full(Fw);
PHIp = full(PHIp);
PHIw = full(PHIw);



%% compute dual variables by solving an overdetermined linear equation system in a backward recursion manner
eta = zeros(Dim.eta, nStages);
lambda = zeros(Dim.lambda, nStages);
% initialize lambdaNext_n
lambdaNext_n = zeros(Dim.lambda, 1);

for n = nStages: -1 : 1
    % load
    Lx_n = Lx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Lu_n = Lu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Lp_n = Lp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Lw_n = Lw(:, 1 + (n - 1) * Dim.w : n * Dim.w);

    Gx_n = Gx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Gu_n = Gu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Gp_n = Gp(:, 1 + (n - 1) * Dim.p : n * Dim.p);

    Cx_n = Cx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Cu_n = Cu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Cp_n = Cp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Cw_n = Cw(:, 1 + (n - 1) * Dim.w : n * Dim.w);
    C_grad_n = [Cx_n, Cu_n, Cp_n, Cw_n];
    
    Fx_n = Fx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Fu_n = Fu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Fp_n = Fp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Fw_n = Fw(:, 1 + (n - 1) * Dim.w : n * Dim.w);   
    F_grad_n = [Fx_n, Fu_n, Fp_n, Fw_n];

    PHIp_n = PHIp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    PHIw_n = PHIw(:, 1 + (n - 1) * Dim.w : n * Dim.w);     
    
    % solve overdetermined linear equation system
    sigma_n = Var.sigma(:, n);
    gamma_n = Var.gamma(:, n);
    
    A_n = [C_grad_n', F_grad_n'];
    b_n = [Lx_n' - Gx_n' * sigma_n                      + lambdaNext_n;...
           Lu_n' - Gu_n' * sigma_n;...
           Lp_n' - Gp_n' * sigma_n  - PHIp_n' * gamma_n;...
           Lw_n'                    - PHIw_n' * gamma_n];
    EtaLambda_n = lsqminnorm(A_n,-b_n);
    if norm(EtaLambda_n, Inf) > LAMBDA_threshold
        EtaLambda_n = zeros(Dim.eta + Dim.lambda, 1);
    end
    eta_n    = EtaLambda_n(1 : Dim.eta, :);
    lambda_n = EtaLambda_n(Dim.eta + 1 : end, :);
    
    % record and prepare for next iteration
    eta(:, n) = eta_n;
    lambda(:, n) = lambda_n;    
    lambdaNext_n = lambda_n;    
end

end

