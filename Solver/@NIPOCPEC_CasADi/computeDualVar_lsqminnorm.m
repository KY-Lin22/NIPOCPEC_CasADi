function [eta, lambda] = computeDualVar_lsqminnorm(self, Var, s)
OCPEC = self.OCPEC;
Option = self.Option;

Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
LAMBDA_threshold = Option.FRP.LAMBDA_threshold;

% Jacobian evaluation
Jac = self.JacobianEvaluation(Var, s, 'Regular', []);

%% compute dual variables by solving an overdetermined linear equation system in a backward recursion manner
eta = zeros(Dim.eta, nStages);
lambda = zeros(Dim.lambda, nStages);
% initialize lambdaNext_n
lambdaNext_n = zeros(Dim.lambda, 1);

for n = nStages: -1 : 1
    % load
    Lx_n = Jac.Lx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Lu_n = Jac.Lu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Lp_n = Jac.Lp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Lw_n = Jac.Lw(:, 1 + (n - 1) * Dim.w : n * Dim.w);

    Gx_n = Jac.Gx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Gu_n = Jac.Gu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Gp_n = Jac.Gp(:, 1 + (n - 1) * Dim.p : n * Dim.p);

    Cx_n = Jac.Cx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Cu_n = Jac.Cu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Cp_n = Jac.Cp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    Cw_n = Jac.Cw(:, 1 + (n - 1) * Dim.w : n * Dim.w);
    C_grad_n = [Cx_n, Cu_n, Cp_n, Cw_n];
    
    Fx_n = Jac.Fx(:, 1 + (n - 1) * Dim.x : n * Dim.x);
    Fu_n = Jac.Fu(:, 1 + (n - 1) * Dim.u : n * Dim.u);
    Fp_n = Jac.Fp(:, 1 + (n - 1) * Dim.p : n * Dim.p); 
    F_grad_n = [Fx_n, Fu_n, Fp_n, zeros(Dim.lambda, Dim.w)];

    PHIp_n = Jac.PHIp(:, 1 + (n - 1) * Dim.p : n * Dim.p);
    PHIw_n = Jac.PHIw(:, 1 + (n - 1) * Dim.w : n * Dim.w);     
    
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

