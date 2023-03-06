function [eta, lambda] = computeDualVar_lsqminnorm(self, Var, s)
OCPEC = self.OCPEC;
FunObj = self.FunObj;
Option = self.Option;

Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
LAMBDA_threshold = Option.FRP.LAMBDA_threshold;

% Jacobian evaluation
Jac = self.JacobianEvaluation(Var, s, 'Regular', []);

%% compute dual variables by solving an overdetermined linear equation system in a backward recursion manner
% compute A_lsqminnorm and b_lsqminnorm
[A_lsqminnorm, b_lsqminnorm] =...
    FunObj.A_b_lsqminnorm(Var.sigma, Var.gamma,...
    Jac.Lx, Jac.Lu, Jac.Lp, Jac.Lw,...
    Jac.Gx, Jac.Gu, Jac.Gp,...
    Jac.Cx, Jac.Cu, Jac.Cp, Jac.Cw,...
    Jac.Fx, Jac.Fu, Jac.Fp,...
    Jac.PHIp, Jac.PHIw);
A_lsqminnorm = full(A_lsqminnorm);
b_lsqminnorm = full(b_lsqminnorm);
% initialize lambdaNext_n
lambdaNext_n = zeros(Dim.lambda, 1);
% initialize eta and lambda
eta = zeros(Dim.eta, nStages);
lambda = zeros(Dim.lambda, nStages);

for n = nStages: -1 : 1  
    % solve overdetermined linear equation system
    A_n = A_lsqminnorm(:, (Dim.eta + Dim.lambda) * (n - 1) + 1 : (Dim.eta + Dim.lambda) * n);
    b_n = b_lsqminnorm(:, n);
    b_n(1 : Dim.x, :) = b_n(1 : Dim.x, :) + lambdaNext_n;
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

