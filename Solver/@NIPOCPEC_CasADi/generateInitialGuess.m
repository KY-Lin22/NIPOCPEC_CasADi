function generateInitialGuess(self)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;

Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
% init
Var = struct('x', [], 'u', [], 'p', [], 'w', [],...
    'sigma', [], 'eta', [], 'lambda', [], 'gamma', []);

disp('Generating Initial Guess...')
%% generate initial guess for primal variable
% x
Var.x = zeros(Dim.x, nStages);
% u
Var.u = randn(Dim.u, nStages);
% p
p = zeros(Dim.p, nStages);
for i = 1 : Dim.p
    if (OCPEC.lbp(i) == 0) && (OCPEC.ubp(i) == Inf)
        % nonlinear complementary problem
        p(i, :) = ones(1, nStages); % p > 0
    else
        % box constraint variation inequality
        p(i, :) = repmat(1/2*(OCPEC.lbp(i) + OCPEC.ubp(i)), 1, nStages); % l < = p < = u
    end
end
Var.p = p;
% w
K = self.FunObj.K(Var.x, Var.u, Var.p);
Var.w = full(K);

%% generate initial guess for dual variable
Var.sigma  = ones(Dim.sigma, nStages); % sigma >= 0
Var.eta = randn(Dim.eta, nStages);
Var.lambda = randn(Dim.lambda, nStages);
Var.gamma = ones(Dim.gamma, nStages); % gamma >= 0

%% save initial guess
save('Gen_InitialGuess.mat', 'Var');
disp('Done!')
end

