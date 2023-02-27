clear all
clc

addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

Dim.sigma = 10;
Dim.eta = 10;
Dim.lambda = 15;
Dim.gamma = 10;
Dim.x = 10;
nStages = 100;

% SX fun
sigma_SX = SX.sym('sigma_SX', Dim.sigma, 1);
eta_SX = SX.sym('eta_SX', Dim.eta, 1);
lambda_SX = SX.sym('lambda_SX', Dim.lambda, 1);
gamma_SX = SX.sym('gamma_SX', Dim.gamma, 1);

Lx_SX = SX.sym('Lx_SX', 1, Dim.x);
Gx_SX = SX.sym('Gx_SX', Dim.sigma, Dim.x);
Cx_SX = SX.sym('Cx_SX', Dim.eta, Dim.x);
Fx_SX = SX.sym('Fx_SX', Dim.lambda, Dim.x);
PHIx_SX = SX.sym('PHIx_SX', Dim.gamma, Dim.x);

LAGx_SX = Lx_SX - sigma_SX'*Gx_SX + eta_SX'*Cx_SX + lambda_SX'*Fx_SX - gamma_SX'*PHIx_SX;
LAGx_SX_Fun = Function('LAGx_SX', {sigma_SX, eta_SX, lambda_SX, gamma_SX, Lx_SX, Gx_SX, Cx_SX, Fx_SX, PHIx_SX}, {LAGx_SX},...
    {'sigma', 'eta', 'lambda', 'gamma', 'Lx', 'Gx', 'Cx', 'Fx', 'PHIx'}, {'LAGx'});

% SX fun Map
LAGx_SX_Fun_Map = LAGx_SX_Fun.map(nStages);

% MX fun
sigma_MX = MX.sym('sigma_MX', Dim.sigma, 1);
eta_MX = MX.sym('eta_MX', Dim.eta, 1);
lambda_MX = MX.sym('lambda_MX', Dim.lambda, 1);
gamma_MX = MX.sym('gamma_MX', Dim.gamma, 1);

Lx_MX = MX.sym('Lx_MX', 1, Dim.x);
Gx_MX = MX.sym('Gx_MX', Dim.sigma, Dim.x);
Cx_MX = MX.sym('Cx_MX', Dim.eta, Dim.x);
Fx_MX = MX.sym('Fx_MX', Dim.lambda, Dim.x);
PHIx_MX = MX.sym('PHIx_MX', Dim.gamma, Dim.x);

LAGx_MX = Lx_MX - sigma_MX'*Gx_MX + eta_MX'*Cx_MX + lambda_MX'*Fx_MX - gamma_MX'*PHIx_MX;
LAGx_MX_Fun = Function('LAGx_MX', {sigma_MX, eta_MX, lambda_MX, gamma_MX, Lx_MX, Gx_MX, Cx_MX, Fx_MX, PHIx_MX}, {LAGx_MX},...
    {'sigma', 'eta', 'lambda', 'gamma', 'Lx', 'Gx', 'Cx', 'Fx', 'PHIx'}, {'LAGx'});

% MX fun map
LAGx_MX_Fun_Map = LAGx_MX_Fun.map(nStages);
%% test
sigma = randn(Dim.sigma, nStages);
eta = randn(Dim.eta, nStages);
lambda = randn(Dim.lambda, nStages);
gamma = randn(Dim.gamma, nStages);
Lx = randn(1, Dim.x * nStages);
Gx = randn(Dim.sigma, Dim.x * nStages);
Cx = randn(Dim.eta, Dim.x * nStages);
Fx = randn(Dim.lambda, Dim.x * nStages);
PHIx = randn(Dim.gamma, Dim.x * nStages);

% SX fun for
LAGx_SX_Value = zeros(1, Dim.x * nStages);
tic
for n = 1 : nStages
    LAGx_SX_Value_n = LAGx_SX_Fun(sigma(:, n), eta(:, n), lambda(:, n), gamma(:, n),...
        Lx(:, (n - 1) * Dim.x + 1 : n * Dim.x), Gx(:, (n - 1) * Dim.x + 1 : n * Dim.x), Cx(:, (n - 1) * Dim.x + 1 : n * Dim.x),...
        Fx(:, (n - 1) * Dim.x + 1 : n * Dim.x), PHIx(:, (n - 1) * Dim.x + 1 : n * Dim.x));
    LAGx_SX_Value(:, (n - 1) * Dim.x + 1 : n * Dim.x) = full(LAGx_SX_Value_n);
end
toc

% SX fun map
tic
LAGx_SX_Map_Value = LAGx_SX_Fun_Map(sigma, eta, lambda, gamma, Lx, Gx, Cx, Fx, PHIx);
LAGx_SX_Map_Value = full(LAGx_SX_Map_Value);
toc

% MX fun for
LAGx_MX_Value = zeros(1, Dim.x * nStages);
tic
for n = 1 : nStages
    LAGx_MX_Value_n = LAGx_MX_Fun(sigma(:, n), eta(:, n), lambda(:, n), gamma(:, n),...
        Lx(:, (n - 1) * Dim.x + 1 : n * Dim.x), Gx(:, (n - 1) * Dim.x + 1 : n * Dim.x), Cx(:, (n - 1) * Dim.x + 1 : n * Dim.x),...
        Fx(:, (n - 1) * Dim.x + 1 : n * Dim.x), PHIx(:, (n - 1) * Dim.x + 1 : n * Dim.x));
    LAGx_MX_Value(:, (n - 1) * Dim.x + 1 : n * Dim.x) = full(LAGx_MX_Value_n);
end
toc

% MX fun map
tic
LAGx_MX_Map_Value = LAGx_MX_Fun_Map(sigma, eta, lambda, gamma, Lx, Gx, Cx, Fx, PHIx);
LAGx_MX_Map_Value = full(LAGx_MX_Map_Value);
toc

% direct computation
LAGx_MX_Value = zeros(1, Dim.x * nStages);
tic
for n = 1 : nStages
    LAGx_MX_Value_n = Lx(:, (n - 1) * Dim.x + 1 : n * Dim.x)...
        - sigma(:, n)'  * Gx(:, (n - 1) * Dim.x + 1 : n * Dim.x)...
        + eta(:, n)'    * Cx(:, (n - 1) * Dim.x + 1 : n * Dim.x)...
        + lambda(:, n)' * Fx(:, (n - 1) * Dim.x + 1 : n * Dim.x)...
        - gamma(:, n)'  * PHIx(:, (n - 1) * Dim.x + 1 : n * Dim.x);    
    LAGx_MX_Value(:, (n - 1) * Dim.x + 1 : n * Dim.x) = full(LAGx_MX_Value_n);
end
toc