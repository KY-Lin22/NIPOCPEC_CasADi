function FunObj = createFunObj(self)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
OCPEC = self.OCPEC;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;

%% function evaluation
% f
f_FunObj = Function('f', {OCPEC.x, OCPEC.u, OCPEC.p}, {OCPEC.f},...
    {'x', 'u', 'p'}, {'f'});
FunObj.f = f_FunObj.map(nStages);
% K
K_FunObj = Function('K', {OCPEC.x, OCPEC.u, OCPEC.p}, {OCPEC.K}, ...
    {'x', 'u', 'p'}, {'K'});
FunObj.K = K_FunObj.map(nStages);

% L_T
FunObj.L_T = Function('L_T', {OCPEC.x, OCPEC.p}, {OCPEC.L_T},...
    {'x', 'p'}, {'L_T'});
% L_S
L_S_FunObj = Function('L_S', {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w}, {OCPEC.L_S},...
    {'x', 'u', 'p', 'w'}, {'L_S'});
FunObj.L_S = L_S_FunObj.map(nStages);
% G
G_FunObj = Function('G',{OCPEC.x, OCPEC.u, OCPEC.p}, {OCPEC.G},...
    {'x', 'u', 'p'}, {'G'});
FunObj.G = G_FunObj.map(nStages);
% C
C_FunObj = Function('C', {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w}, {OCPEC.C},...
    {'x', 'u', 'p', 'w'}, {'C'});
FunObj.C = C_FunObj.map(nStages);
% F
F_FunObj = Function('F', {OCPEC.xPrev, OCPEC.x, OCPEC.u, OCPEC.p}, {OCPEC.F},...
    {'xPrev', 'x', 'u', 'p'}, {'F'});
FunObj.F = F_FunObj.map(nStages);
% PHI
PHI_FunObj = Function('PHI', {OCPEC.p, OCPEC.w, OCPEC.s}, {OCPEC.PHI},...
    {'p', 'w', 's'}, {'PHI'});
FunObj.PHI = PHI_FunObj.map(nStages);

%% Jacobian evaluation
% f_grad
fx_formula = jacobian(OCPEC.f, OCPEC.x);
fu_formula = jacobian(OCPEC.f, OCPEC.u);
fp_formula = jacobian(OCPEC.f, OCPEC.p);
f_grad_FunObj = Function('f_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p},...
    {fx_formula, fu_formula, fp_formula},...
    {'x', 'u', 'p'},...
    {'fx', 'fu', 'fp'});
FunObj.f_grad = f_grad_FunObj.map(nStages);

% L_T_grad
Lx_T_formula = jacobian(OCPEC.L_T, OCPEC.x);
Lp_T_formula = jacobian(OCPEC.L_T, OCPEC.p);
FunObj.L_T_grad = Function('L_T_grad',...
    {OCPEC.x, OCPEC.p},...
    {Lx_T_formula, Lp_T_formula},...
    {'x', 'p'},...
    {'Lx_T','Lp_T'});
% L_S_grad
Lx_S_formula = jacobian(OCPEC.L_S, OCPEC.x);
Lu_S_formula = jacobian(OCPEC.L_S, OCPEC.u);
Lp_S_formula = jacobian(OCPEC.L_S, OCPEC.p);
Lw_S_formula = jacobian(OCPEC.L_S, OCPEC.w);
L_S_grad_FunObj = Function('L_S_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w},...
    {Lx_S_formula, Lu_S_formula, Lp_S_formula, Lw_S_formula},...
    {'x', 'u', 'p', 'w'},...
    {'Lx_S', 'Lu_S', 'Lp_S', 'Lw_S'});
FunObj.L_S_grad = L_S_grad_FunObj.map(nStages);
% G_grad
Gx_formula = jacobian(OCPEC.G, OCPEC.x);
Gu_formula = jacobian(OCPEC.G, OCPEC.u);
Gp_formula = jacobian(OCPEC.G, OCPEC.p);
G_grad_FunObj = Function('G_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p},...
    {Gx_formula, Gu_formula, Gp_formula},...
    {'x', 'u', 'p'},...
    {'Gx', 'Gu', 'Gp'});
FunObj.G_grad = G_grad_FunObj.map(nStages);
% C_grad
Cx_formula = jacobian(OCPEC.C, OCPEC.x);
Cu_formula = jacobian(OCPEC.C, OCPEC.u);
Cp_formula = jacobian(OCPEC.C, OCPEC.p);
Cw_formula = jacobian(OCPEC.C, OCPEC.w);
C_grad_FunObj = Function('C_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w},...
    {Cx_formula, Cu_formula, Cp_formula, Cw_formula},...
    {'x', 'u', 'p', 'w'},...
    {'Cx', 'Cu', 'Cp', 'Cw'});
FunObj.C_grad = C_grad_FunObj.map(nStages);
% F_grad
Fx_formula = fx_formula * OCPEC.timeStep - SX.eye(Dim.x);
Fu_formula = fu_formula * OCPEC.timeStep;
Fp_formula = fp_formula * OCPEC.timeStep;
F_grad_FunObj = Function('F_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p},...
    {Fx_formula, Fu_formula, Fp_formula},...
    {'x', 'u', 'p'},...
    {'Fx', 'Fu', 'Fp'});
FunObj.F_grad = F_grad_FunObj.map(nStages);
% PHI_grad
PHIp_formula = jacobian(OCPEC.PHI, OCPEC.p);
PHIw_formula = jacobian(OCPEC.PHI, OCPEC.w);
PHI_grad_FunObj = Function('PHI_grad',...
    {OCPEC.p, OCPEC.w, OCPEC.s},...
    {PHIp_formula, PHIw_formula},...
    {'p', 'w', 's'},...
    {'PHIp', 'PHIw'});
FunObj.PHI_grad = PHI_grad_FunObj.map(nStages);

%% Hessian evaluation
% terminal cost function hessian
[L_T_hessian_formula, ~] = hessian(OCPEC.L_T, [OCPEC.x; OCPEC.u; OCPEC.p; OCPEC.w]);
FunObj.L_T_hessian = Function('L_T_hessian',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w},...
    {L_T_hessian_formula},...
    {'x', 'u', 'p', 'w'},...
    {'L_T_hessian'});
% stage cost function hessian
[L_S_hessian_formula, ~] = hessian(OCPEC.L_S, [OCPEC.x; OCPEC.u; OCPEC.p; OCPEC.w]);
L_S_hessian_FunObj = Function('L_S_hessian',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w},...
    {L_S_hessian_formula},...
    {'x', 'u', 'p', 'w'},...
    {'L_S_hessian'});
FunObj.L_S_hessian = L_S_hessian_FunObj.map(nStages);

% Hamiltonian stage Hessian
sigma  = SX.sym('sigma', Dim.sigma, 1); % dual variable for G
eta    = SX.sym('eta', Dim.eta, 1); % dual variable for C
lambda = SX.sym('lambda', Dim.lambda, 1); % dual variable for F, 
gamma  = SX.sym('gamma', Dim.gamma, 1); % dual variable for PHI

HAM_formula = OCPEC.L_S - sigma'*OCPEC.G + eta'*OCPEC.C + lambda'*(OCPEC.f*OCPEC.timeStep-OCPEC.x) - gamma'*OCPEC.PHI;
[HAM_hessian_formula, ~] = hessian(HAM_formula, [OCPEC.x; OCPEC.u; OCPEC.p; OCPEC.w]);
HAM_hessian_FunObj = Function('HAM_hessian',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w, sigma, eta, lambda, gamma, OCPEC.s},...
    {HAM_hessian_formula},...
    {'x', 'u', 'p', 'w', 'sigma', 'eta', 'lambda', 'gamma', 's'},...
    {'HAM_hessian'});
FunObj.HAM_hessian = HAM_hessian_FunObj.map(nStages);

%% FB function and Jacobian
% basic FB function and Jacobian
a = SX.sym('a', 1, 1); % dual variable
b = SX.sym('b', 1, 1); % inequality
z = SX.sym('z', 1, 1); % perturbed parameter for FB function
sqrt_abz = sqrt(a^2 + b^2 + z^2);
PSI = sqrt_abz - a - b;
PSIa = a/sqrt_abz - 1;
PSIb = b/sqrt_abz - 1;
FB = Function('FB', {a, b, z}, {PSI}, {'a', 'b', 'z'}, {'PSI'});
FB_grad = Function('FB_grad',{a, b, z}, {PSIa, PSIb}, {'a', 'b', 'z'}, {'PSIa', 'PSIb'});

% FB function for G 
G = SX.sym('G', Dim.sigma, 1);
PSIg_formula = SX.sym('PSIg_formula', Dim.sigma, 1);
PSIgSigma_diagVec_formula = SX.sym('PSIgSigma_diagVec_formula', Dim.sigma, 1);% diagonal elements vector of PSIgSigma
PSIgG_diagVec_formula = SX.sym('PSIgG_diagVec_formula', Dim.sigma, 1);% diagonal elements vector of PSIgG
for i = 1 : Dim.sigma
    PSIg_formula(i, 1) = FB(sigma(i, 1), G(i, 1), z);
    [PSIgSigma_diagVec_formula(i, 1), PSIgG_diagVec_formula(i, 1)] = FB_grad(sigma(i, 1), G(i, 1), z);
end
FB_G_FunObj = Function('FB_G',...
    {sigma, G, z}, {PSIg_formula},...
    {'sigma', 'G', 'z'}, {'PSIg'});
FunObj.FB_G = FB_G_FunObj.map(nStages);
FB_G_grad_FunObj = Function('FB_G_grad',...
    {sigma, G, z}, {PSIgSigma_diagVec_formula, PSIgG_diagVec_formula},...
    {'sigma', 'G', 'z'}, {'PSIgSigma_diagVec', 'PSIgG_diagVec'});
FunObj.FB_G_grad = FB_G_grad_FunObj.map(nStages);

% FB function for PHI 
PHI = SX.sym('PHI', Dim.gamma, 1);
PSIphi_formula = SX.sym('PSIphi_formula', Dim.gamma, 1);
PSIphiGamma_diagVec_formula = SX.sym('PSIphiGamma_diagVec_formula', Dim.gamma, 1);% diagonal elements vector of PSIphiGamma
PSIphiPHI_diagVec_formula = SX.sym('PSIphiPHI_diagVec_formula', Dim.gamma, 1);% diagonal elements vector of PSIphiPHI
for i = 1 : Dim.gamma
    PSIphi_formula(i, 1) = FB(gamma(i, 1), PHI(i, 1), z);
    [PSIphiGamma_diagVec_formula(i, 1), PSIphiPHI_diagVec_formula(i, 1)] = FB_grad(gamma(i, 1), PHI(i, 1), z);
end
FB_PHI_FunObj = Function('FB_PHI',...
    {gamma, PHI, z}, {PSIphi_formula},...
    {'gamma', 'PHI', 'z'}, {'PSIphi'});
FunObj.FB_PHI = FB_PHI_FunObj.map(nStages);
FB_PHI_grad_FunObj = Function('FB_PHI_grad',...
    {gamma, PHI, z}, {PSIphiGamma_diagVec_formula, PSIphiPHI_diagVec_formula},...
    {'gamma', 'PHI', 'z'}, {'PSIphiGamma_diagVec', 'PSIphiPHI_diagVec'});
FunObj.FB_PHI_grad = FB_PHI_grad_FunObj.map(nStages);

%% KKT Residual: G_Fsb, PHI_Fsb, HAM_grad(HxT, HuT, HpT, HwT)
% G_Fsb 
nu_G = SX.sym('nu_G', 1, 1); % nonsingular regularization parameter nu_G

PSIg = SX.sym('PSIg', Dim.sigma, 1);
PSIgSigma_diagVec = SX.sym('PSIgSigma_diagVec', Dim.sigma, 1);% diagonal elements vector of PSIgSigma
PSIgG_diagVec = SX.sym('PSIgG_diagVec', Dim.sigma, 1);% diagonal elements vector of PSIgG

G_Fsb = -PSIg./(PSIgG_diagVec - nu_G * ones(Dim.sigma, 1));
G_Fsb_FunObj = Function('G_Fsb', {PSIg, PSIgG_diagVec, nu_G}, {G_Fsb},...
    {'PSIg', 'PSIgG_diagVec','nu_G'}, {'G_Fsb'});
FunObj.G_Fsb = G_Fsb_FunObj.map(nStages);

% PHI_Fsb
PSIphi = SX.sym('PSIphi', Dim.gamma, 1);
PSIphiGamma_diagVec = SX.sym('PSIphiGamma_diagVec', Dim.gamma, 1);% diagonal elements vector of PSIphiGamma
PSIphiPHI_diagVec = SX.sym('PSIphiPHI_diagVec', Dim.gamma, 1);% diagonal elements vector of PSIphiPHI

PHI_Fsb = -PSIphi./(PSIphiPHI_diagVec - nu_G * ones(Dim.gamma, 1));
PHI_Fsb_FunObj = Function('PHI_Fsb', {PSIphi, PSIphiPHI_diagVec, nu_G}, {PHI_Fsb},...
    {'PSIphi', 'PSIphiPHI_diagVec', 'nu_G'}, {'PHI_Fsb'});
FunObj.PHI_Fsb = PHI_Fsb_FunObj.map(nStages);

% Hamiltonian stage Jacobian: HxT, HuT, HpT, HwT
Lx = SX.sym('Lx', 1, Dim.x);
Lu = SX.sym('Lu', 1, Dim.u);
Lp = SX.sym('Lp', 1, Dim.p);
Lw = SX.sym('Lw', 1, Dim.w);

Gx = SX.sym('Gx', Dim.sigma, Dim.x);
Gu = SX.sym('Gu', Dim.sigma, Dim.u);
Gp = SX.sym('Gp', Dim.sigma, Dim.p);

Cx = SX.sym('Cx', Dim.eta, Dim.x);
Cu = SX.sym('Cu', Dim.eta, Dim.u);
Cp = SX.sym('Cp', Dim.eta, Dim.p);
Cw = SX.sym('Cw', Dim.eta, Dim.w);

Fx = SX.sym('Fx', Dim.lambda, Dim.x);
Fu = SX.sym('Fu', Dim.lambda, Dim.u);
Fp = SX.sym('Fp', Dim.lambda, Dim.p);

PHIp = SX.sym('PHIp', Dim.gamma, Dim.p);
PHIw = SX.sym('PHIw', Dim.gamma, Dim.w);

HxT = Lx' - Gx'*sigma + Cx'*eta + Fx'*lambda;
HuT = Lu' - Gu'*sigma + Cu'*eta + Fu'*lambda;
HpT = Lp' - Gp'*sigma + Cp'*eta + Fp'*lambda - PHIp'*gamma;
HwT = Lw'             + Cw'*eta +            - PHIw'*gamma;

HAM_grad_FunObj = Function('HAM_grad',...
    {sigma, eta, lambda, gamma,...
    Lx, Lu, Lp, Lw,...
    Gx, Gu, Gp,...
    Cx, Cu, Cp, Cw,...
    Fx, Fu, Fp,...
    PHIp, PHIw},...
    {HxT, HuT, HpT, HwT},...
    {'sigma', 'eta', 'lambda', 'gamma',...
    'Lx', 'Lu', 'Lp', 'Lw',...
    'Gx', 'Gu', 'Gp',...
    'Cx', 'Cu', 'Cp', 'Cw',...
    'Fx', 'Fu', 'Fp',...
    'PHIp', 'PHIw'},...
    {'HxT', 'HuT', 'HpT', 'HwT'});
FunObj.HAM_grad = HAM_grad_FunObj.map(nStages);

%% KKT diagnal matrix J 
nu_J = SX.sym('nu_J', 1, 1); % nonsingular regularization parameter nu_J
nu_H = SX.sym('nu_H', 1, 1); % nonsingular regularization parameter nu_H

d = - (PSIgSigma_diagVec - nu_G * ones(Dim.sigma, 1))./(PSIgG_diagVec - nu_G * ones(Dim.sigma, 1));
nu_J_vec = - nu_J * ones(Dim.eta + Dim.lambda, 1);
e = -(PSIphiGamma_diagVec - nu_G * ones(Dim.gamma, 1))./(PSIphiPHI_diagVec - nu_G * ones(Dim.gamma, 1));
diagVec = [d; nu_J_vec; e];
G_grad = [Gx, Gu, Gp, SX(Dim.sigma, Dim.w)];
C_grad = [Cx, Cu, Cp, Cw];
F_grad = [Fx, Fu, Fp, SX(Dim.lambda, Dim.w)];
PHI_grad = [SX(Dim.gamma, Dim.x + Dim.u), PHIp, PHIw];
Hessian = SX.sym('Hessian', Dim.Z, Dim.Z);

J = SX(Dim.Y, Dim.Y);
for i = 1 : Dim.Node(4)
    J(i, i) = diagVec(i);
end
J(1 : Dim.Node(4), Dim.Node(4) + 1 : Dim.Node(8)) = [-G_grad; C_grad; F_grad; -PHI_grad];
J(Dim.Node(4) + 1 : Dim.Node(8), :) = [-G_grad', C_grad', F_grad', -PHI_grad', Hessian + nu_H * SX.eye(Dim.Z)];

J_FunObj = Function('J',...
    {nu_G, nu_J, nu_H,...
    PSIgSigma_diagVec, PSIgG_diagVec, PSIphiGamma_diagVec, PSIphiPHI_diagVec,...
    Gx, Gu, Gp,...
    Cx, Cu, Cp, Cw,...
    Fx, Fu, Fp,...
    PHIp, PHIw,...
    Hessian},...
    {J},...
    {'nu_G', 'nu_J', 'nu_H',...
    'PSIgSigma_diagVec', 'PSIgG_diagVec', 'PSIphiGamma_diagVec', 'PSIphiPHI_diagVec',...
    'Gx', 'Gu', 'Gp',...
    'Cx', 'Cu', 'Cp', 'Cw',...
    'Fx', 'Fu', 'Fp',...
    'PHIp', 'PHIw',...
    'Hessian'},...
    {'J'});
FunObj.J = J_FunObj.map(nStages);

%% Line Search
% total cost directional derivative
dx = SX.sym('dx', Dim.x, 1);
du = SX.sym('du', Dim.u, 1);
dp = SX.sym('dp', Dim.p, 1);
dw = SX.sym('dw', Dim.w, 1);

totalCostDD = Lx * dx + Lu * du + Lp * dp + Lw * dw;
totalCostDD_FunObj = Function('totalCostDD', {Lx, Lu, Lp, Lw, dx, du, dp, dw}, {totalCostDD},...
    {'Lx', 'Lu', 'Lp', 'Lw', 'dx', 'du', 'dp', 'dw'}, {'totalCostDD'});
FunObj.totalCostDD = totalCostDD_FunObj.map(nStages);

%% Feasibility Restoration Phase
ZRef = SX.sym('ZRef', Dim.Z, 1); % reference primal variable in FRP cost function
ZWeight = SX.sym('ZWeight', Dim.Z, 1); % weight matrix in FRP cost function
Z = [OCPEC.x; OCPEC.u; OCPEC.p; OCPEC.w];
% cost function
FRP_L = 0.5 * (Z - ZRef)' * diag(ZWeight) * (Z - ZRef);
FRP_L_FunObj = Function('FRP_L',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w, ZRef, ZWeight},...
    {FRP_L},...
    {'x', 'u', 'p', 'w', 'ZRef', 'ZWeight'},...
    {'FRP_L'});
FunObj.FRP_L = FRP_L_FunObj.map(nStages);

% Jacobian
FRP_Lx = jacobian(FRP_L, OCPEC.x);
FRP_Lu = jacobian(FRP_L, OCPEC.u);
FRP_Lp = jacobian(FRP_L, OCPEC.p);
FRP_Lw = jacobian(FRP_L, OCPEC.w);

FRP_L_grad_FunObj = Function('FRP_L_grad',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w, ZRef, ZWeight},...
    {FRP_Lx, FRP_Lu, FRP_Lp, FRP_Lw},...
    {'x', 'u', 'p', 'w', 'ZRef', 'ZWeight'},...
    {'FRP_Lx', 'FRP_Lu', 'FRP_Lp', 'FRP_Lw'});
FunObj.FRP_L_grad = FRP_L_grad_FunObj.map(nStages);

% Hessian
[FRP_L_hessian, ~] = hessian(FRP_L, [OCPEC.x; OCPEC.u; OCPEC.p; OCPEC.w]);
FRP_L_hessian_FunObj = Function('FRP_L_hessian',...
    {OCPEC.x, OCPEC.u, OCPEC.p, OCPEC.w, ZRef, ZWeight},...
    {FRP_L_hessian},...
    {'x', 'u', 'p', 'w', 'ZRef', 'ZWeight'},...
    {'FRP_L_hessian'});
FunObj.FRP_L_hessian = FRP_L_hessian_FunObj.map(nStages);
end

