clear all
clc
delete Gen_InitialGuess.mat

%% Problem Formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
timeStep = 0.01;
nStages = 100;
InitState = [-1/2; -1];
EndState = [0; 0];
RefState = [0; 0];
x_Dim = 2;
u_Dim = 1;
p_Dim = 1;
x = SX.sym('x', x_Dim, 1);
u = SX.sym('u', u_Dim, 1);
p = SX.sym('p', p_Dim, 1);

% cost function
xWeight = [20; 20];
uWeight = 1;

L_T = 0.5*(x - EndState)'*diag(xWeight)*(x - EndState);
L_S = 0.5*(x - RefState)'*diag(xWeight)*(x - RefState) + 0.5*u'*diag(uWeight)*u;

% inequality constraint
x_Max = [2; 2];
x_Min = [-2; -2];
u_Max = 2;
u_Min = -2;
G = [x_Max - x;...
    x - x_Min;...
    u_Max - u;...
    u - u_Min];

% equality constraint
C = [];

% dynamics
f = [1, -3; -8, 10] * x + [-3; -1] * p + [4; 8] * u; 

% equilibrium constraint
K = [1, -3] * x + 5 * p + 3 * u;
lbp = -1;
ubp = 1;

% formulate OCPEC
OCPEC.x = x;
OCPEC.u = u;
OCPEC.p = p;
OCPEC.L_T = L_T;
OCPEC.L_S = L_S;
OCPEC.G = G;
OCPEC.C = C;
OCPEC.f = f;
OCPEC.K = K;
OCPEC.lbp = lbp;
OCPEC.ubp = ubp;
OCPEC.timeStep = timeStep;
OCPEC.nStages = nStages;
OCPEC.InitState = InitState;

%% create Solver
% create solver object
solver = NIPOCPEC_CasADi(OCPEC);

solver.showInfo();

solver.generateInitialGuess();
Gen_InitialGuess = load('Gen_InitialGuess.mat');
Var_Init = Gen_InitialGuess.Var;

%% solving MPEC
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 500;
solver.Option.Tolerance.KKT_Error_Total = 1e-2;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-4;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-4;

solver.Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton'
solver.Option.RegularParam.nu_J = 1e-7;
solver.Option.RegularParam.nu_G = 1e-7;
solver.Option.RegularParam.nu_H = 0;
solver.Option.linearSystemSolver = 'linsolve_Sym_dense'; % 'linsolve_Sym_dense', 'mldivide_dense', 'mldivide_sparse'

solver.Option.employSecondOrderCorrection = false;
solver.Option.LineSearch.stepSize_Min = 0.001;
solver.Option.employFeasibilityRestorationPhase = true;

solver.Option.zInit = 1e-1; 
solver.Option.zEnd  = 1e-3;
solver.Option.sInit = 1e-1;
solver.Option.sEnd  = 1e-3;

tic
[solution, Info] = solver.solveOCPEC(Var_Init);
toc

%% show result
% iteration process information
solver.showResult(Info)

timeAxis = 0 : timeStep : nStages * timeStep;
K_value = solver.FunObj.K(solution.x, solution.u, solution.p); 
K_value = full(K_value);

figure(111)
subplot(3,1,1)
plot(timeAxis, [InitState(1), solution.x(1, :)], 'r',...
     timeAxis, [InitState(2), solution.x(2, :)], 'g', 'LineWidth',1.2)
legend('x1', 'x2')
xlabel('time(s)')
title('system state')

subplot(3,1,2)
plot(timeAxis(2:end), solution.u(1,:), 'LineWidth', 1.2)
xlabel('time(s)')
title('control')

subplot(3,1,3)
plot(timeAxis(2:end), solution.p(1, :), 'k',...
     timeAxis(2:end), K_value(1, :), 'b', 'LineWidth', 1.2)
legend('p', 'K') 
xlabel('time(s)')
title('equilibrium dynamics')
