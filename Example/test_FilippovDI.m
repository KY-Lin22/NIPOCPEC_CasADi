clear all
clc
delete Gen_InitialGuess.mat

%% Problem Formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
% FilippovDI: an example from Stewart's paper:
% Optimal control of systems with discontinuous differential equations

timeStep = 0.02;
nStages = 100;
InitState = -1;
EndState = 5/3;
RefState = 0;
x_Dim = 1;
u_Dim = 0;
p_Dim = 1;
x = SX.sym('x', x_Dim, 1);
u = SX.sym('u', u_Dim, 1);
p = SX.sym('p', p_Dim, 1); % variable y in the Stewart's paper

% cost function
L_T = (x - EndState)^2;
L_S = (x - RefState)^2;

% inequality constraint
x_Max = 10000;
x_Min = -10000;
G = [x_Max - x;...
    x - x_Min];

% equality constraint
C = [];

% dynamics
f1 = 1;% switch function > 0
f2 = 3; % switch function < 0
f = f1*(1 - p) + f2*p;

% equilibrium constraint
K = x;
lbp = 0;
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
solver.Option.Tolerance.KKT_Error_Total = 1e-4;
solver.Option.Tolerance.KKT_Error_Feasibility = 1e-6;
solver.Option.Tolerance.KKT_Error_Stationarity = 1e-6;

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

% plot solution
timeAxis = 0 : timeStep : nStages * timeStep;

f_value = solver.FunObj.f(solution.x, solution.u, solution.p);
f_value = full(f_value);

figure(112)
subplot(3,2,1)
plot(timeAxis, [InitState, solution.x], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time(s)')
title('system state')

subplot(3,2,2)
plot(timeAxis(2:end), f_value, 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(3,2,3)
plot(timeAxis(2:end), solution.p, 'k', 'LineWidth', 1.2)
legend('y') 
xlabel('time(s)')
title('smoothing function')

subplot(3,2,4)
plot(timeAxis(2:end), solution.x, 'g', 'LineWidth', 1.2)
legend('z := x') 
xlabel('time(s)')
title('switch function')

subplot(3,2,5)
plot(timeAxis(2:end), solution.x, 'k',...
    timeAxis(2:end), solution.p, 'b', 'LineWidth', 1.2)
legend('z(K)', 'y(p)') 
xlabel('time(s)')
title('checking BVI')