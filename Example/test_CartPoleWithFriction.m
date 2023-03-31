clear all
clc
delete Gen_InitialGuess.mat

%% Problem Formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 400;
x_Dim = 4;
u_Dim = 1;
p_Dim = 1;
InitState = [1; 0/180*pi; 0; 0];
EndState = [1; 180/180*pi; 0; 0];
RefState = [1; 180/180*pi; 0; 0];
x = SX.sym('x', x_Dim, 1);
u = SX.sym('u', u_Dim, 1);
p = SX.sym('p', p_Dim, 1);
% cost function
xWeight_T = [1; 100; 10; 20];
xWeight_S = [1; 100; 1; 1];
uWeight = 1;
L_T = 0.5*(x - EndState)'*diag(xWeight_T)*(x - EndState);
L_S = 0.5*(x - RefState)'*diag(xWeight_S)*(x - RefState) + 0.5*u'*diag(uWeight)*u;
% inequality constraint
x_Max = [5; 240/180*pi; 20; 20];
x_Min = [0; -240/180*pi; -20; -20];
u_Max = 30;
u_Min = -30;
G = [x_Max - x;...
    x - x_Min;...
    u_Max - u;...
    u - u_Min];
% equality constraint
C = [];
% dynamics
mass = [1; 0.1];
linkLength = 1;
g = 9.8;
M = [mass(1) + mass(2),                 mass(2) * linkLength * cos(x(2));...
     mass(2) * linkLength * cos(x(2)),  mass(2) * linkLength^2];
C_Mat = [0,   -mass(2) * linkLength * x(4) * sin(x(2));...
     0,   0]; 
G_Mat = [0;...
     -mass(2) * g * linkLength * sin(x(2))];
Bu = [u(1);...
      0]; 
P = [p(1);...
     0]; % friction bewteen cart and ground  
H = G_Mat + Bu + P - C_Mat * [x(3); x(4)];
f = [x(3:4);...
    inv(M)*H];% xDot = f(x, tau, p)
% equilibrium constraint
K = x(3);
lbp = -2;
ubp = 2;
% formulate OCPEC
OCPEC = struct('x', x, 'u', u, 'p', p,...
    'L_T', L_T, 'L_S', L_S,...
    'G', G, 'C', C, 'f', f, 'K', K,...
    'lbp', lbp, 'ubp', ubp,...
    'timeStep', timeStep, 'nStages', nStages, 'InitState', InitState);

%% create Solver
% create solver object
solver = NIPOCPEC_CasADi(OCPEC);

%% setting solver option
solver.Option.printLevel = 2;
solver.Option.maxIterNum = 2000;
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

solver.Option.sInit = 1e-3;
solver.Option.sEnd  = 1e-3;

solver.showInfo();

%% solve
test_Num = 20;
solver.Option.printLevel = 0;
TestRecord.InitialGuess = cell(test_Num, 1);
TestRecord.solution = cell(test_Num, 1);
successCase = 0;
TestRecord.iterNum = zeros(test_Num, 1);
TestRecord.totalTime = zeros(test_Num, 1);
TestRecord.cost = zeros(test_Num, 1);
TestRecord.eqCstr = zeros(test_Num, 1);
TestRecord.ineqCstr = zeros(test_Num, 1);
TestRecord.compCstr = zeros(test_Num, 1);

for i = 1 : test_Num    
    solver.generateInitialGuess();
    Gen_InitialGuess = load('Gen_InitialGuess.mat');

    [solution, Info] = solver.solveOCPEC(Gen_InitialGuess.Var);
    TestRecord.InitialGuess{i, 1} = Gen_InitialGuess.Var;
    TestRecord.solution{i, 1} = solution;

    if Info.iterProcess.terminalStatus == 1   
        successCase = successCase + 1;
        TestRecord.iterNum(successCase, 1) = Info.iterProcess.iterNum;
        TestRecord.totalTime(successCase, 1) = Info.iterProcess.Time.total;
        
        TestRecord.cost(successCase, 1) = Info.solutionMsg.totalCost;    
        
        TestRecord.eqCstr(successCase, 1) = max([Info.solutionMsg.r_eq_C, Info.solutionMsg.r_eq_F]);
        TestRecord.ineqCstr(successCase, 1) = max([Info.solutionMsg.r_ineq_G, Info.solutionMsg.r_eqlb_ineq]);
        TestRecord.compCstr(successCase, 1)  = Info.solutionMsg.r_eqlb_comp;
    end
    disp(['success / Test No.: ', num2str(successCase), ' / ', num2str(i)])
end  
save('Test_NIP_Data.mat', 'TestRecord');
% show result 
disp('Test Result')
disp(['success/total: ', num2str(successCase), '/', num2str(test_Num)])
disp(['time per iter: ', num2str(1000 * sum(TestRecord.totalTime) /sum(TestRecord.iterNum), '%10.3f'), ' ms/Iter' ])
disp(['iterations: ', num2str(sum(TestRecord.iterNum) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(TestRecord.iterNum(1 : successCase, 1))), '(max); ',...
    num2str(min(TestRecord.iterNum(1 : successCase, 1))), '(min)'])
disp(['totalTime [s]: ', num2str(sum(TestRecord.totalTime) / successCase, '%10.3f'), '(mean); ',...
    num2str(max(TestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(max); ',...
    num2str(min(TestRecord.totalTime(1 : successCase, 1)), '%10.3f'), '(min)'])
disp(['cost: ', num2str(sum(TestRecord.cost) / successCase , '%10.3f'), '(mean); ',...
    num2str(max(TestRecord.cost(1 : successCase, 1)), '%10.3f'), '(max); ',...
    num2str(min(TestRecord.cost(1 : successCase, 1)), '%10.3f'),'(min)'])
disp(['eqCstr: ', num2str(sum(TestRecord.eqCstr) / successCase,'%10.3e'), '(mean); ',...
    num2str(max(TestRecord.eqCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(TestRecord.eqCstr(1 : successCase, 1)), '%10.3e'),'(min)'])
disp(['ineqCstr: ', num2str(sum(TestRecord.ineqCstr) / successCase, '%10.3e'), '(mean); ',...
    num2str(max(TestRecord.ineqCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(TestRecord.ineqCstr(1 : successCase, 1)), '%10.3e'),'(min)'])
disp(['compCstr: ', num2str(sum(TestRecord.compCstr) / successCase, '%10.3e'), '(mean); ',...
    num2str(max(TestRecord.compCstr(1 : successCase, 1)), '%10.3e'), '(max); ',...
    num2str(min(TestRecord.compCstr(1 : successCase, 1)), '%10.3e'),'(min)'])

%% show result
% iteration process information
solver.showResult(Info)

timeAxis = 0 : timeStep : nStages * timeStep;
K_value = solver.FunObj.K(solution.x, solution.u, solution.p); 
K_value = full(K_value);

figure(111)
subplot(4,1,1)
plot(timeAxis, [InitState(1), solution.x(1, :)], 'r',...
     timeAxis, [InitState(2), solution.x(2, :)], 'g','LineWidth',1.2)
legend('cart', 'pole')
xlabel('time [s]')
title('position')

subplot(4,1,2)
plot(timeAxis, [InitState(3), solution.x(3, :)], 'r',...
     timeAxis, [InitState(4), solution.x(4, :)], 'g', 'LineWidth',1.2)
xlabel('time [s]')
title('velocity')

subplot(4,1,3)
plot(timeAxis(2:end), solution.u(1,:), 'LineWidth', 1.2)
xlabel('time [s]')
title('control')

subplot(4,1,4)
plot(timeAxis(2:end), solution.p(1, :), 'k',...
     timeAxis(2:end), K_value(1, :), 'b', 'LineWidth', 1.2)
legend('friction', 'cart vel') 
xlabel('time [s]')
title('equilibrium dynamics')
