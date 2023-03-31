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

solver.Option.sInit = 1e-8;
solver.Option.sEnd  = 1e-8;

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
