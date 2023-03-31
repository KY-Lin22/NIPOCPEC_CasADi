clc
clear all
%% problem formulation
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

timeStep = 0.01;
nStages = 100;
sInit = 1e-8;
sEnd = 1e-8; 

zEnd = 1e-3;

x_Dim = 2;
tau_Dim = 1;
p_Dim = 1;
w_Dim = p_Dim; % auxiliary var
s_Dim = 1; % regularization parameter
x = SX.sym('x', x_Dim);
tau = SX.sym('tau', tau_Dim);
p = SX.sym('p', p_Dim);
w = SX.sym('w', w_Dim);
s = SX.sym('s', s_Dim);
InitState = [-1/2; -1];
RefState = [0; 0];
EndState = [0; 0];
% dynamics
f = [1, -3; -8, 10] * x + [-3; -1] * p + [4; 8] * tau; % xDot = f(x, tau, p)
f_Fun = Function('f_Fun', {x, tau, p}, {f}, {'x', 'tau', 'p'}, {'f'});
eqlbm.l = -1;
eqlbm.u = 1;
eqlbm.K = [1, -3] * x + 5 * p + 3 * tau;
K_Fun = Function('K_Fun', {x, tau, p}, {eqlbm.K}, {'x', 'tau', 'p'}, {'K'});
% reformulate equilibrium dynamics as a set of inequality and equality constriants using Scholtes reformulation
BVI = [p - eqlbm.l;...
    eqlbm.u - p;...%  p in [l, u]   
    w - eqlbm.K;...% auxiliary variable
    s - (p - eqlbm.l) * w;...
    s + (eqlbm.u - p) * w];% regularization
BVI_Fun = Function('BVI_Fun', {x, tau, p, w, s}, {BVI}, {'x', 'tau', 'p', 'w', 's'}, {'BVI'});
lbg_BVI = [0; 0; 0; 0; 0];
ubg_BVI = [Inf; Inf; 0; Inf; Inf];
% cost function
xWeight = [20; 20];
tauWeight = 1;
L_stageCost = 0.5*(x - RefState)'*diag(xWeight)*(x - RefState) + 0.5*tau'*diag(tauWeight)*tau;
L_stageCost_Fun = Function('L_stageCost_Fun', {x, tau, p, w}, {L_stageCost}, {'x', 'tau', 'p', 'w'}, {'L_stageCost'});
L_terminalCost = 0.5*(x - EndState)'*diag(xWeight)*(x - EndState);
L_terminalCost_Fun = Function('L_terminalCost_Fun', {x, p}, {L_terminalCost}, {'x', 'p'}, {'L_terminalCost'});
% formulate NLP 
X = SX.sym('X', x_Dim, nStages); % state variable
TAU = SX.sym('TAU', tau_Dim, nStages); % control variable
P = SX.sym('P', p_Dim, nStages); % algebraic variable
W = SX.sym('W', w_Dim, nStages); % auxilary variable
S = SX.sym('S', s_Dim, 1);% regularizaton parameter
x_Max = [2; 2];
x_Min = [-2; -2];
tau_Max = 2;
tau_Min = -2;
lbx = -Inf * ones(x_Dim + tau_Dim + p_Dim + w_Dim, nStages);
ubx = Inf * ones(x_Dim + tau_Dim + p_Dim + w_Dim, nStages);
lbx(1 : x_Dim + tau_Dim, :) = repmat([x_Min; tau_Min], 1, nStages);
ubx(1 : x_Dim + tau_Dim, :) = repmat([x_Max; tau_Max], 1, nStages);

L = 0; % init cost function
g_Dim = size([f; BVI], 1);
g = SX.sym('g', g_Dim, nStages); % constraint function
lbg = zeros(g_Dim, nStages);
ubg = zeros(g_Dim, nStages);
lbg(size(f, 1) + 1 : end, :) = repmat(lbg_BVI, 1, nStages);
ubg(size(f, 1) + 1 : end, :) = repmat(ubg_BVI, 1, nStages);
for n = 1 : nStages
    if n == 1
        x_nPrev = InitState;
    else
        x_nPrev = X(:, n - 1);
    end    
    x_n = X(:, n);
    tau_n = TAU(:, n);
    p_n = P(:, n);
    w_n = W(:, n);
    % cost function
    L_n = L_stageCost_Fun(x_n, tau_n, p_n, w_n);
    L = L + L_n*timeStep;    
    if n == nStages
        L_terminal = L_terminalCost_Fun(x_n, p_n);
        L = L + L_terminal;
    end 
    % discretize dynamics by implicit Euler method
    F_n = x_nPrev - x_n + timeStep * f_Fun(x_n, tau_n, p_n);   
    % reformulated equilibrium constraint
    BVI_n = BVI_Fun(x_n, tau_n, p_n, w_n, S); 
    % constraint function
    g(:, n) = [F_n;...
        BVI_n];     
end
OCPEC = struct('x', reshape([X; TAU; P; W], (x_Dim + tau_Dim + p_Dim + w_Dim) * nStages, 1),...
    'f', L,...    
    'g', reshape(g, g_Dim * nStages, 1),...
    'p', S);
args.lbx = reshape(lbx, [], 1);
args.ubx = reshape(ubx, [], 1);
args.lbg = reshape(lbg, [], 1);
args.ubg = reshape(ubg, [], 1);
% option (solver.print_options)
Option = struct;
Option.print_time = false;
Option.ipopt.max_iter = 2000;
Option.ipopt.tol = 1e-4;
Option.ipopt.mu_target = 0.5 * (zEnd)^2;
Option.ipopt.print_level = 0;
% create solver
solver = nlpsol('solver', 'ipopt', OCPEC, Option);

%% homotopy
test_Num = 20;
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
    % initialize x0 and s  
    x0 = zeros(x_Dim + tau_Dim + p_Dim + w_Dim, nStages);
    x0(x_Dim + 1 : x_Dim + tau_Dim, :) = randn(tau_Dim, nStages);    
    args.x0 = reshape(x0, (x_Dim + tau_Dim + p_Dim + w_Dim) * nStages, 1);
    args.p = sInit;% 

    [solution, Info] = homotopy_solve(solver, args, sEnd);
    TestRecord.InitialGuess{i, 1} = args.x0;
    TestRecord.solution{i, 1} = solution;

    if Info.homotopy_status == 1   
        successCase = successCase + 1;
        TestRecord.iterNum(successCase, 1) = Info.iterNum;
        TestRecord.totalTime(successCase, 1) = Info.totalTime;       
        TestRecord.cost(successCase, 1) = full(solution.f);    
        
        % compute constraint violation
        XTAUPW_Opt = reshape(full(solution.x), (x_Dim + tau_Dim + p_Dim + w_Dim), nStages);
        x_Opt   = XTAUPW_Opt(1                           : x_Dim, :);
        tau_Opt = XTAUPW_Opt(1 + x_Dim                   : x_Dim + tau_Dim, :);
        p_Opt   = XTAUPW_Opt(1 + x_Dim + tau_Dim         : x_Dim + tau_Dim + p_Dim, :);
        w_Opt   = XTAUPW_Opt(1 + x_Dim + tau_Dim + p_Dim : end, :);
        cstr = reshape(full(solution.g), g_Dim, nStages);
        EQ = cstr(5, :); % w - K = 0
        DynF = cstr(1 : x_Dim, :); % 1 : 2 
        G = [[x_Opt; tau_Opt] - repmat([x_Min; tau_Min], 1, nStages);...
            repmat([x_Max; tau_Max], 1, nStages) - [x_Opt; tau_Opt]];         
        K_Fun_map = K_Fun.map(nStages);
        K_value = K_Fun_map(x_Opt, tau_Opt, p_Opt);
        K_value = full(K_value);
        lpu = cstr(x_Dim + 1: x_Dim + 2 * p_Dim, :);
        G_residual = min([zeros(2 * (x_Dim + tau_Dim) * nStages, 1), reshape(G, [], 1)], [], 2);
        pK_residual = min([zeros(2 * p_Dim * nStages, 1), reshape(lpu, [], 1)], [], 2); % 3 : 4
        
        Complementarity_pK = zeros(p_Dim, nStages);
        for n = 1 : nStages
            for j = 1 : p_Dim
                l_Vio = max(0, eqlbm.l(j) - p_Opt(j, n));
                K_l_VioScale = min(1, max(0, p_Opt(j, n) - eqlbm.l(j)));
                l_ComlVio = max(l_Vio, K_l_VioScale * max(K_value(j, n), 0));
                u_Vio = max(0, p_Opt(j, n) - eqlbm.u(j));
                K_u_VioScale = min(1, max(0, eqlbm.u(j) - p_Opt(j, n)));
                u_ComlVio = max(u_Vio, K_u_VioScale * max(-K_value(j, n), 0));
                Complementarity_pK(j, n) = max(l_ComlVio, u_ComlVio);
            end
        end
    
        TestRecord.eqCstr(successCase, 1) = max([norm(reshape(EQ, [], 1), Inf), norm(reshape(DynF, [], 1), Inf)]);
        TestRecord.ineqCstr(successCase, 1) = max(norm(G_residual, Inf), norm(reshape(pK_residual, [], 1), Inf));
        TestRecord.compCstr(successCase, 1)  = norm(reshape(Complementarity_pK, [], 1), Inf);
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

%%

% timeAxis = 0 : timeStep : nStages * timeStep;
% 
% K_Fun_map = K_Fun.map(nStages);
% K_value = K_Fun_map(x_Opt, tau_Opt, p_Opt);
% K_value = full(K_value);
% figure(111)
% subplot(3,1,1)
% plot(timeAxis, [InitState(1), x_Opt(1, :)], 'r',...
%      timeAxis, [InitState(2), x_Opt(2, :)], 'g', 'LineWidth',1.2)
% legend('x1', 'x2')
% xlabel('time(s)')
% title('system state')
% 
% subplot(3,1,2)
% plot(timeAxis(2:end), tau_Opt(1,:), 'LineWidth', 1.2)
% xlabel('time(s)')
% title('control')
% 
% subplot(3,1,3)
% plot(timeAxis(2:end), p_Opt(1, :), 'k',...
%      timeAxis(2:end), K_value(1, :), 'b', 'LineWidth', 1.2)
% legend('p', 'K') 
% xlabel('time(s)')
% title('equilibrium dynamics')