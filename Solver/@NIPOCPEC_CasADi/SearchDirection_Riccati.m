function [dY, Info] = SearchDirection_Riccati(self, KKT_Residual, KKT_Matrix)
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here
TimeStart = tic;

OCPEC = self.OCPEC;
Option = self.Option;
Dim = OCPEC.Dim;
N = OCPEC.nStages;

% regroup KKT residual
T = [KKT_Residual.G_Fsb;...
     KKT_Residual.C_Fsb;...
     KKT_Residual.F_Fsb;...
     KKT_Residual.PHI_Fsb;...
     KKT_Residual.HxTlambdaNext;...
     KKT_Residual.HuT;...     
     KKT_Residual.HpT;...
     KKT_Residual.HwT];
%
J  = KKT_Matrix.J;
BL = KKT_Matrix.BL; 
BU = KKT_Matrix.BU;

%% backward recursion for Frak_A and Frak_b
Frak_A = cell(1, N);
Frak_b = zeros(Dim.Y, N);
% initialization: Frak_A_N and Frak_b_N
J_N = J(:, Dim.Y * (N - 1) + 1 : Dim.Y * N);
T_N = T(:, N);
Frak_A{N} = J_N;
Frak_b(:, N) = T_N;
for n = N - 1 : -1 : 1
    % load 
    J_n = J(:, Dim.Y * (n - 1) + 1 : Dim.Y * n);
    T_n = T(:, n);
    Frak_A_nNext = Frak_A{n + 1};
    Frak_b_nNext = Frak_b(:, n + 1);
    % compute Frak_A_n and Frak_b_n
    switch Option.linearSystemSolver
        case 'linsolve_Sym_dense'
            % solve linear system using linsolve(dense) with 'SYM' option
            invFrak_A_nNext_BL = zeros(Dim.Y, Dim.Y);
            opts.SYM = true;
            invFrak_A_nNext_BL(:, Dim.Node(4) + 1 : Dim.Node(5)) = linsolve(Frak_A_nNext, BL(:, Dim.Node(4) + 1 : Dim.Node(5)), opts);
            invFrak_A_nNext_Frak_b_nNext = linsolve(Frak_A_nNext, Frak_b_nNext, opts);
            Frak_A_n = J_n - BU * invFrak_A_nNext_BL;
            Frak_b_n = T_n - BU * invFrak_A_nNext_Frak_b_nNext;
        case 'mldivide_dense'
            % solve linear system using mldivide(dense)
            Frak_A_n = J_n - BU * (Frak_A_nNext\BL);
            Frak_b_n = T_n - BU * (Frak_A_nNext\Frak_b_nNext);
        case 'mldivide_sparse'
            % solve linear system using mldivide(sparse)
            Frak_A_n = sparse(J_n) - sparse(BU) * (sparse(Frak_A_nNext)\BL);
            Frak_b_n = T_n - sparse(BU) * (sparse(Frak_A_nNext)\Frak_b_nNext);            
        otherwise
            error('specified method is not supported')
    end
    % record
    Frak_A{n} = Frak_A_n;
    Frak_b(:, n) = Frak_b_n;
end

%% forward recursion for dY
dY = zeros(Dim.Y, N);
dY_nPrev = zeros(Dim.Y, 1);
for n = 1 : N
    %
    Frak_A_n = Frak_A{n};
    Frak_b_n = Frak_b(:, n);
    % compute search direction dY_n
    switch Option.linearSystemSolver
        case 'linsolve_Sym_dense'
            % solve linear system using linsolve(dense) with 'SYM' option
            opts.SYM = true;
            dY(:, n) = linsolve(Frak_A_n, -(Frak_b_n + BL * dY_nPrev), opts);
        case 'mldivide_dense'
            % solve linear system using mldivide(dense)
            dY(:, n) = - Frak_A_n\(Frak_b_n + BL * dY_nPrev);
        case 'mldivide_sparse'
            % solve linear system using mldivide(sparse)
            dY(:, n) = - sparse(Frak_A_n)\(Frak_b_n + sparse(BL) * dY_nPrev);
        otherwise
            error('specified method is not supported')
    end
    dY_nPrev = dY(:, n);
end

TimeElapsed = toc(TimeStart);
Info.Time = TimeElapsed;
end

