function Option = createOption(self)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Basic Options
Option.printLevel = 2; % 0: print brief results;  
                       % 1: print detailed results
                       % 2: print detailed results and each iteration's information
Option.maxIterNum = 500;

Option.Tolerance.KKT_Error_Total = 1e-2;
Option.Tolerance.KKT_Error_Feasibility = 1e-3;
Option.Tolerance.KKT_Error_Stationarity = 1e-3;

%% Options for Function and Jacobian Evaluation
Option.HessianApproximation = 'CostFunction'; %  'Exact', 'CostFunction', 'GaussNewton',
% singularity regular parameter for KKT matrix 
Option.RegularParam.nu_J = 1e-7; % rank-deficiency of equality-type constraint
Option.RegularParam.nu_G = 1e-8; % non-negative definite of diagonal matrix related to inequality-type constraint
Option.RegularParam.nu_H = 0; % non-positive definite of Hessian matrix

%% Options for Search Direction Evaluation
Option.linearSystemSolver = 'linsolve_Sym_dense'; % 'linsolve_Sym_dense', 'mldivide_dense', 'mldivide_sparse'

%% Options for Line Search
Option.employSecondOrderCorrection = true;

Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 0.01;
Option.LineSearch.stepSize_DecayRate = 0.7;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4

%% Options for Feasibility Restoration Phase
Option.employFeasibilityRestorationPhase = true;

Option.FRP.maxIterNum = 15;
Option.FRP.nu_Z = 1e-6; % Z weight matrix scaling parameter         
Option.FRP.nu_M = 0.9; % desired reduction scale of feasibility

Option.FRP.betaInit = 10;
Option.FRP.rho = 0.5;
Option.FRP.stepSize_Min = 0.0001;
Option.FRP.stepSize_DecayRate = 0.5;

Option.FRP.employLeastSquareMinNorm = true;
Option.FRP.LAMBDA_threshold = 1000;

%% Options for Perturbed Parameter
% init and end perturbed parameter z in smoothing FB function
Option.zInit = 0.1; 
Option.zEnd = 0.001;

% init and end perturbed parameter s in perturbed equilibrium dynamics
Option.sInit = 0.1; 
Option.sEnd  = 0.001;

% update s and z
Option.kappa_s_times = 0.2;
Option.kappa_s_exp = 1.5;
Option.kappa_z_times = 0.2;
Option.kappa_z_exp = 1.5;

% threshold weight
Option.kappa_F = 2;
end

