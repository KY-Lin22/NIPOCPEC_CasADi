function showInfo(self)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%% problem information
OCPEC = self.OCPEC;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;
timeStep = OCPEC.timeStep;
disp('*--------------------------- OCPEC Problem Information ----------------------------*')
% time
disp(['timeHorizon: ...', num2str(nStages * timeStep), ' s'])
disp(['- Stage: .......', num2str(nStages),])
disp(['- timeStep: ....', num2str(timeStep),' s '])
% number of variable
disp(['Number of Total Variables: ..............', num2str(Dim.Y), ' * ', num2str(nStages), ' = ', num2str(Dim.Y * nStages)])
disp(['Number of Primal Variable: ..............', num2str(Dim.Z), ' * ', num2str(nStages), ' = ', num2str(Dim.Z * nStages)])
disp(['- x(optimal variable): ..................', num2str(Dim.x), ' * ', num2str(nStages), ' = ', num2str(Dim.x * nStages), '; '])
disp(['- u(control variable): ..................', num2str(Dim.u), ' * ', num2str(nStages), ' = ', num2str(Dim.u * nStages), '; '])
disp(['- p(algebraic state): ...................', num2str(Dim.p), ' * ', num2str(nStages), ' = ', num2str(Dim.p * nStages), '; '])
disp(['- w(auxiliary variable): ................', num2str(Dim.w), ' * ', num2str(nStages), ' = ', num2str(Dim.w * nStages), '; '])

disp(['Number of Dual Variable(constraint): ....', num2str(Dim.LAMBDA), ' * ', num2str(nStages), ' = ', num2str(Dim.LAMBDA * nStages)])
disp(['- sigma(inequality): ....................', num2str(Dim.sigma),  ' * ', num2str(nStages), ' = ', num2str(Dim.sigma * nStages), '; '])
disp(['- eta(equality): ........................', num2str(Dim.eta),    ' * ', num2str(nStages), ' = ', num2str(Dim.eta * nStages), '; '])
disp(['- lambda(state equation): ...............', num2str(Dim.lambda), ' * ', num2str(nStages), ' = ', num2str(Dim.lambda * nStages), '; '])
disp(['- gamma(equilibrium constraint): ........', num2str(Dim.gamma),  ' * ', num2str(nStages), ' = ', num2str(Dim.gamma * nStages)])
disp('')

%% solver option information
Option = self.Option;
disp('*----------------------- NIPOCPEC Solver Option Information -----------------------*')
disp('1. Basic Options')
disp(['- maxIterNum: .............................', num2str(Option.maxIterNum)])
disp(['- KKT Error Tolerance(Total): .............', num2str(Option.Tolerance.KKT_Error_Total)])
disp(['                     (Feasibility): .......', num2str(Option.Tolerance.KKT_Error_Feasibility)])
disp(['                     (Stationarity): ......', num2str(Option.Tolerance.KKT_Error_Stationarity)])
disp('2. Options for Function and Jacobian Evaluation')
disp(['- Hessian Approximation Method: ...........', Option.HessianApproximation])
disp(['- Singularity Regular Parameter(J): .......', num2str(Option.RegularParam.nu_J)])
disp(['                               (G): .......', num2str(Option.RegularParam.nu_G)])
disp(['                               (H): .......', num2str(Option.RegularParam.nu_H)])
disp('3. Options for Search Direction Evaluation')
disp(['- Solve Linear System Method: .............', Option.linearSystemSolver])
disp('4. Options for Line Search')
disp(['- Employ Second Order Correction: .........', mat2str(Option.employSecondOrderCorrection)]);
disp(['- stepSize_Min: ...........................', num2str(Option.LineSearch.stepSize_Min)])
disp('5. Options for Feasibility Restoration Phase')
disp(['- Employ Feasibility Restoration Phase: ...', mat2str(Option.employFeasibilityRestorationPhase)])
disp(['- maxIterNum: .............................', num2str(Option.FRP.maxIterNum)]);
disp(['- stepSize_Min: ...........................', num2str(Option.FRP.stepSize_Min)]);
disp('6. Options for Perturbed Parameter')
disp(['- sInit: ..................................', num2str(Option.sInit)])
disp(['- sEnd: ...................................', num2str(Option.sEnd)])
disp(['- zInit: ..................................', num2str(Option.zInit)])
disp(['- zEnd: ...................................', num2str(Option.zEnd)])
disp('*----------------------------------------------------------------------------------*')

end

