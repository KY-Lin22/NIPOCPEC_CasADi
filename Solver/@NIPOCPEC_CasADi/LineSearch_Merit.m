function [Var_LS, Info] = LineSearch_Merit(self, Var, Fun, Jac, beta, s, z,...
            KKT_Residual, KKT_Matrix, dY_k, mode, FRP)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
TimeStart = tic;

OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;

Dim = OCPEC.Dim;
employSOC = Option.employSecondOrderCorrection;

nu_D = Option.LineSearch.nu_D;
stepSize_Init = 1;
switch mode
    case 'Regular'
        rho = Option.LineSearch.rho;
        stepSize_Min = Option.LineSearch.stepSize_Min;
        stepSize_DecayRate = Option.LineSearch.stepSize_DecayRate;
    case 'FRP'
        rho = Option.FRP.rho;
        stepSize_Min = Option.FRP.stepSize_Min;
        stepSize_DecayRate = Option.FRP.stepSize_DecayRate;      
end

%% Some Evaluation Quantities of Previous Iterate
% extract search direction
dsigma_k  = dY_k(              1 : Dim.Node(1), :);
deta_k    = dY_k(Dim.Node(1) + 1 : Dim.Node(2), :);
dlambda_k = dY_k(Dim.Node(2) + 1 : Dim.Node(3), :);
dgamma_k  = dY_k(Dim.Node(3) + 1 : Dim.Node(4), :);
dx_k      = dY_k(Dim.Node(4) + 1 : Dim.Node(5), :);
du_k      = dY_k(Dim.Node(5) + 1 : Dim.Node(6), :);
dp_k      = dY_k(Dim.Node(6) + 1 : Dim.Node(7), :);
dw_k      = dY_k(Dim.Node(7) + 1 : Dim.Node(8), :);

% cost and its directional derivative
totalCost = sum(Fun.L);
totalCostDD = FunObj.totalCostDD(Jac.Lx, Jac.Lu, Jac.Lp, Jac.Lw, dx_k, du_k, dp_k, dw_k);
totalCostDD = sum(full(totalCostDD));

% constraint violation (L1 norm)
totalCstrVio_L1Norm = norm(reshape(Fun.PSIg, [], 1), 1) + norm(reshape(Fun.C, [], 1), 1)...
    + norm(reshape(Fun.F, [], 1), 1) + norm(reshape(Fun.PSIphi, [], 1), 1); 

% penalty parameter
beta_Trial = totalCostDD/((1 - rho) * totalCstrVio_L1Norm);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end

% merit and its directional derivative 
merit = totalCost + beta_k * totalCstrVio_L1Norm;
meritDD = totalCostDD - beta_k * totalCstrVio_L1Norm;

%% Backtracking Line Search
hasFoundNewIterate = false;

while ~hasFoundNewIterate
    %% Step 1: estimate trail stepsize, var and merit 
    stepSize_trial = max([stepSize_Init, stepSize_Min]);
    
    Var_trial.sigma  = Var.sigma  + stepSize_trial * dsigma_k;
    Var_trial.eta    = Var.eta    + stepSize_trial * deta_k;
    Var_trial.lambda = Var.lambda + stepSize_trial * dlambda_k;
    Var_trial.gamma  = Var.gamma  + stepSize_trial * dgamma_k;
    Var_trial.x      = Var.x      + stepSize_trial * dx_k;
    Var_trial.u      = Var.u      + stepSize_trial * du_k;
    Var_trial.p      = Var.p      + stepSize_trial * dp_k;
    Var_trial.w      = Var.w      + stepSize_trial * dw_k;    

    Fun_trial = self.FunctionEvaluation(Var_trial, s, z, mode, FRP);
    totalCost_trail = sum(Fun_trial.L);
    totalCstrVio_L1Norm_trial = norm(reshape(Fun_trial.PSIg, [], 1), 1) + norm(reshape(Fun_trial.C, [], 1), 1)...
    + norm(reshape(Fun_trial.F, [], 1), 1) + norm(reshape(Fun_trial.PSIphi, [], 1), 1);
    merit_trial = totalCost_trail + beta_k * totalCstrVio_L1Norm_trial;

    %% Step 2: checking sufficient decrease condition
    if merit_trial <= merit + stepSize_trial * nu_D * meritDD
        % return merit line search Var
        hasFoundNewIterate = true;
        LineSearchFailureFlag = false;
        VarType = 'MLS';         
        
    elseif (strcmp(mode, 'Regular')) && (employSOC) && (stepSize_trial == 1) && (totalCost_trail <= totalCost)
        % estimate second order correction (SOC) Var and merit
        Fun_full = Fun_trial;
        Var_trial = self.SecondOrderCorrection(Var, Jac, KKT_Residual, KKT_Matrix, Fun_full);        
        
        Fun_trial = self.FunctionEvaluation(Var_trial, s, z, 'Regular', []);
        totalCost_trail = sum(Fun_trial.L);
        totalCstrVio_L1Norm_trial = norm(reshape(Fun_trial.PSIg, [], 1), 1) + norm(reshape(Fun_trial.C, [], 1), 1)...
            + norm(reshape(Fun_trial.F, [], 1), 1) + norm(reshape(Fun_trial.PSIphi, [], 1), 1);
        merit_trial = totalCost_trail + beta_k * totalCstrVio_L1Norm_trial;
        
        if merit_trial <= merit + stepSize_trial * nu_D * meritDD
            % return SOC Var
            hasFoundNewIterate = true;
            LineSearchFailureFlag = false;
            VarType = 'SOC';             
        else
            % discard this SOC step and resume backtracking line search with smaller stepsize
            stepSize_Init = stepSize_DecayRate * stepSize_Init;
        end
    else
        % need to estimate a smaller stepsize
        stepSize_Init = stepSize_DecayRate * stepSize_Init;        
    end
    
    %% Step 3: checking min stepsize
    if (stepSize_trial == stepSize_Min)&&(~hasFoundNewIterate)
        LineSearchFailureFlag = true;
        VarType = 'MLS'; % return Var obtained by the stepSize_Min   
        break
    end    
  
end
%% Record Information
TimeElapsed = toc(TimeStart);

Info.VarType = VarType;
Var_LS = Var_trial;
Info.Fun = Fun_trial;
Info.merit = merit_trial;
Info.beta = beta_k;
Info.stepSize = stepSize_trial;

Info.FailureFlag = LineSearchFailureFlag;
Info.Time = TimeElapsed;

end

