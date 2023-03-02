function Hessian = HessianEvaluation(self, Var, Jac, s, mode, FRP)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

OCPEC = self.OCPEC;
FunObj = self.FunObj;
Dim = OCPEC.Dim;
nStages = OCPEC.nStages;

switch mode
    case 'Regular'
        %% regular iteration routine
        % terminal 
        L_T_hessian = FunObj.L_T_hessian(Var.x(:, end), Var.u(:, end), Var.p(:, end), Var.w(:, end));        
        % stage
        switch self.Option.HessianApproximation
            case 'Exact'
                HAM_hessian = FunObj.HAM_hessian(Var.x, Var.u, Var.p, Var.w,...
                    Var.sigma, Var.eta, Var.lambda, Var.gamma, s);  
                Hessian = full(HAM_hessian);
            case 'CostFunction'
                L_S_hessian = FunObj.L_S_hessian(Var.x, Var.u, Var.p, Var.w);
                Hessian = full(L_S_hessian);
            case 'GaussNewton'
                error('specified method to compute Hessian is not supported')
            otherwise
                error('specified method to compute Hessian is not supported')
        end
        % merge
        Hessian(:, (nStages - 1) * Dim.Z + 1 : end) = Hessian(:, (nStages - 1) * Dim.Z + 1 : end) + full(L_T_hessian);
        
    case 'FRP'
        %% FRP iteration routine
        FRP_L_hessian = FunObj.FRP_L_hessian(Var.x, Var.u, Var.p, Var.w, FRP.ZRef, FRP.ZWeight);
        Hessian = full(FRP_L_hessian);
end

end

