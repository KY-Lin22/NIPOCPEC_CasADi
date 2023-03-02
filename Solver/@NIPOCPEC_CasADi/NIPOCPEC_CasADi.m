classdef NIPOCPEC_CasADi < handle
    %Implementation of our proposed non-interior-point algorithm to solve 
    %the optimal control problem with equilibrium constraint:
    % min L_T(x(t), p(t)) + \int_0^T L_S(x(t), u(t), p(t))
    %s.t. G(x(t), u(t), p(t)) >= 0
    %     C(x(t), u(t), p(t)) = 0
    %     \dot(x) = f(x(t), u(t), p(t))
    %     p(t) \in SOL([lbp, ubp], K(x(t), u(t), p(t)))
    
    properties
        OCPEC % struct, problem data and symbolic representation of discretized(and reformulated) OCPEC, with field:
              %        'timeStep', 'nStages', 'InitState',        
              %        'Dim'(struct, variable dimension),                       
              %        'x', 'u', 'p', 'w'(SX symbolic variable, auxilary variable for function K),     
              %        'L_T', 'L_S',
              %        'G', 'C', 
              %        'f', 'xPrev'(SX symbolic variable, previous state), 'F'(discretization state equation),
              %        'K', 'lbp', 'ubp', 's'(SX symbolic variable, perturbed parameter for BVI), 
              %        'PHI'(reformulate BVI as inequalities PHI >= 0 using Scholtes regularization method),                 
        Option % struct, solver option
        
        FunObj % struct, CasADi function object
                % with field 'f', 'K',
                %            'L_T', 'L_S', 'G', 'C', 'F', 'PHI',
                %            'f_grad', 
                %            'L_T_grad', 'L_S_grad', 'G_grad', 'C_grad', 'F_grad', 'PHI_grad',
                %            'L_T_hessian', 'L_S_hessian', 'HAM_hessian', 
                %            'FB_G', 'FB_G_grad', 'FB_PHI','FB_PHI_grad',
                %            'G_Fsb', 'PHI_Fsb', 'HAM_grad'
                %            'J'
                %            'FRP_L', 'FRP_L_grad', 'FRP_L_hessian'                
    end
    %% Constructor Method for NIPOCPEC_CasADi     
    methods
        function self = NIPOCPEC_CasADi(OCPEC)
            %NIPOCPEC_CasADi: Construct an instance of this class
            % Syntax:
            %          self = NIPOCPEC_CasADi(OCPEC)
            % Argument:
            %          OCPEC: struct, containing OCPEC problem data and symbolic representation, with field: 
            %                'timeStep' -- double, discretization time step
            %                'nStages'  -- double, number of discretized stages
            %                'InitState'-- double, initial state            
            %                'x'   -- SX symbolic variable, state variable
            %                'u'   -- SX symbolic variable, control variable
            %                'p'   -- SX symbolic variable, algebraic variable
            %                'L_T' -- SX symbolic variable, terminal cost function L_T(x,p)
            %                'L_S' -- SX symbolic variable, stage cost function L_S(x,u,p)
            %                'G'   -- SX symbolic variable, inequality constraint G(x,u,p) >= 0
            %                'C'   -- SX symbolic variable, equality constraint C(x,u,p) = 0
            %                'f'   -- SX symbolic variable, state equation \dot(x) = f(x,u,p) 
            %                'K'   -- SX symbolic variable, function used to define Box constraint Variational Inequality(BVI) for p
            %                'lbp' -- double, lower bounds used to define Box constraint Variational Inequality(BVI) for p
            %                'ubp' -- double, upper bounds used to define Box constraint Variational Inequality(BVI) for p          
            % Output:
            %          self: instance of this class
            % import CasADi to workspace
            import casadi.*
            %% check and poblish input
            OCPEC = self.checkInput(OCPEC);
            
            %% initialize properties: OCPEC (reformulate and discretize OCPEC)
            % problem data
            Dim = struct('x', size(OCPEC.x, 1), 'u', size(OCPEC.u, 1), 'p', size(OCPEC.p, 1), 'w', size(OCPEC.p, 1));
            self.OCPEC = struct('timeStep', OCPEC.timeStep, 'nStages', OCPEC.nStages, 'InitState', OCPEC.InitState,...
                'Dim', Dim);     
            
            % symbolic representation of problem primal variable
            self.OCPEC.x = OCPEC.x;
            self.OCPEC.u = OCPEC.u;
            self.OCPEC.p = OCPEC.p;
            w = SX.sym('w', Dim.w, 1);
            self.OCPEC.w = w;
            
            % symbolic representation of problem function: L_T, L_S
            self.OCPEC.L_T = OCPEC.L_T;

            pWeight = 0.001 * eye(Dim.p);
            wWeight = 0.001 * eye(Dim.w);                    
            L_S = OCPEC.timeStep * (OCPEC.L_S ...
                + 0.5 * OCPEC.p' * pWeight * OCPEC.p ...
                + 0.5 * w' * wWeight * w); 
            self.OCPEC.L_S = L_S;
            
            % symbolic representation of problem function: G, C
            self.OCPEC.G = OCPEC.G;
            
            C = [OCPEC.C;...
                w - OCPEC.K];            
            self.OCPEC.C = C;
            
            % problem reformulation: f and F            
            self.OCPEC.f = OCPEC.f; 
            
            xPrev = SX.sym('xPrev', Dim.x, 1);
            self.OCPEC.xPrev = xPrev;
            F = xPrev - OCPEC.x + OCPEC.f * OCPEC.timeStep;
            self.OCPEC.F = F;
            
             % problem reformulation: PHI
            self.OCPEC.K   = OCPEC.K;
            self.OCPEC.lbp = OCPEC.lbp;
            self.OCPEC.ubp = OCPEC.ubp;
            s = SX.sym('s', 1, 1);
            self.OCPEC.s = s;
            
            NCP_num = 0;
            for i = 1 : Dim.p
                if (OCPEC.lbp(i) == 0) && (OCPEC.ubp(i) == Inf)
                    NCP_num = NCP_num + 1;
                end
            end
            BVI_num = Dim.p - NCP_num;            
            PHI = SX.sym('PHI', 3 * NCP_num + 4 * BVI_num, 1);
            PHI_Counter = 0;
            for i = 1 : Dim.p
                if (OCPEC.lbp(i) == 0) && (OCPEC.ubp(i) == Inf)
                    % nonlinear complementary problem
                    PHI(PHI_Counter + 1 : PHI_Counter + 3, 1) = ...
                        [OCPEC.p(i);...
                        w(i);...
                        s - OCPEC.p(i) * w(i)];                    
                    PHI_Counter = PHI_Counter + 3;
                else
                    % box constraint variation inequality
                    PHI(PHI_Counter + 1 : PHI_Counter + 4, 1) = ...
                        [OCPEC.p(i) - OCPEC.lbp(i);...
                        OCPEC.ubp(i) - OCPEC.p(i);...
                        s - (OCPEC.p(i) - OCPEC.lbp(i)) * w(i);...
                        s + (OCPEC.ubp(i) - OCPEC.p(i)) * w(i)];              
                    PHI_Counter = PHI_Counter + 4;
                end
            end            
            self.OCPEC.PHI = PHI;
            
            % dim for dual variable
            self.OCPEC.Dim.sigma  = size(self.OCPEC.G, 1);
            self.OCPEC.Dim.eta    = size(self.OCPEC.C, 1);
            self.OCPEC.Dim.lambda = self.OCPEC.Dim.x;
            self.OCPEC.Dim.gamma  = size(self.OCPEC.PHI, 1);
            
            % dim sum and node
            self.OCPEC.Dim.Z = self.OCPEC.Dim.x + self.OCPEC.Dim.u + self.OCPEC.Dim.p + self.OCPEC.Dim.w;
            self.OCPEC.Dim.LAMBDA = self.OCPEC.Dim.sigma + self.OCPEC.Dim.eta + self.OCPEC.Dim.lambda + self.OCPEC.Dim.gamma;
            self.OCPEC.Dim.Y = self.OCPEC.Dim.Z + self.OCPEC.Dim.LAMBDA;
            self.OCPEC.Dim.Node = cumsum([self.OCPEC.Dim.sigma, self.OCPEC.Dim.eta, self.OCPEC.Dim.lambda, self.OCPEC.Dim.gamma,...
                self.OCPEC.Dim.x, self.OCPEC.Dim.u, self.OCPEC.Dim.p, self.OCPEC.Dim.w]);                       
            
            %% initialize properties: Option
            self.Option = self.createOption();            
            
            %% initialize properties: FunObj
            self.FunObj = self.createFunObj();           
            
        end

    end
    
    %% Other Methods for NIPOCPEC_CasADi    
    methods
        %%
        OCPEC = checkInput(self, OCPEC)        
        
        Option = createOption(self)
        
        FunObj = createFunObj(self)
        
        showInfo(self)
        
        generateInitialGuess(self)
        
        [solution, Info] = solveOCPEC(self, Var_Init)
        
        showResult(self, Info)
        
        %% Methods in solveOCPEC method
        % function, jacobian and Hessian evaluation
        Fun = FunctionEvaluation(self, Var, s, z, mode, FRP)        
        
        Jac = JacobianEvaluation(self, Var, Fun, s, z, mode, FRP)
        
        Hessian = HessianEvaluation(self, Var, Jac, s, mode, FRP)
        
        % KKT evaluation
        [KKT_Residual, KKT_Error] = computeKKT_Residual_Error(self, Var, Fun, Jac)        
        
        KKT_Matrix = computeKKT_Matrix(self, Jac, Hessian)
        
        % evaluate search direction
        [dY, Info] = SearchDirection_Riccati(self, KKT_Residual, KKT_Matrix)    
        
        % Line Search  
        [Var_LS, Info] = LineSearch_Merit(self, Var, Fun, Jac, beta, s, z,...
            KKT_Residual, KKT_Matrix, dY_k, mode, FRP)        
        
        % Second Order Correction
        Var_SOC = SecondOrderCorrection(self, Var, Jac, KKT_Residual, KKT_Matrix, Fun_full)        
        
        % Feasibility Restoration Phase
        [Var_FRP, Info] = FeasibilityRestorationPhase_MinDeviation(self, Var_Ref, Fun_Ref, Jac_Ref, s, z)
        
        [eta, lambda] = computeDualVar_lsqminnorm(self, Var, s)
        
         % compute perturbed parameter s and z
        [s_k, z_k, Fun_k] = computePerturedParam(self, Var_k, Fun_k, s, z)
        
        % examine solution
        Info = solutionExaminer(self, solution, Record)          
    end
     
end

