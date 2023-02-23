function OCPEC = checkInput(self, OCPEC)
%check and polish the input of NIPOCPEC_CasADi: OCPEC struct
%   
%% check Input: OCPEC
% check timeStep
if isempty(OCPEC.timeStep)
    error('please specify the discretization time step in OCPEC.timeStep')
end

% check nStages
if isempty(OCPEC.nStages)
    error('please specify the number of discretization stages in OCPEC.nStages')
end

% check InitState
if isempty(OCPEC.InitState)
    error('please specify initial state in OCPEC.InitState')
else
    if ~all(size(OCPEC.InitState) == size(OCPEC.x))
        error('OCPEC.InitState should have the same dimension as OCPEC.x')
    end
end

% check x
if isempty(OCPEC.x)
    error('please specify state variable in OCPEC.x')
else
    if size(OCPEC.x, 2) ~= 1
        error('OCPEC.x should be a column vector')
    end
end

% check u
if isempty(OCPEC.u)
    OCPEC.u = casadi.SX.sym('u', 0, 1);
else
    if size(OCPEC.u, 2) ~= 1
        error('OCPEC.u should be a column vector')
    end
end

% check p
if isempty(OCPEC.p)
    OCPEC.p = casadi.SX.sym('p', 0, 1);
else
    if size(OCPEC.p, 2) ~= 1
        error('OCPEC.p should be a column vector')
    end
end

% check L_T
if isempty(OCPEC.L_T)
    OCPEC.L_T = casadi.SX.sym('L_T', 0, 1);
else
    if ~all(size(OCPEC.L_T) == [1, 1])
        error('OCPEC.L_T should be a scale function')
    end
end

% check L_S
if isempty(OCPEC.L_S)
    error('please specify stage cost function in OCPEC.L_S')
else
    if ~all(size(OCPEC.L_S) == [1, 1])
        error('OCPEC.L_S should be a scale function')
    end
end

% check G
if isempty(OCPEC.G)
    OCPEC.G = casadi.SX.sym('G', 0, 1);
else
    if size(OCPEC.G, 2) ~= 1
        error('OCPEC.G should be a column function')
    end
end

% check C
if isempty(OCPEC.C)
    OCPEC.C = casadi.SX.sym('C', 0, 1);
else
    if size(OCPEC.C, 2) ~= 1
        error('OCPEC.C should be a column function')
    end
end

% check f
if isempty(OCPEC.f)
    error('please specify the state equation in OCPEC.f')
else
    if ~all(size(OCPEC.f) == size(OCPEC.x))
        error('OCPEC.f should have the same dimension as OCPEC.x')
    end
end

% check K
if isempty(OCPEC.K)
    if size(OCPEC.p, 1) == 0
        OCPEC.K = casadi.SX.sym('K', 0, 1);
    else
        error('please specify the function K in OCPEC.K')
    end
else
    if ~all(size(OCPEC.K) == size(OCPEC.p))
        error('OCPEC.K should have the same dimension as OCPEC.p')
    end
end

% check lbp
if isempty(OCPEC.lbp)
    if size(OCPEC.p, 1) == 0
        OCPEC.lbp = zeros(0, 1);
    else
        error('please specify the lower bounds in OCPEC.lbp')
    end   
else   
    OCPEC.lbp = double(OCPEC.lbp);
    if ~all(size(OCPEC.lbp) == size(OCPEC.p))
        error('OCPEC.lbp should have the same dimension as OCPEC.p')
    end
end

% check ubp
if isempty(OCPEC.ubp)
    if size(OCPEC.p, 1) == 0
        OCPEC.ubp = zeros(0, 1);
    else
        error('please specify the upper bounds in OCPEC.ubp')
    end   
else    
    OCPEC.ubp = double(OCPEC.ubp);
    if ~all(size(OCPEC.ubp) == size(OCPEC.p))
        error('OCPEC.ubp should have the same dimension as OCPEC.p')
    end
end

% check lbp <= ubp
if ~all(OCPEC.lbp <= OCPEC.ubp)
    error('OCPEC.lbp should less than or equal to OCPEC.ubp')
end
end

