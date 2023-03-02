function [s_k, z_k, Fun_k] = computePerturedParam(self, Var_k, Fun_k, s, z)
%UNTITLED23 Summary of this function goes here
%   Detailed explanation goes here
OCPEC = self.OCPEC;
Option = self.Option;
FunObj = self.FunObj;

nStages = OCPEC.nStages;

Tol = Option.Tolerance;
sEnd = Option.sEnd;
zEnd = Option.zEnd;

kappa_F = Option.kappa_F;
kappa_s_times = Option.kappa_s_times;
kappa_s_exp = Option.kappa_s_exp;
kappa_z_times = Option.kappa_z_times;
kappa_z_exp = Option.kappa_z_exp;

% threshold 
threshold_sz = max([s, z, Tol.KKT_Error_Total]);

% constraint violation (L_Inf norm)
totalCstrVio_L_InfNorm = norm(reshape([Fun_k.PSIg; Fun_k.C; Fun_k.F; Fun_k.PSIphi], [], 1), Inf); 

%% update s and z
if totalCstrVio_L_InfNorm < kappa_F * threshold_sz
    % close to the optimal solution, need to set them a smaller value
    s_k_trail = min([kappa_s_times .* s, s.^kappa_s_exp]);
    s_k = max([s_k_trail, sEnd]);
    z_k_trail = min([kappa_z_times .* z, z.^kappa_z_exp]);
    z_k = max([z_k_trail, zEnd]);
else
    % far way from the optimal solution, set them the previous value
    s_k = s;
    z_k = z;
end

%% update Fun_k
s_k_Repmat = repmat(s_k, 1, nStages);
z_k_Repmat = repmat(z_k, 1, nStages);

if (s_k == s) && (z_k == z)
    % both s and z do not update, hence no function need to be updated
    
elseif (s_k ~= s) && (z_k == z)
    % only s is update, hence update function about PHI
    % PHI
    PHI_k = FunObj.PHI(Var_k.x, Var_k.u, Var_k.p, Var_k.w, s_k_Repmat);
    Fun_k.PHI = full(PHI_k);
    % FB for PHI
    PSIphi_k = FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k_Repmat);
    Fun_k.PSIphi = full(PSIphi_k);      
       
elseif (s_k == s) && (z_k ~= z)
    % only z is update, hence update function about FB
    % FB for G
    PSIg_k = FunObj.FB_G(Var_k.sigma, Fun_k.G, z_k_Repmat);
    Fun_k.PSIg = full(PSIg_k);    
    % FB for PHI
    PSIphi_k = FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k_Repmat);
    Fun_k.PSIphi = full(PSIphi_k);       
else
    % both s and z update, hence update function about PHI and FB
    % PHI
    PHI_k = FunObj.PHI(Var_k.x, Var_k.u, Var_k.p, Var_k.w, s_k_Repmat);
    Fun_k.PHI = full(PHI_k);
    % FB for G
    PSIg_k = FunObj.FB_G(Var_k.sigma, Fun_k.G, z_k_Repmat);
    Fun_k.PSIg = full(PSIg_k);    
    % FB for PHI
    PSIphi_k = FunObj.FB_PHI(Var_k.gamma, Fun_k.PHI, z_k_Repmat);
    Fun_k.PSIphi = full(PSIphi_k);             
end

end

