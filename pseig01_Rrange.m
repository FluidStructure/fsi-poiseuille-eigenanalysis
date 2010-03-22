Rrange = linspace(5000,5000,1)

# Source the file with all the octave functions
source('octaveFunctions.m')

% Set up the variables and influence coefficients for a nominal case
setupMatrices

for R = Rrange
    disp(['CALCULATING EIGENVALUES FOR R = ' num2str(R)])

    Rt=R;load VARS.mat;R=Rt;clear Rt
    R
    mu = Uinf*rho*delta/R;      % Flow viscosity
    nu = mu./rho;               % Flow kinematic viscosity
    save -v7 VARS.mat

    % Recalculate the matrices
    evaluateMatrices
    % Recalculate the eigenvalues
    calculateEigenmodes
end

