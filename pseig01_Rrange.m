Rrange = linspace(5000,5000,1)

% Set up the variables and influence coefficients for a nominal case
pseig01_a

for R = Rrange
    disp(['CALCULATING EIGENVALUES FOR R = ' num2str(R)])

    Rt=R;load VARS.mat;R=Rt;clear Rt
    R
    mu = Uinf*rho*delta/R;      % Flow viscosity
    nu = mu./rho;               % Flow kinematic viscosity
    save -v7 VARS.mat

    % Recalculate the matrices
    pseig01_b
    % Recalculate the eigenvalues
    pseig01_c
end

