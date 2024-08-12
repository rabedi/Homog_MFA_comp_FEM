function [E, nu, lambda] = Kmu2_Enulambda(K, mu, mode)
% mode 3: 3D
% mode 21: plane strain
% mode 22: plane stress
if nargin  < 3
    mode = 21;
end
if (mode == 3) % 3D
    denomInv = 1.0 / (3.0 * K + mu);
    E = 9 * K * mu * denomInv;
    nu = 0.5 * (3.0 * K - 2.0 * mu) * denomInv;
    lambda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);
elseif (mode == 21) % plane strain
    nu = 0.5 * (K - mu) / K;
    E = mu * (3.0 * K - mu) / K;
    lambda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);
elseif (mode == 22) % plane stress
    nu = (K - mu) / (K + mu);
    E = 4.0 * K * mu / (K + mu);
    lambda = E * nu / (1.0 - nu * nu);
end