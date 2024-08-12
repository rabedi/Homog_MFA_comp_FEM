function [K, mu, lambda] = Enu2_Kmulambda(E, nu, mode)
% mode 3: 3D
% mode 21: plane strain
% mode 22: plane stress
if nargin  < 3
    mode = 21;
end
mu = E / 2.0 / (1.0 + nu);
if (mode == 3) % 3D
    K = E / 3.0 / (1.0 - 2.0 * nu);
    lambda = K - 2.0 * mu / 3.0;
elseif (mode == 21) % plane strain
    K = E / 2.0 / (1.0 + nu) / (1.0 - 2.0 * nu);
    lambda = K - 2.0 * mu / 3.0;
elseif (mode == 22) % plane stress
    K = E / 2.0 / (1.0 - nu);
    lambda = K - mu;
end