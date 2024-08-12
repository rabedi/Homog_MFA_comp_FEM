function [s1, s2] = Kmu2_s1s2(K, mu, mode)
% s1 and s2 are bulk and shear components of Eshelby tensor

% mode 3: 3D
% mode 21: plane strain
% mode 22: plane stress
if nargin  < 3
    mode = 21;
end
if (mode == 3) % 3D
    s1 = K / (K + 4.0 / 3.0 * mu);
    s2 = (6.0 * (K + 2.0 * mu)) / (5.0 * (3.0 * K + 4.0 * mu));
else % 2D cases are the same
    s1 = K / (K + mu);
    s2 = (K + 2.0 * mu) / (2.0 * (K + mu));
end