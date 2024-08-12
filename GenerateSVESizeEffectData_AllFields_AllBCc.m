function GenerateSVESizeEffectData_AllFields_AllBCc(epsilon, use_minMax_over_meanpmstd, dim, logbase, cI)

BCnames = {'disp', 'mixed', 'trac', 'MT', 'SC', 'Diff'};
if nargin < 1
    epsilon = 0.01;
end
if (nargin < 2)
    % to generate mean, min/max (value below 1) or mean, mean-/+ std (value
    % below 0)
    use_minMax_over_meanpmstd = 1;
end
if (nargin < 3)
    % domain dimension
    dim = 2;
end
if (nargin < 4)
    logbase = 10;
end
if (nargin < 5)
    cI = 0.23905;
end


%%% K, mu, E, nu, lambda
preName = 'C';
fieldNames = {'K', 'mu', 'E', 'nu', 'lambda'};
BCnames = {'disp', 'mixed', 'trac', 'MT', 'SC', 'Diff'};
nf = length(fieldNames);
use_cov_over_stds = ones(nf, 1);
use_cov_over_stds(4) = 0;

for fi = 1:nf
    [alphas, As, rsqs, lVs, lVar2s, lVar2_fits, nExistings, nRequireds, middle_low_highs] = ...
    GenerateSVESizeEffectData_OneField_AllBCc(preName, fieldNames{fi}, BCnames, epsilon, use_cov_over_stds(fi), use_minMax_over_meanpmstd, dim, logbase, cI);
end


%%% Anisotropy indices
preName = 'As';
fieldNames = {'ANZ', 'ARO'};
BCnames = {'disp', 'mixed', 'trac'};
nf = length(fieldNames);
use_cov_over_stds = zeros(nf, 1);

for fi = 1:nf
    [alphas, As, rsqs, lVs, lVar2s, lVar2_fits, nExistings, nRequireds, middle_low_highs] = ...
    GenerateSVESizeEffectData_OneField_AllBCc(preName, fieldNames{fi}, BCnames, epsilon, use_cov_over_stds(fi), use_minMax_over_meanpmstd, dim, logbase, cI);
end


%%%% other scalars
preName = 'scalars';
fieldNames = {'vf'};
BCnames = {'none'};
nf = length(fieldNames);
use_cov_over_stds = zeros(nf, 1);

for fi = 1:nf
    [alphas, As, rsqs, lVs, lVar2s, lVar2_fits, nExistings, nRequireds, middle_low_highs] = ...
    GenerateSVESizeEffectData_OneField_AllBCc(preName, fieldNames{fi}, BCnames, epsilon, use_cov_over_stds(fi), use_minMax_over_meanpmstd, dim, logbase, cI);
end