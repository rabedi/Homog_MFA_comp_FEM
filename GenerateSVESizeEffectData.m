function GenerateSVESizeEffectData(nameWOExt, epsilon, use_cov_over_std, use_minMax_over_meanpmstd, dim, logbase)
if nargin < 1
    nameWOExt = 'C_K_BC_Diff';
end

if nargin < 2
    epsilon = 0.01;
end

if (nargin < 3)
    % to investigate convergence of a property to variation based size, we use
    % Coefficient of variation for dimensional parameters, e.g. K
    % use std for nondimensional ones (e.g. Anisotropy indices, Poisson
    % ratio)
    use_cov_over_std = 1;
end

if (nargin < 4)
    % to generate mean, min/max (value below 1) or mean, mean-/+ std (value
    % below 0)
    use_minMax_over_meanpmstd = 1;
end

if (nargin < 5)
    % domain dimension
    dim = 2;
end

if (nargin < 6)
    logbase = 10;
end

% columns
% SVESizeInv	n	mean	std	coef	min	max	mean-std	mean+std
filename = [nameWOExt, '.txt'];
A = readmatrix(filename);

lengths = A(:,1);
nSVEsz = length(lengths);
ovlength =  1 ./ lengths;
lovlenth = log(ovlength) / log(logbase);
middle_low_high{1} = A(:,3);
if (use_minMax_over_meanpmstd)
    middle_low_high{2} = A(:,6);
    middle_low_high{3} = A(:,7);
else
    middle_low_high{2} = A(:,8);
    middle_low_high{3} = A(:,9);
end
    
if (use_cov_over_std)
    varVals = A(:,5);
else
    varVals = A(:,4);
end
lvarVals = log(varVals) / log(logbase);

nExisting = A(:,2);
lV = dim * lovlenth;
lVar2 = 2 * lvarVals;
[alpha, A, lV, lVar2, ]
% variation based computation
p = polyfit(lV, lVar2, 1); % Fits a first-degree polynomial (line) to the data
% Generate points for the regression line
y_fit = polyval(p, lV); % Compute corresponding y values
slope = p(1);
intercept = p(2);
alpha = -slope; % / 2;
A = power(logbase, intercept);

% best fit
% lVar2 = A - alpha log(volume),
% lVar2 = epsilon^2 for RVE
% -> 
%log eps^2 = A - p log(RVE_vol) ->
% RVE_vol = logbase^(A - 2 * log(eps) / alpha)
% RVE_edge = sqrt(RVE_vol)
logepsilon = log(epsilon) / log(logbase);
RVE_vol = power(logbase, (intercept - 2.0 * logepsilon) / alpha);
RVE_edge = power(RVE_vol, 1.0 / dim);

% nRequired
% Var (or COV^2) / n = epsilon^2, use the linear regression line for this
% std / sqrt(n) = epsilon -> Var / n = epsilon^2, take the log ->
% log(Var) - log(n) = 2 * log10(epsilon)
% log(Var) is the best fit for any given size ->

% the same argument works for quantities wherein std rather than cov is
% used (e.g. vf and Anisotropy index)
tlogeps = 2.0 * logepsilon;
for szi = 1:nSVEsz
    logVar = y_fit(szi);
    logn = logVar - tlogeps;
    nRequired(szi) = ceil(power(logbase, logn));
end

