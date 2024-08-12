function [alphas, As, rsqs, lVs, lVar2s, lVar2_fits, nExistings, nRequireds, middle_low_highs] = ...
GenerateSVESizeEffectData_OneField_AllBCc(preName, fieldName, BCnames, epsilon, use_cov_over_std, use_minMax_over_meanpmstd, dim, logbase, cI)

if nargin < 1
    preName = 'C';
end
if nargin < 2
    fieldName = 'K';
end
if nargin < 3
    BCnames = {'disp', 'mixed', 'trac', 'MT', 'SC', 'Diff'};
end
if nargin < 4
    epsilon = 0.01;
end
if (nargin < 5)
    % to investigate convergence of a property to variation based size, we use
    % Coefficient of variation for dimensional parameters, e.g. K
    % use std for nondimensional ones (e.g. Anisotropy indices, Poisson
    % ratio)
    use_cov_over_std = 1;
end
if (nargin < 6)
    % to generate mean, min/max (value below 1) or mean, mean-/+ std (value
    % below 0)
    use_minMax_over_meanpmstd = 1;
end
if (nargin < 7)
    % domain dimension
    dim = 2;
end
if (nargin < 8)
    logbase = 10;
end
if (nargin < 9)
    cI = 0.23905;
end

numBC = length(BCnames);
hasBC = ((numBC > 1) || (strcmp(BCnames{1}, 'none') == 0));
BCname_leg = {}; 
hasLeg = 0;
if (hasBC)
    BCprename = 'BC_';
    for BCi = 1:numBC
        BCname_leg{BCi} = BCnames{BCi};
        nameWOExt{BCi} = [preName, '_', fieldName, '_BC_', BCnames{BCi}];
    end
else
    nameWOExt{1} = [preName, '_', fieldName];
end

for BCi = 1:numBC
    [alphas(BCi), As(BCi), rsqs(BCi), RVE_edges(BCi), RVE_vols(BCi), lVs{BCi}, lVar2s{BCi}, lVar2_fits{BCi}, nExistings{BCi}, nRequireds{BCi}, middle_low_highs{BCi}] = ...
    GenerateSVESizeEffectData_OneField_OneBC(nameWOExt{BCi}, epsilon, use_cov_over_std, use_minMax_over_meanpmstd, dim, logbase);
end

fileNameBase0 = [preName, '_', fieldName];
fileName = [fileNameBase0, '_VariationBasedRVESizeStat.txt'];
fid = fopen(fileName, 'w');
fprintf(fid, 'BC\talpha\tA\tR2\tRVE_edge\tRVE_vol\n');
for BCi = 1:numBC
    fprintf(fid, '%s\t%g\t%g\t%g\t%g\t%g\n', BCnames{BCi}, alphas(BCi), As(BCi), rsqs(BCi), RVE_edges(BCi), RVE_vols(BCi));
end
fclose(fid);

for BCi = 1:numBC
    fileName = [fileNameBase0, '_BC_', BCnames{BCi}, '_EachSizeDat.txt'];
    fid = fopen(fileName, 'w');
    lV = lVs{BCi};
    lEdge = 1.0 / dim * lV;
    edgeInv = power(logbase, -lEdge);
    lVar2 = lVar2s{BCi};
    lVar2_fit = lVar2_fits{BCi};
    nExisting = nExistings{BCi};
    nRequired = nRequireds{BCi};
    middle = middle_low_highs{BCi}{1};
    low = middle_low_highs{BCi}{2};
    high = middle_low_highs{BCi}{3};
    sz = length(lV);
    fprintf(fid, 'sizeIndex\tlogEdge\tedgeInv\tlogVol\tlogVar2\tlogVar2_fit\tExisting\tnRequired\tmeanVal\tlowVal\thighVal\n');
    for j = 1:sz
        fprintf(fid, '%d\t', j);
        fprintf(fid, '%g\t', lEdge(j));
        fprintf(fid, '%g\t', edgeInv(j));
        fprintf(fid, '%g\t', lV(j));
        fprintf(fid, '%g\t', lVar2(j));
        fprintf(fid, '%g\t', lVar2_fit(j));
        fprintf(fid, '%g\t', nExisting(j));
        fprintf(fid, '%g\t', nRequired(j));
        fprintf(fid, '%g\t', middle(j));
        fprintf(fid, '%g\t', low(j));
        fprintf(fid, '%g\n', high(j));
    end
    fclose(fid);
end

% plotting variation based results
labsz = 25; % x, y label font size
lfs = 13;
xlab = '$$ \mathrm{log}(V_{\mathrm{SVE})} $$';
K_mu_E_nu_lambda_index = -1;
if (strcmp(fieldName, 'mu') == 1)
    ylab = ['\', fieldName];
    K_mu_E_nu_lambda_index = 1;
elseif (strcmp(fieldName, 'nu') == 1)
    ylab = ['\', fieldName];
    K_mu_E_nu_lambda_index = 4;
elseif  (strcmp(fieldName, 'lambda') == 1)
    ylab = ['\', fieldName];
    K_mu_E_nu_lambda_index = 5;
elseif (strcmp(fieldName, 'K') == 1)
    ylab = [fieldName];
    K_mu_E_nu_lambda_index = 1;
elseif (strcmp(fieldName, 'E') == 1)
    ylab = [fieldName];
    K_mu_E_nu_lambda_index = 3;
elseif (strcmp(fieldName, 'ANZ') == 1)
    ylab = 'A^{\mathrm{NZ}}';
elseif (strcmp(fieldName, 'ARO') == 1)
    ylab = 'A^{\mathrm{RO}}';
elseif (strcmp(fieldName, 'vf') == 1)
    ylab = 'v_f';
else
    fprintf(1, 'fieldName = %s -> ylab should be added\n', fieldName); 
end 
    
BCname_expanded_leg = BCname_leg;
ne = 0;
ye = cell(0);
lce = cell(0);
lse = cell(0);
if (K_mu_E_nu_lambda_index > 0)
    [res, nms, resVoigt, resResuss, resHSp, resHSm] = SiCB4(cI);
    valVoigt = resVoigt(K_mu_E_nu_lambda_index);
    valReuss = resResuss(K_mu_E_nu_lambda_index);
    valHSp = resHSp(K_mu_E_nu_lambda_index);
    valHSm = resHSm(K_mu_E_nu_lambda_index);
    fidvr = fopen('_summaryVoigtReus_HS.txt', 'w');
    fprintf(fidvr, 'model\tK\tmu\tE\tnu\tlambda\n');
    fprintf(fidvr, 'Voigt\t%g\t%g\t%g\t%g\t%g\t%g', resVoigt(1), resVoigt(2), resVoigt(3), resVoigt(4), resVoigt(5));
    fprintf(fidvr, '\n');
    fprintf(fidvr, 'Reuss\t%g\t%g\t%g\t%g\t%g\t%g', resResuss(1), resResuss(2), resResuss(3), resResuss(4), resResuss(5));
    fprintf(fidvr, '\n');
    fprintf(fidvr, 'HS+\t%g\t%g\t%g\t%g\t%g\t%g', resHSp(1), resHSp(2), resHSp(3), resHSp(4), resHSp(5));
    fprintf(fidvr, '\n');
    fprintf(fidvr, 'HS-\t%g\t%g\t%g\t%g\t%g\t%g', resHSm(1), resHSm(2), resHSm(3), resHSm(4), resHSm(5));
    fprintf(fidvr, '\n');
    fclose(fidvr);
    
    ne = 4;
    xe = lVs{1};
    ye{1} = xe * 0 + valVoigt;
    ye{2} = xe * 0 + valReuss;
    ye{3} = xe * 0 + valHSp;
    ye{4} = xe * 0 + valHSm;
    BCname_expanded_leg{numBC + 1} = 'Voigt';
    BCname_expanded_leg{numBC + 2} = 'Reuss';
    BCname_expanded_leg{numBC + 3} = 'HS+';
    BCname_expanded_leg{numBC + 4} = 'HS-';
    lce{1} = 'k';
    lce{2} = 'k';
    lce{3} = [0.7, 0.7, 0.7];
    lce{4} = [0.7, 0.7, 0.7];
    lse{1} = '--';
    lse{2} = '-.';
    lse{3} = '--';
    lse{4} = '-.';
end

lc = getColors(1, 0, 0);
fg = figure(1);
fileNameBase = ['plot_Variation_4_field_', fieldName];
set(fg,'defaultLegendAutoUpdate', 'off');
for BCi = 1:numBC
    x = lVs{BCi};
    y = middle_low_highs{BCi}{1};
    plot(x, y, 'Color', lc{BCi}, 'LineStyle', '-', 'LineWidth', 2);
    hold on;
end

for ie = 1:ne
    plot(xe , ye{ie}, 'Color', lce{ie}, 'LineStyle', lse{ie}, 'LineWidth', 2);
    hold on;
end    
if (length(BCname_expanded_leg) > 1)
    lg = legend(BCname_expanded_leg, 'FontSize', lfs, 'Interpreter', 'latex');
    legend('boxoff');
end
for BCi = 1:numBC
    x = lVs{BCi};
    y = middle_low_highs{BCi}{2};
    plot(x, y, 'Color', lc{BCi}, 'LineStyle', '--', 'LineWidth', 2);
    hold on;
end
for BCi = 1:numBC
    x = lVs{BCi};
    y = middle_low_highs{BCi}{3};
    plot(x, y, 'Color', lc{BCi}, 'LineStyle', '-.', 'LineWidth', 2);
    hold on;
end
xh = get(gca, 'XLabel');
set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
set(yh, 'String', ['$$ ', ylab, ' $$'], 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

fileNameBase = [fileNameBase0, '_mean_lh_vs_vol'];
print('-dpng', [fileNameBase, '.png']);
savefig([fileNameBase, '.fig']);

% plotting RVE size
fg = figure(2);
fileNameBase = ['plot_RVE_size_4_field_', fieldName];

set(fg,'defaultLegendAutoUpdate', 'off');
for BCi = 1:numBC
    x = lVs{BCi};
    y = lVar2s{BCi};
    plot(x, y, 'Color', lc{BCi}, 'LineStyle', '-', 'LineWidth', 2);
    hold on;
end
if (length(BCname_leg) > 1)
    lg = legend(BCname_leg, 'FontSize', lfs, 'Interpreter', 'latex');
    legend('boxoff');
end
for BCi = 1:numBC
    x = lVs{BCi};
    y = lVar2_fits{BCi};
    plot(x, y, 'Color', lc{BCi}, 'LineStyle', '--', 'LineWidth', 2);
    hold on;
end
y = x * 0 + 2 * log(epsilon) / log(logbase);

plot(x, y, 'Color', [0.7, 0.7, 0.7], 'LineStyle', '-', 'LineWidth', 1);
hold on;

y = x * 0 + 2.2 * log(epsilon) / log(logbase);
plot(x, y, 'Color', 'none', 'LineStyle', '-', 'LineWidth', 1);
hold on;
y = x * 0 + 1.8 * log(epsilon) / log(logbase);
plot(x, y, 'Color', 'none', 'LineStyle', '-', 'LineWidth', 1);
hold on;

xh = get(gca, 'XLabel');
set(xh, 'String', xlab, 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
yh = get(gca, 'YLabel');
ad = '\mathrm{COV}^2';
if (use_cov_over_std == 0)
    ad = 'VAR';
end
ylabc = ['$$ \mathrm{log}_{10}(', ad, '(', ylab, ' ))$$'];
set(yh, 'String', ylabc, 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');

fileNameBase = [fileNameBase0, '_RVE_size_var_based'];
print('-dpng', [fileNameBase, '.png']);
savefig([fileNameBase, '.fig']);

close('all');
fclose('all');

