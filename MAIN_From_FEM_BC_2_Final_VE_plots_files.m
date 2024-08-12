SVESizeInvs = [16, 8, 4, 2, 1];
noVE = 30;
epsilon = 0.01;
use_minMax_over_meanpmstd = 1;
dim = 2;
logbase = 10;
cI = 0.23905;

% reading inputs
if 0
fprintf(1, 'ReadSVEResults_AllBC_AllSizes ...\n');
ReadSVEResults_AllBC_AllSizes(SVESizeInvs, noVE);
fprintf(1, 'CreateAllStatFiles ...\n');
CreateAllStatFiles(SVESizeInvs);
fprintf(1, 'CreateAllFieldStats4_plots_RVESizeAnalysis ...\n');
CreateAllFieldStats4_plots_RVESizeAnalysis(SVESizeInvs);
end
fprintf(1, 'GenerateSVESizeEffectData_AllFields_AllBCc ...\n');
GenerateSVESizeEffectData_AllFields_AllBCc(epsilon, use_minMax_over_meanpmstd, dim, logbase, cI);

