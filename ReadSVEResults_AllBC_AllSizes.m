function [vfs, Cs, KR2Ds, AROs, cntrs, fnouts, BCnames, subfolderName] = ReadSVEResults_AllBC_AllSizes(SVESizeInvs, noVE)

if (nargin < 1)
    SVESizeInvs = [16, 8, 4, 2, 1];
end

if (nargin < 2)
    noVE = 30;
end

nSVEsz = length(SVESizeInvs);
for SVEszi = 1:nSVEsz
    [vfs{SVEszi}, Cs{SVEszi}, KR2Ds{SVEszi}, AROs{SVEszi}, cntrs{SVEszi}, fnouts{SVEszi}, BCnames, subfolderName{SVEszi}] = ReadSVEResults_AllBC_OneSize(SVESizeInvs(SVEszi), noVE);
end

