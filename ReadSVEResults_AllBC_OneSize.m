function [vfs, Cs, KR2Ds, AROs, cntrs, fnouts, BCnames, subfolderName] = ReadSVEResults_AllBC_OneSize(SVESizeInv, noVE)
% SVESizeInv 1, 2, 4, 8, 16
% noVE = 30
if (nargin < 1)
    SVESizeInv = 2;
end
if (nargin < 2)
    noVE = 30;
end

BCnamesA = {'disp', 'trac', 'mixed'};
% BC disp, trac, mixed
for BCi = 1:length(BCnamesA)
    BC = BCnamesA{BCi};
    [vfs, Cs{BCi}, KR2Ds{BCi}, AROs{BCi}, cntrs(BCi), fnouts{BCi}, subfolderName] = ...
        ReadSVEResults_OneBC_OneSize(BC, SVESizeInv, noVE);
end
num_vf = length(vfs);
fprintf(1, '\tCalculating MFA SVESizeInv = %d\t ...', SVESizeInv);
for i = 1:num_vf
    [res, nms, resVoigt, resResuss, resHSp, resHSm] = SiCB4(vfs(i));
    resi{i} = res;
end
fprintf(1, '\tfinished\n');
BCnamesB = nms;
szA = length(BCnamesA);
szB = length(BCnamesB);
szNames = szA + szB;
BCnames = BCnamesA;
for BCi = 1:szB
    BCii = BCi + szA;
    BCnames{BCii} = BCnamesB{BCi};

    tmpMat = zeros(num_vf, 5);
    for i = 1:num_vf
        tmpV = resi{i}{BCi};
        for j = 1:5
            tmpMat(i, j) = tmpV(j);
        end
    end
    KR2Ds{BCii} = tmpMat;
    fnout = [subfolderName, '_', BCnamesB{BCi}];
    fid = fopen([fnout, '_KMuENuLambda.txt'], 'w');
    for i = 1:num_vf
        tmpV = tmpMat(i, :);
        fprintf(fid, '%g\t%g\t%g\t%g\t%g\n', tmpV(1), tmpV(2), tmpV(3), tmpV(4), tmpV(5));
    end
end
