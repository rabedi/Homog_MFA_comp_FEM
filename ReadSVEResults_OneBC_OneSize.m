function [vfs, Cs, KR2Ds, AROs, cntr, fnout, subfolderName] = ReadSVEResults_OneBC_OneSize(BC, SVESizeInv, noVE)

% BC disp, trac, mixed
% SVESizeInv 1, 2, 4, 8, 16
% noVE = 30
if (nargin < 1)
    BC = 'disp';
end
if (nargin < 2)
    SVESizeInv = 2;
end
if (nargin < 3)
    noVE = 30;
end

[vfs, Cs, cntr, fnout, subfolderName] = ReadSVEResults(BC, SVESizeInv, noVE);
KR2Ds = zeros(cntr, 5);
AROs = zeros(cntr, 2);
for i = 1:cntr
    [tmpKmu, tmpAs, C, CIso] = Compute_k_A_from_C(Cs{i});    
    for j = 1:5
        KR2Ds(i, j) = tmpKmu(j);
    end
    for j = 1:2
        AROs(i, j) = tmpAs(j);
    end
end

fid = fopen([fnout, '_KMuENuLambda.txt'], 'w');
for i = 1:cntr
    tmpV = KR2Ds(i, :);
    fprintf(fid, '%g\t%g\t%g\t%g\t%g\n', tmpV(1), tmpV(2), tmpV(3), tmpV(4), tmpV(5));
end

fid = fopen([fnout, '_As.txt'], 'w');
for i = 1:cntr
    tmpV = AROs(i, :);
    fprintf(fid, '%g\t%g\n', tmpV(1), tmpV(2));
end

fid = fopen([fnout, '_Cs.txt'], 'w');
for i = 1:cntr
    C = Cs{i};
    fprintf(fid, '%g', C(1));
    for j = 2:9
        fprintf(fid, '\t%g', C(j));
    end
    fprintf(fid, '\n');
end

fid = fopen([subfolderName, '_vf.txt'], 'w');
for i = 1:cntr
    fprintf(fid, '%g\n', vfs(i));
end
fclose(fid);

