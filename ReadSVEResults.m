function [vfs, Cs, cntr, fnout, subfolderName] = ReadSVEResults(BC, SVESizeInv, noVE)
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

folder = ['csv_results/homo_data_', BC];
SVEsn = num2str(SVESizeInv);
subfolderName = ['Res', SVEsn, 'x', SVEsn];
fnBase = [folder, '/', subfolderName, '/img']; 
cntr = 0;
Cs = cell(0);
vfs = [];
for I = 1:noVE
    filename = [fnBase, num2str(I), '.csv'];
    A = readmatrix(filename);
    [m, n] = size(A);
    for i = 1:m
        cntr = cntr + 1;
        Cs{cntr} = A(i, 6:14);
        vfs(cntr) = A(i, 5);
    end
end
fnout = [subfolderName, '_', BC];