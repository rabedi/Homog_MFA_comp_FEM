function CreateAllFieldStats4_plots_RVESizeAnalysis(SVESizeInvs)

if (nargin < 1)
    SVESizeInvs = [16, 8, 4, 2, 1];
end

BCs = {'disp', 'mixed', 'trac', 'MT', 'SC', 'Diff'};
nBCs = length(BCs);


nSVEsz = length(SVESizeInvs);
names = {'n', 'mean', 'std', 'coefVar', 'min', 'max', 'mean-std', 'mean+std'};
nnames = length(names);

fieldNames = {'K', 'mu', 'E', 'nu', 'lambda'};
nFields = length(fieldNames);
%%%% stiffnessess

for BCi = 1:nBCs
    BC = BCs{BCi};
    for SVEszi = 1:nSVEsz
        SVEsz = SVESizeInvs(SVEszi);
        SVEszs = num2str(SVEsz);
        ResName = ['Res', SVEszs, 'x', SVEszs];
        fileNameWExt = [ResName, '_', BC, '_KMuENuLambda_stat.txt'];
        B = readmatrix(fileNameWExt);
        [m, n] = size(B);
        matVals{SVEszi} = B(:, 2:n);
    end

    for fi = 1:nFields
        fields = fieldNames{fi};
        fileName = ['C_', fields, '_BC_', BC, '.txt'];
        fid = fopen(fileName, 'w');
    
        fprintf(fid, 'SVESizeInv');
        for ni = 1:nnames
            fprintf(fid, '\t%s', names{ni});
        end
        for SVEszi = 1:nSVEsz
            vals = matVals{SVEszi}(:, fi);
            fprintf(fid, '\n%d', SVESizeInvs(SVEszi));
            for ni = 1:nnames
                fprintf(fid, '\t%g', vals(ni));
            end
        end
        fclose(fid);
    end
end


fieldNames = {'ARO', 'ANZ'};
nFields = length(fieldNames);
%%%% Anisotropy indices
matVals = cell(0);

for BCi = 1:3
    BC = BCs{BCi};
    for SVEszi = 1:nSVEsz
        SVEsz = SVESizeInvs(SVEszi);
        SVEszs = num2str(SVEsz);
        ResName = ['Res', SVEszs, 'x', SVEszs];
        fileNameWExt = [ResName, '_', BC, '_As_stat.txt'];
        B = readmatrix(fileNameWExt);
        [m, n] = size(B);
        matVals{SVEszi} = B(:, 2:n);
    end

    for fi = 1:nFields
        fields = fieldNames{fi};
        fileName = ['As_', fields, '_BC_', BC, '.txt'];
        fid = fopen(fileName, 'w');
    
        fprintf(fid, 'SVESizeInv');
        for ni = 1:nnames
            fprintf(fid, '\t%s', names{ni});
        end
        for SVEszi = 1:nSVEsz
            vals = matVals{SVEszi}(:, fi);
            fprintf(fid, '\n%d', SVESizeInvs(SVEszi));
            for ni = 1:nnames
                fprintf(fid, '\t%g', vals(ni));
            end
        end
        fclose(fid);
    end
end


fieldNames = {'vf'};
nFields = length(fieldNames);
%%%% Anisotropy indices
matVals = cell(0);

for SVEszi = 1:nSVEsz
    SVEsz = SVESizeInvs(SVEszi);
    SVEszs = num2str(SVEsz);
    ResName = ['Res', SVEszs, 'x', SVEszs];
    fileNameWExt = [ResName, '_vf_stat.txt'];
    B = readmatrix(fileNameWExt);
    [m, n] = size(B);
    matVals{SVEszi} = B(:, 2:n);
end

for fi = 1:nFields
    fields = fieldNames{fi};
    fileName = ['scalars_', fields, '.txt'];
    fid = fopen(fileName, 'w');

    fprintf(fid, 'SVESizeInv');
    for ni = 1:nnames
        fprintf(fid, '\t%s', names{ni});
    end
    for SVEszi = 1:nSVEsz
        vals = matVals{SVEszi}(:, fi);
        fprintf(fid, '\n%d', SVESizeInvs(SVEszi));
        for ni = 1:nnames
            fprintf(fid, '\t%g', vals(ni));
        end
    end
    fclose(fid);
end
