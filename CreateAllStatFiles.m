function CreateAllStatFiles(SVESizeInvs)

if (nargin < 1)
    SVESizeInvs = [16, 8, 4, 2, 1];
end

BCs = {'disp', 'mixed', 'trac', 'MT', 'SC', 'Diff'};
nBCs = length(BCs);

nSVEsz = length(SVESizeInvs);
for SVEszi = 1:nSVEsz
    SVEsz = SVESizeInvs(SVEszi);
    SVEszs = num2str(SVEsz);
    ResName = ['Res', SVEszs, 'x', SVEszs];
    for BCi = 1:nBCs
        BC = BCs{BCi};
        fileNameWOExt = [ResName, '_', BC, '_KMuENuLambda'];
        CreateStatOneFile(fileNameWOExt);
        if (BCi < 4)
            fileNameWOExt = [ResName, '_', BC, '_As'];
            CreateStatOneFile(fileNameWOExt);
        end
    end
    fileNameWOExt = [ResName, '_vf'];
    CreateStatOneFile(fileNameWOExt);
end
    


