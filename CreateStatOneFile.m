function CreateStatOneFile(fileNameInputWOExt)
if (nargin < 1)
    fileNameInputWOExt = 'Res1x1_Diff_KMuENuLambda';
end

fileNameIn = [fileNameInputWOExt, '.txt'];
A = readmatrix(fileNameIn);
[m, n] = size(A);
for ci = 1:n
    vec = A(:, ci);
    meanV(ci) = mean(vec);
    minV(ci) = min(vec);
    maxV(ci) = max(vec);
    stdV(ci) = std(vec);
    coefvV(ci) = stdV / meanV;
end
fileNameOut = [fileNameInputWOExt, '_stat.txt'];
fid = fopen(fileNameOut, 'w');
fprintf(fid, 'n');
for ci = 1:n
    fprintf(fid, '\t%d', m);
end

fprintf(fid, '\nmean');
for ci = 1:n
    fprintf(fid, '\t%d', meanV(ci));
end
fprintf(fid, '\nstd');
for ci = 1:n
    fprintf(fid, '\t%d', stdV(ci));
end
fprintf(fid, '\ncoefVar');
for ci = 1:n
    fprintf(fid, '\t%d', coefvV(ci));
end
fprintf(fid, '\nmin');
for ci = 1:n
    fprintf(fid, '\t%d', minV(ci));
end
fprintf(fid, '\nmax');
for ci = 1:n
    fprintf(fid, '\t%d', maxV(ci));
end
fprintf(fid, '\nmean-std');
for ci = 1:n
    fprintf(fid, '\t%d', meanV(ci) - stdV(ci));
end
fprintf(fid, '\nmean+std');
for ci = 1:n
    fprintf(fid, '\t%d', meanV(ci) + stdV(ci));
end
fclose(fid);
