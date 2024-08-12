function [KmuEnulamba, As, C, CIso] = Compute_k_A_from_C(CIn)
printResults = 0;
if (nargin < 1)
    CIn = [4, 1, 0; 1, 5, 0; 0, 0, 2];
    CIn = [72.4482	14.439	0.346849;	14.439	67.6221	0.555398;	0.346854	0.555396	26.485];
    printResults = 1;
end
[m, n] = size(CIn);
if (n == 9)
    C = zeros(3, 3);
    cntr = 0;
    for i = 1:3
        for j = 1:3
            cntr = cntr + 1;
            C(i, j) = CIn(cntr);
        end
    end
else
    C = CIn;
end
CIso = zeros(3, 3);
C1122 = C(1, 1) + C(2, 2);
C12 = C(1, 2);
C33 = C(3, 3);
CIso(1, 1) = 0.375 * C1122 + 0.25 * C12 + 0.5 * C33;
CIso(2, 2) = CIso(1, 1);
CIso(1, 2) = 0.125 * C1122 + 0.75 * C12 - 0.5 * C33;
CIso(2, 1) = CIso(1, 2);
CIso(3, 3) = 0.5 * (CIso(1, 1) - CIso(1, 2));

D = inv(C);
D1122 = D(1, 1) + D(2, 2);
D12 = D(1, 2);
D33 = D(3, 3);
DIso(1, 1) = 0.375 * D1122 + 0.25 * D12 + 0.125 * D33;
DIso(2, 2) = DIso(1, 1);
DIso(1, 2) = 0.125 * D1122 + 0.75 * D12 - 0.125 * D33;
DIso(2, 1) = DIso(1, 2);
DIso(3, 3) = 2.0 * (DIso(1, 1) - DIso(1, 2));

KV2D = 0.5 * (CIso(1, 1) + CIso(1, 2));
GV2D = CIso(3, 3);

KR2D = 0.5 / (DIso(1, 1) + DIso(1, 2));
GR2D = 1.0 / DIso(3, 3);

ARO = KV2D / KR2D + 2.0 * GV2D / GR2D - 3.0;


delC11 = C(1, 1) - CIso(1, 1);
delC22 = C(2, 2) - CIso(1, 1);
delC12 = C12 - CIso(1, 2);
delC13 = C(1, 3);
delC23 = C(2, 3);
delC33 = C33 - CIso(3, 3);
normC = sqrt(C(1, 1) * C(1, 1) + C(2, 2) * C(2, 2) + 2.0 * C12 * C12 + 4.0 * ...
    (C(1, 3) * C(1, 3) + C(2, 3) * C(2, 3) + C(3, 3) * C(3, 3)));
normdelC = sqrt(delC11 * delC11 + delC22 * delC22 + 2.0 * delC12 * delC12 + 4.0 * (delC13 * delC13 + delC23 * delC23 + delC33 * delC33));
ANZ = normdelC / normC;
mode = 21;
[E, nu, lambda] = Kmu2_Enulambda(KR2D, GR2D, mode);
KmuEnulamba = [KR2D, GR2D, E, nu, lambda];
As = [ARO, ANZ];

if (printResults)
    KV = 0.25 * (C1122 + 2.0 * C12);
    muV = 0.125 * (C1122 - 2.0 * C12 + 4.0 * C33);
    KR = 1.0 / (D1122 + 2.0 * D12);
    muR = 2.0 / (D1122 - 2.0 * D12 + D33);
    ARO_easyCalc = KV / KR + 2.0 * muV / muR - 3;
    fprintf(1, 'ARO calclated from (3) = %g\n', ARO_easyCalc);
    fprintf(1, 'ARO calclated from longer path (isotropic C first calculated) = %g', ARO);
    fprintf(1, '\nANZ = %g\n\n', ANZ);
end