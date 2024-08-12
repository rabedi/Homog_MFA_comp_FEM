function [res, nms, resVoigt, resResuss, resHSp, resHSm] = SiCB4(cI)
if nargin < 1
    cI = 0.23905;
end
imf = IsoMeanField;
parasM = [45.50, 0.18]; % SiC
parasI = [455.0, 0.18]; % SiC
mode = 21;
input_type = 1;
[K_mu_E_nu_lambdas_AllMethods, objout] = imf.Compute_PropLastcI(cI, parasM, parasI, mode, input_type);
resVoigt = K_mu_E_nu_lambdas_AllMethods{1};
resResuss = K_mu_E_nu_lambdas_AllMethods{2};
resHSp = K_mu_E_nu_lambdas_AllMethods{3};
resHSm = K_mu_E_nu_lambdas_AllMethods{4};


res{1} = K_mu_E_nu_lambdas_AllMethods{7};
res{2} = K_mu_E_nu_lambdas_AllMethods{8};
res{3} = K_mu_E_nu_lambdas_AllMethods{9};
nms{1} = imf.names{7};
nms{2} = imf.names{8};
nms{3} = imf.names{9};
