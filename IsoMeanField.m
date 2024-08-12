classdef IsoMeanField
    properties
        verbose = 1; % outputs on screen are one
        cIVerbose = 0.02; % inclusion volume fraction for which we want the verbose output
        % mode 3: 3D
        % mode 21: plane strain
        % mode 22: plane stress
        mode = 21;
        %   1: E and nu provided
        %   2: K and mu provided
        input_type = 1;
        % number of inclusions (number of materials = number of inclusions
        % + 1)
        num_inc = 1;
        % relative volume fracture of inclusion phases among themselves:
        % for 1 inclusion type, this is 1
        % for two inclusions with 50% each inclusion type, this is [0.5,
        % 0.5], if the first one is 60% ,this is [0.6, 0.4]
        xis = 1;

        matrix_pars; % [1x2] = stores (K, mu) OR (E, nu) of matrix based on input_type
        inclusion_paras; % [num_inc x 2] stores (K, mu) OR (E, nu) of each inclusion

        %%% computed
        % this is a (num_mat x 5) matrix, first row is the matrix material, next rows are inclusions
        % based on input_type, output of 5 columns only 2 are read. Other 3 are computed using functions
        % Kmu2_Enulambda o Enu2_Kmulambda
        K_mu_E_nu_lambdas;
        s12_M; % this are Eshelby 1 and 2 for material 1 (matrix), precomputed


        %%% indices and names for output of different methods
        i_Voigt = 1;
        i_Reuss = 2;
        i_HS_plus = 3;
        i_HS_minus = 4;
        i_DD_epsilon = 5;
        i_DD_sigma = 6;
        i_MT_epsilon = 7;
        i_SC = 8;
        i_Diff_epsilon = 9;
        i_Diff_sigma = 10;
        i_MT_sigma = 11;
        names = {'Voigt', 'Reuss', 'HS+', 'HS-', 'DD-espilon', 'DD-sigma', 'MT', 'SC', 'Diff', 'Diff-sigma', 'MT-sigma'};

        % number of C types computed
        i_max = 11;

        % for SC the initial implementation did not have the correction
        % term. The flag below allows it to be with the correction (i.e.
        % like HS). If 0, it does not have <A>^-1 term (correction)
        % UPDATE: 03/14/2024 -> Do not change this to 0. It seems the
        % update make the convergence smoother and even some cases prevent
        % nonconvergence
        use_corrections4SC = 1;
    end
    methods
        function objout = Read_IsoMeanField(obj, fid)
            buf = fscanf(fid, '%s', 1);
            obj.mode = fscanf(fid, '%d', 1);

            buf = fscanf(fid, '%s', 1);
            obj.input_type = fscanf(fid, '%d', 1);

            buf = fscanf(fid, '%s', 1);
            obj.matrix_pars(1) = fscanf(fid, '%f', 1);
            obj.matrix_pars(2) = fscanf(fid, '%f', 1);

            buf = fscanf(fid, '%s', 1);
            obj.num_inc = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 3);
            obj.xis = zeros(obj.num_inc, 1);
            obj.inclusion_paras = zeros(obj.num_inc, 2);
            for i = 1:obj.num_inc
                obj.xis(i) = fscanf(fid, '%f', 1);
                obj.inclusion_paras(i, 1) = fscanf(fid, '%f', 1);
                obj.inclusion_paras(i, 2) = fscanf(fid, '%f', 1);
            end
            objout = Initialize(obj);
        end
        function objout = Initialize(obj)
            [m, n] = size(obj.inclusion_paras);
            obj.num_inc = m;
            if (size(obj.xis) ~= obj.num_inc)
                xis = obj.xis
                num_inc = obj.num_inc
                fprintf(1, 'error in size of xis (should be equal to num_inc)\n');
                pause
            end
            sum_xis = sum(obj.xis);
            if (sum_xis < 1e-4)
                obj.xis(1) = 1.0;
                sum_xis = sum(obj.xis);
            end
            obj.xis = 1.0 / sum_xis  * obj.xis;

            num_mat = obj.num_inc + 1;
            vals = zeros(num_mat, 2);
            vals(1, 1) = obj.matrix_pars(1);
            vals(1, 2) = obj.matrix_pars(2);

            for i = 1:obj.num_inc
                ip1 = i + 1;
                vals(ip1, 1) = obj.inclusion_paras(i, 1);
                vals(ip1, 2) = obj.inclusion_paras(i, 2);
            end

            %   1: E and nu provided
            if (obj.input_type == 1)
                for i = 1:num_mat
                    E = vals(i, 1);
                    nu = vals(i, 2);
                    [K, mu, lambda] = Enu2_Kmulambda(E, nu, obj.mode);
                    obj.K_mu_E_nu_lambdas(i, 1) = K;
                    obj.K_mu_E_nu_lambdas(i, 2) = mu;
                    obj.K_mu_E_nu_lambdas(i, 3) = E;
                    obj.K_mu_E_nu_lambdas(i, 4) = nu;
                    obj.K_mu_E_nu_lambdas(i, 5) = lambda;
                end
                %   2: K and mu provided
            elseif (obj.input_type == 2)
                for i = 1:num_mat
                    K = vals(i, 1);
                    mu = vals(i, 2);
                    [E, nu, lambda] = Kmu2_Enulambda(K, mu, obj.mode);
                    obj.K_mu_E_nu_lambdas(i, 1) = K;
                    obj.K_mu_E_nu_lambdas(i, 2) = mu;
                    obj.K_mu_E_nu_lambdas(i, 3) = E;
                    obj.K_mu_E_nu_lambdas(i, 4) = nu;
                    obj.K_mu_E_nu_lambdas(i, 5) = lambda;
                end
            end
            K = obj.K_mu_E_nu_lambdas(1, 1);
            mu = obj.K_mu_E_nu_lambdas(1, 2);
            [obj.s12_M(1), obj.s12_M(2)] = Kmu2_s1s2(K, mu, obj.mode);

            objout = obj;
        end

        % cI is the absolute volume fraction, matrix will be cM = 1 - cI,
        % inclusion i absolute volume fraction will be xis(i) * cI,
        % absolute volume fractions will be listed as [cM, cInclusion1,
        % cInclusion2, ...] = [1 - cI, cI * xis(1), ...]
        function fs = GetAbsoluteVolumeFractions(obj, cI)
            fs = zeros(obj.num_inc + 1, 1);
            fs(1) = 1.0 - cI;
            for i = 1:obj.num_inc
                fs(i + 1) = cI * obj.xis(i);
            end
        end
        % this is the same as above for all methods (DD, MT, ...) except
        % Differential wherein this value is for inclusions only and takes
        % the value xi(i)/(1 - cI)
        function effective_fs = GetEffectiveVolumeFractions(obj, cI, isDifferential)
            if (isDifferential == 0)
                effective_fs = obj.GetAbsoluteVolumeFractions(cI);
                return;
            end
            effective_fs = zeros(obj.num_inc + 1, 1);
            effective_fs(1) = 1.0 - cI; % will not be used anyway
            cImax = min(cI, 1.0 - 1e-6);
            denomInv = 1.0 / (1.0 - cImax);
            for i = 1:obj.num_inc
                effective_fs(i + 1) = denomInv * obj.xis(i);
            end
        end
        % As and Bs are [num_mat x 2] matrices, row i -> material i
        % (numbering starting from matrix), column 1 -> bulk, column 2 ->
        % shear

        % Kmu0 (2 x 1) are bulk and shear properties of reference material.
        % 0 is:
        %:  DD, MT  : matrix
        %   HS      : max and min K/mu for HS +/- (or in general any 0)
        %   SC      : an ongoing homogenized K, mu that will be updated
        %   Diff    : current matrix (C*(cI))

        % mat0 is matrix says whether 0 is matrix, which obviously only
        % true for DD and MT

        % computeB: computes B too
        function [As, Bs] = ComputeConcentrations_AB(obj, Kmu0, mat0_isMatrix, computeB)
            nmat = obj.num_inc + 1;
            As = ones(nmat, 2);
            if (computeB)
                Bs = As;
            else
                Bs = [];
            end
            % st: starting point for computing A's (if 0 is m first A and B
            % are 1)
            st = 2;
            s12 = obj.s12_M;
            if (~mat0_isMatrix)
                [s1, s2] = Kmu2_s1s2(Kmu0(1), Kmu0(2), obj.mode);
                s12 = [s1, s2];
                st = 1;
            end
            for i = st:nmat
                for beta = 1:2 % mode: beta = 1 (bulk), 2 (shear)
                    % ratio of material i stiffness (K or mu) to
                    % corresponding value of reference material
                    r_i0_beta = obj.K_mu_E_nu_lambdas(i, beta) / Kmu0(beta);
                    Ai0 = 1.0 / (1.0 + s12(beta) * (r_i0_beta - 1.0));
                    As(i, beta) = Ai0;
                    if (computeB)
                        Bs(i, beta) = r_i0_beta * Ai0;
                    end
                end
            end
        end

        % cI: volume fraction of all inclusions combined

        % computeD: compute compliance (sigma based approaches, only make
        % sense for DD and MT)

        % compute_corrections: computes corrections (besides base terms
        % below) -> corrections are computed for MT and HS

        % Kmu0: 0 material (for Ai,0 and Bi,0 values)
        % mat0_isMatrix: says if 0 is matrix

        % KmuB: is KmuB, they are KmuMatrix unless dealing with
        % Differential method where in these values need to be provided
        % isDifferential: say is the method is differential: In addition to above change, it also does not add CB and DB terms below

        % C_term = CB + (sum_i=2 fi(Ci - CB)Ai,0) correctionC
        % correctionC =  (sum_i=1 fi Ai,0)^-1
        %       C_term_wo_correction without correction (DD, SC, Diff)
        %       C_term_w_correction with correction (MT, HS)
        %   CB term is not added for differential

        % D_term = DB + (sum_i=2 fi(Di - DB)Bi,0) correctionB
        % correctionD =  (sum_i=1 fi Bi,0)^-1
        %       D_term_wo_correction without correction (DD, SC, Diff)
        %       D_term_w_correction with correction (MT, HS)
        %   DB term is not added for differential

        function [C_term_wo_correction, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI, Kmu0, KmuB, mat0_isMatrix, isDifferential, computeD, compute_corrections)
            C_term_wo_correction = zeros(1, 2);
            C_term_w_correction = zeros(1, 2);
            D_term_wo_correction = zeros(1, 2);
            D_term_w_correction = zeros(1, 2);
            if (~isDifferential)
                KmuB(1) = obj.K_mu_E_nu_lambdas(1, 1);
                KmuB(2) = obj.K_mu_E_nu_lambdas(1, 2);
            end
            [As, Bs] = ComputeConcentrations_AB(obj, Kmu0, mat0_isMatrix, computeD);
            % compute volume fractions that are used in formulas below
            effective_fs = GetEffectiveVolumeFractions(obj, cI, isDifferential);
            nmat = obj.num_inc + 1;
            for beta = 1:2 % mode
                CB = KmuB(beta);
                for i = 2:nmat
                    f = effective_fs(i);
                    Ai = As(i, beta);
                    Ci_CB = obj.K_mu_E_nu_lambdas(i, beta) - CB;
                    C_term_wo_correction(beta) = C_term_wo_correction(beta) + f * Ci_CB * Ai;
                end
            end
            if (compute_corrections)
                for beta = 1:2
                    correction = 0;
                    for i = 1:nmat
                        f = effective_fs(i);
                        Ai = As(i, beta);
                        correction = correction + f * Ai;
                    end
                    correction = 1.0 / correction;
                    C_term_w_correction(beta) = C_term_wo_correction(beta) * correction;
                    if (~isDifferential)
                        C_term_w_correction(beta) = C_term_w_correction(beta) + KmuB(beta);
                    end
                end
            end
            if (~isDifferential)
                for beta = 1:2
                    C_term_wo_correction(beta) = C_term_wo_correction(beta) + KmuB(beta);
                end
            end

            if (computeD == 0)
                return;
            end

            KmuBInv = 1.0 ./ KmuB;
            for beta = 1:2 % mode
                DB = KmuBInv(beta);
                for i = 2:nmat
                    f = effective_fs(i);
                    Bi = Bs(i, beta);
                    Di_DB = 1.0 / obj.K_mu_E_nu_lambdas(i, beta) - DB;
                    D_term_wo_correction(beta) = D_term_wo_correction(beta) + f * Di_DB * Bi;
                end
            end
            if (compute_corrections)
                for beta = 1:2
                    correction = 0;
                    for i = 1:nmat
                        f = effective_fs(i);
                        Bi = Bs(i, beta);
                        correction = correction + f * Bi;
                    end
                    correction = 1.0 / correction;
                    D_term_w_correction(beta) = D_term_wo_correction(beta) * correction;
                    if (~isDifferential)
                        D_term_w_correction(beta) = D_term_w_correction(beta) + KmuBInv(beta);
                    end
                end
            end
            if (~isDifferential)
                for beta = 1:2
                    D_term_wo_correction(beta) = D_term_wo_correction(beta) + KmuBInv(beta);
                end
            end
        end

        %%%%% calcuating different type of stiffness
        % Voigt and Reuss limits
        function [C_Reuss, C_Voigt] = Compute_Reuss_Voigt(obj, cI)
            C_Reuss = zeros(1, 2);
            C_Voigt = zeros(1, 2);
            % compute volume fractions that are used in formulas below
            effective_fs = GetEffectiveVolumeFractions(obj, cI, 0);
            nmat = obj.num_inc + 1;
            for beta = 1:2
                for i = 1:nmat
                    f = effective_fs(i);
                    C = obj.K_mu_E_nu_lambdas(i, beta);
                    C_Reuss(beta) = C_Reuss(beta) + f / C;
                    C_Voigt(beta) = C_Voigt(beta) + f * C;
                end
                C_Reuss(beta) = 1.0 / C_Reuss(beta);
            end
        end

        % Hashin-Shtrikman limits
        function [C_HS_minus, C_HS_plus] = Compute_HS_limits(obj, cI)
            Kmu0s = cell(2, 1);
            for mp = 1:2
                Kmu0s{mp} = zeros(1, 2);
            end
            KmuB = zeros(1, 2);

            % finding min and max of K, mu of materials
            for beta = 1:2
                Cs = obj.K_mu_E_nu_lambdas(:, beta);
                Csm = min(Cs);
                CsM = max(Cs);
                Kmu0s{1}(beta) = Csm;
                Kmu0s{2}(beta) = CsM;
                KmuB(beta) = obj.K_mu_E_nu_lambdas(1, beta);
            end
            if ((obj.verbose) && (abs(cI - obj.cIVerbose) < 1e-3))
                fprintf(1, 'cI = %g\n', cI);
                fprintf(1, 'Matrix     : KM = %f, muM = %f, EM = %f, nuM = %f\n', obj.K_mu_E_nu_lambdas(1, 1), obj.K_mu_E_nu_lambdas(1, 2), obj.K_mu_E_nu_lambdas(1, 3), obj.K_mu_E_nu_lambdas(1, 4));
                fprintf(1, 'Inclusion 1: KI = %f, muI = %f, EI = %f, muI = %f\n', obj.K_mu_E_nu_lambdas(2, 1), obj.K_mu_E_nu_lambdas(2, 2), obj.K_mu_E_nu_lambdas(2, 3), obj.K_mu_E_nu_lambdas(2, 4));
                fprintf(1, 'Ratios:\n');
                kM2KI = obj.K_mu_E_nu_lambdas(1, 1) / obj.K_mu_E_nu_lambdas(2, 1);
                muM2muI = obj.K_mu_E_nu_lambdas(1, 2) / obj.K_mu_E_nu_lambdas(2, 2);
                fprintf(1, 'KM/KI = %f, muM/muI = %f\n', kM2KI, muM2muI);
                fprintf(1, 'HS-(min) material properties:  : K0-HS- = %f, mu0-HS- = %f\n', Kmu0s{1}(1), Kmu0s{1}(2));
                fprintf(1, 'HS+(max) material properties:  : K0-HS+ = %f, mu0-HS+ = %f\n', Kmu0s{2}(1), Kmu0s{2}(2));
            end

            isDifferential = 0;
            computeD = 0;
            compute_corrections = 1;
            mat0_isMatrix = 0;
            for mp = 1:2
                Kmu0 = Kmu0s{mp};
                [C_term_wo_correction, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI, Kmu0, KmuB, mat0_isMatrix, isDifferential, computeD, compute_corrections);
                if (mp == 1)
                    C_HS_minus = C_term_w_correction;
                else
                    C_HS_plus = C_term_w_correction;
                end
            end
        end

        % DD and Mori-Tanaka, epsilon based, sigma based
        function [C_DD_epsilon, C_DD_sigma, C_MT_epsilon, C_MT_sigma] = Compute_DD_MT(obj, cI)
            KmuB = zeros(1, 2);
            for beta = 1:2
                KmuB(beta) = obj.K_mu_E_nu_lambdas(1, beta);
            end
            Kmu0 = KmuB;
            mat0_isMatrix = 1;
            isDifferential = 0;
            computeD = 1;
            compute_corrections = 1;
            [C_term_wo_correction, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI, Kmu0, KmuB, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            C_DD_epsilon = C_term_wo_correction;
            C_DD_sigma = 1.0 ./ D_term_wo_correction;
            C_MT_epsilon = C_term_w_correction;
            C_MT_sigma = 1.0 ./ D_term_w_correction;
        end

        % Self-Consistent (SC) with provided initial guess
        function [C_SC, converged] = Compute_CS(obj, cI, maxNumIteration, relTol, C_SC_ini)
            if (nargin < 3)
                maxNumIteration = 200;
            end
            if (nargin < 4)
                relTol = 1e-3;
            end
            if (nargin < 5)
                C_SC_ini = [];
            end
            if (isempty(C_SC_ini))
                [C_HS_minus, C_HS_plus] = Compute_HS_limits(obj, cI);
                C_SC_ini = 0.5 * (C_HS_minus + C_HS_plus);
            end
            KmuB = zeros(1, 2);
            for beta = 1:2
                KmuB(beta) = obj.K_mu_E_nu_lambdas(1, beta);
            end
            mat0_isMatrix = 0;
            isDifferential = 0;
            computeD = 0;

            compute_corrections = obj.use_corrections4SC;

            converged = 0;
            iter = 0;
            Kmu0 = C_SC_ini;
            while ((converged == 0) && (iter < maxNumIteration))
                [Kmu0_new, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI, Kmu0, KmuB, mat0_isMatrix, isDifferential, computeD, compute_corrections);
                if (obj.use_corrections4SC)
                    Kmu0_new = C_term_w_correction;
                end
                converged = 1;
                for beta = 1:2
                    valNew = Kmu0_new(beta);
                    valOld = Kmu0(beta);
                    relError = abs(valNew - valOld) / max(valNew, valOld);
                    if (relError > relTol)
                        converged = 0;
                        break;
                    end
                end
                Kmu0 = Kmu0_new;
                iter = iter + 1;
            end
            C_SC = Kmu0_new;
            %% This part is not needed and is just included to demonstrate that for SC average of Ai and Bi is 1.
            % Comment out to show this ...
            % compute sum of fi Ai,0_dil ...
            if (0)
                computeB = 1;
                [As, Bs] = ComputeConcentrations_AB(obj, Kmu0_new, mat0_isMatrix, computeB);
                effective_fs = GetEffectiveVolumeFractions(obj, cI, isDifferential);
                nmat = obj.num_inc + 1;
                As = ones(nmat, 2);

                sum_fA = effective_fs(1) * As(1,:);
                sum_fB = effective_fs(1) * Bs(1,:);
                for i = 2:nmat
                    sum_fA = sum_fA + effective_fs(i) * As(i,:);
                    sum_fB = sum_fB + effective_fs(i) * Bs(i,:);
                end
                fprintf(1, '%g\t%g\t%g\t%g\t%g\n', cI, sum_fA(1), sum_fA(2), sum_fB(1), sum_fB(2));
            end
        end

        % Differential: Computes differential C at cI+delcI from it's
        % previous value at cI, using RK4 method
        function C_Diff_cI_plus_delcI = Compute_Differential_C_based_cI_plus_del_cI_Aux(obj, cI, delcI, C_Diff_cI)
            % dC/dc = f(c, C), C = stiffness, c = concentration of inclusion phase (cI) ->
            % C(c+delC) = C(c) + delc/6 * (k1 + 2k2 + 2k3 + k4)
            % k1 = f(c, C)
            % k2 = f(c + delc / 2, C + delc k1 / 2)
            % k3 = f(c + delc / 2, C + delc k2 / 2)
            % k4 = f(c + delc    , C + delc k3    )
            half_delcI = 0.5 * delcI;
            compute_corrections = 0;
            computeD = 0;
            isDifferential = 1;
            mat0_isMatrix = 0;

            % calculating k1
            cI2Use = cI;
            Kmu0B = C_Diff_cI;
            [k1, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k2
            cI2Use = cI + half_delcI;
            Kmu0B = C_Diff_cI + half_delcI * k1;
            [k2, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k3
            % cI2Use = cI + half_delcI; -> not changed
            Kmu0B = C_Diff_cI + half_delcI * k2;
            [k3, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k4
            cI2Use = cI + delcI;
            if (1.0 - cI2Use < 1e-4) % too close to 1.0 division by zero problem in 1/(1 - cI) term
                k4 = k3;
            else
                Kmu0B = C_Diff_cI + delcI * k3;
                [k4, C_term_w_correction, D_term_wo_correction, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            end
            C_Diff_cI_plus_delcI = C_Diff_cI + delcI / 6.0 * (k1 + 2.0 * (k2 + k3) + k4);
        end

        function D_Diff_cI_plus_delcI = Compute_Differential_D_based_cI_plus_del_cI_Aux(obj, cI, delcI, D_Diff_cI)
            % dD/dc = f(c, D), D = compliance, c = concentration of inclusion phase (cI) ->
            % D(c+delC) = D(c) + delc/6 * (k1 + 2k2 + 2k3 + k4)
            % k1 = f(c, D)
            % k2 = f(c + delc / 2, D + delc k1 / 2)
            % k3 = f(c + delc / 2, D + delc k2 / 2)
            % k4 = f(c + delc    , D + delc k3    )
            half_delcI = 0.5 * delcI;
            compute_corrections = 0;
            computeD = 1;
            isDifferential = 1;
            mat0_isMatrix = 0;

            % calculating k1
            cI2Use = cI;
            Kmu0B = 1.0 ./ D_Diff_cI;
            [C_term_wo_correction, C_term_w_correction, k1, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k2
            cI2Use = cI + half_delcI;
            Kmu0B = 1.0 ./ (D_Diff_cI + half_delcI * k1);
            [C_term_wo_correction, C_term_w_correction, k2, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k3
            % cI2Use = cI + half_delcI; -> not changed
            Kmu0B = 1.0 ./ (D_Diff_cI + half_delcI * k2);
            [C_term_wo_correction, C_term_w_correction, k3, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            % calculating k4
            cI2Use = cI + delcI;
            if (1.0 - cI2Use < 1e-4) % too close to 1.0 division by zero problem in 1/(1 - cI) term
                k4 = k3;
            else
                Kmu0B = 1.0 ./ (D_Diff_cI + delcI * k3);
                [C_term_wo_correction, C_term_w_correction, k4, D_term_w_correction] = Compute_C_D_terms(obj, cI2Use, Kmu0B, Kmu0B, mat0_isMatrix, isDifferential, computeD, compute_corrections);
            end
            D_Diff_cI_plus_delcI = D_Diff_cI + delcI / 6.0 * (k1 + 2.0 * (k2 + k3) + k4);
        end

        function [cIs, K_mu_E_nu_lambdas_DiffEpsilon, K_mu_E_nu_lambdas_DiffSigma] = Compute_Differential(obj, delcI, cIMax, compute_sigma_versionToo)
            if (nargin < 2)
                delcI = 0.005;
            end
            if (nargin < 3)
                cIMax = 1.0;
            end
            if (nargin < 4)
                compute_sigma_versionToo = 0;
            end
            K_mu_E_nu_lambdas_DiffSigma = [];
            sz = ceil(cIMax / delcI); % number of segments
            delcI = cIMax / sz;
            cIs = 0:delcI:cIMax;
            sz = sz + 1; % number of points
            K_mu_E_nu_lambdas_DiffEpsilon = zeros(sz, 5);
            K_mu_E_nu_lambdas_DiffEpsilon(1,:) = obj.K_mu_E_nu_lambdas(1, :);
            for i = 2:sz
                cI = cIs(i - 1);
                C_Diff_cI = K_mu_E_nu_lambdas_DiffEpsilon(i - 1,1:2);
                C_Diff_cI_plus_delcI = Compute_Differential_C_based_cI_plus_del_cI_Aux(obj, cI, delcI, C_Diff_cI);
                K = C_Diff_cI_plus_delcI(1);
                mu = C_Diff_cI_plus_delcI(2);
                [E, nu, lambda] = Kmu2_Enulambda(K, mu, obj.mode);
                K_mu_E_nu_lambdas_DiffEpsilon(i,:) = [K, mu, E, nu, lambda];
            end
            if (compute_sigma_versionToo == 0)
                return;
            end

            K_mu_E_nu_lambdas_DiffSigma = zeros(sz, 5);
            K_mu_E_nu_lambdas_DiffSigma(1,:) = obj.K_mu_E_nu_lambdas(1, :);
            for i = 2:sz
                cI = cIs(i - 1);
                C_Diff_cI = K_mu_E_nu_lambdas_DiffSigma(i - 1,1:2);
                D_Diff_cI = 1.0 ./ C_Diff_cI;
                D_Diff_cI_plus_delcI = Compute_Differential_D_based_cI_plus_del_cI_Aux(obj, cI, delcI, D_Diff_cI);
                K = 1.0 / D_Diff_cI_plus_delcI(1);
                mu = 1.0 / D_Diff_cI_plus_delcI(2);
                [E, nu, lambda] = Kmu2_Enulambda(K, mu, obj.mode);
                K_mu_E_nu_lambdas_DiffSigma(i,:) = [K, mu, E, nu, lambda];
            end
        end

        % compute all but Differential
        function K_mus = Compute_All_but_Differential(obj, cI, maxNumIteration, relTol, C_SC_ini, compute_sigma_versionToo)
            if (nargin < 3)
                maxNumIteration = 100;
            end
            if (nargin < 4)
                C_SC_ini = [];
            end
            if (nargin < 5)
                compute_sigma_versionToo = 0;
            end
            %            K_mu_E_nu_lambds_All_but_Diff = cell(obj.i_max, 1);
            K_mus = cell(obj.i_max, 1);

            [K_mus{obj.i_Reuss}, K_mus{obj.i_Voigt}] = Compute_Reuss_Voigt(obj, cI);
            [K_mus{obj.i_HS_minus}, K_mus{obj.i_HS_plus}] = Compute_HS_limits(obj, cI);
            [K_mus{obj.i_DD_epsilon}, K_mus{obj.i_DD_sigma}, K_mus{obj.i_MT_epsilon}, K_mus{obj.i_MT_sigma}] = ...
                Compute_DD_MT(obj, cI);
            if (isempty(C_SC_ini))
                C_SC_ini = K_mus{obj.i_MT_epsilon};
            end
            [K_mus{obj.i_SC}, converged] = Compute_CS(obj, cI, maxNumIteration, relTol, C_SC_ini);
            if (converged == 0)
                fprintf(1, 'For cI = %f convergence not achieved for CS\n', cI);
                %                pause;
            end

            for cm = 1:obj.i_SC
                K = K_mus{cm}(1);
                mu = K_mus{cm}(2);
                [E, nu, lambda] = Kmu2_Enulambda(K, mu, obj.mode);
                K_mus{cm}(3) = E;
                K_mus{cm}(4) = nu;
                K_mus{cm}(5) = lambda;
            end
        end

        function [cIs, K_mu_E_nu_lambdas_AllMethods, namesOut] = Compute_All_Methods_All_cIs(obj, delcI, cIMax, maxNumIteration, relTol, compute_sigma_versionToo)
            use_prev_CS_C_4_CS_ini_guess = 1;
            namesOut = obj.names;
            K_mu_E_nu_lambdas_AllMethods = cell(obj.i_max, 1);
            [cIs, K_mu_E_nu_lambdas_AllMethods{obj.i_Diff_epsilon}, K_mu_E_nu_lambdas_AllMethods{obj.i_Diff_sigma}] = Compute_Differential(obj, delcI, cIMax, compute_sigma_versionToo);
            num_cIs = length(cIs);
            for cm = 1:obj.i_SC
                K_mu_E_nu_lambdas_AllMethods{cm} = zeros(num_cIs, 5);
            end
            C_SC_ini = [];
            for ci = 1:num_cIs
                cI = cIs(ci);
                if (use_prev_CS_C_4_CS_ini_guess && (ci > 1))
                    C_SC_ini(1) = K_mu_E_nu_lambdas_AllMethods{obj.i_SC}(ci - 1, 1);
                    C_SC_ini(2) = K_mu_E_nu_lambdas_AllMethods{obj.i_SC}(ci - 1, 2);
                end
                K_mu_E_nu_lambds_All_but_Diff = Compute_All_but_Differential(obj, cI, maxNumIteration, relTol, C_SC_ini, compute_sigma_versionToo);
                for cm = 1:obj.i_SC
                    tmpv = K_mu_E_nu_lambds_All_but_Diff{cm};
                    for j = 1:5
                        K_mu_E_nu_lambdas_AllMethods{cm}(ci,j) = tmpv(j);
                    end
                end
            end
        end

        function [cIs, K_mu_E_nu_lambdas_AllMethods, namesOut] = Print_Plot_All_Cs(obj, delcI, cIMax, maxNumIteration, relTol, compute_sigma_versionToo, rootName, b_print, b_plot)
            if ((nargin < 3) || (cIMax < 0))
                cIMax = 1.0;
            end
            if ((nargin < 2) || (delcI < 0))
                delcI = 0.01 * cIMax;
            end
            if ((nargin < 4) || (maxNumIteration < 0))
                maxNumIteration = 100;
            end
            if ((nargin < 5) || (relTol < 0))
                relTol = 1e-3;
            end
            if (nargin < 6)
                compute_sigma_versionToo = 0;
            end
            if (nargin < 7)
                rootName = 'none';
            end
            if (nargin < 8)
                b_print = 1;
            end
            if (nargin < 9)
                b_plot = 1;
            end
            if (strcmp(rootName, 'none') == 0)
                [status,msg,msgID] = mkdir(rootName);
                rootName = [rootName, '/'];
            end
            fieldNames = {'K', 'mu', 'E', 'nu', 'lambda'};
            fieldNames_Latex = {'K', '\mu', 'E', '\nu', '\lambda'};
            num_fields = 5;
            [cIs, K_mu_E_nu_lambdas_AllMethods, namesOut] = Compute_All_Methods_All_cIs(obj, delcI, cIMax, maxNumIteration, relTol, compute_sigma_versionToo);

            num_cIs = length(cIs);
            % printing:
            if (b_print)
                for cm = 1:obj.i_max
                    dat = K_mu_E_nu_lambdas_AllMethods{cm};
                    if (isempty(dat))
                        break;
                    end
                    methodName = obj.names{cm};
                    fn = [rootName, methodName, '.txt'];
                    fid = fopen(fn, 'w');
                    for ci = 1:num_cIs
                        vals = dat(ci,:);
                        fprintf(fid, '%g\t%g\t%g\t%g\t%g\t%g\n', cIs(ci), vals(1), vals(2), vals(3), vals(4), vals(5));
                    end
                    fclose(fid);
                end
                % writing the summary files
                ci = num_cIs;
                summaryNames = {'_finalValues', '_finalValuesRel2Matrix'};
                for sm = 1:2
                    fn = [rootName, summaryNames{sm}, '.txt'];
                    fid = fopen(fn, 'w');
                    fprintf(fid, 'methodName\tK\tmu\tE\tnu\tlambda\n');
                    for cm = 1:obj.i_max
                        methodName = obj.names{cm};
                        fprintf(fid, '%s\t', methodName);
                        dat = K_mu_E_nu_lambdas_AllMethods{cm};
                        if (isempty(dat))
                            fprintf(fid, 'NaN\tNaN\tNaN\tNaN\tNaN\n');
                            continue;
                        end
                        vals = dat(ci,:);
                        if (sm == 2)
                            for j = 1:5
                                vals(j) = vals(j) / obj.K_mu_E_nu_lambdas(1, j);
                            end
                        end
                        fprintf(fid, '%g\t%g\t%g\t%g\t%g\n', vals(1), vals(2), vals(3), vals(4), vals(5));
                    end
                    fclose(fid);
                end
                cmMax = cm;
                if (compute_sigma_versionToo == 0)
                    cmMax = min(cmMax, obj.i_Diff_epsilon);
                end
                % plotting
                if (b_plot)
                    labsz = 22;
                    legsz = 14;
                    lineStyles = cell(cmMax, 1);
                    lineColors = cell(cmMax, 1);
                    lws = 2 * ones(cmMax, 1);
                    for cm = 1:cmMax
                        lineStyles{cm} = '-';
                        lineColors{cm} = [0 0 0]; % black
                    end
                    ls_D_m = '--';
        
                    dark_gray = [0.5	0.5	0.5];
                    blue = [0	0	1];
                    red = [1	0	0];
                    teal = [0	1	1];
                    dark_blue = [0	0	0.5];
                    orange = [1	102/255	0];
                    green = [0 135/255 0];
                    yellow = [1	1	0];
                    purple = [0.5	0	1];
                    brown = [0.5	0.25	0];
        
                    % dark_gray, blue, red, teal, dark_blue, orange, green, yellow, purple, brown
        
                    leg = cell(cmMax, 1);
                    lineStyles{obj.i_Reuss} = ls_D_m;
                    leg{obj.i_Reuss} = 'Reuss';
                    leg{obj.i_Voigt} = 'Voigt';
                    lws(obj.i_Reuss) = 3;
                    lws(obj.i_Voigt) = 3;
        
                    lineStyles{obj.i_HS_minus} = ls_D_m;
                    lineColors{obj.i_HS_minus} = dark_gray;
                    lineColors{obj.i_HS_plus} = dark_gray;
                    leg{obj.i_HS_minus} = 'HS-';
                    leg{obj.i_HS_plus} = 'HS+';
                    lws(obj.i_HS_minus) = 3;
                    lws(obj.i_HS_plus) = 3;
        
                    lineStyles{obj.i_DD_sigma} = ls_D_m;
                    lineColors{obj.i_DD_sigma} = blue;
                    lineColors{obj.i_DD_epsilon} = blue;
                    leg{obj.i_DD_sigma} = 'DD: $$\Sigma$$';
                    leg{obj.i_DD_epsilon} = 'DD: $$E$$';
        
                    lineStyles{obj.i_MT_epsilon} = ':';
                    lineColors{obj.i_MT_epsilon} = red;
                    leg{obj.i_MT_epsilon} = 'MT';
                    lws(obj.i_MT_epsilon) = 1.5;
        
                    lineColors{obj.i_SC} = green;
                    leg{obj.i_SC} = 'SC';
        
                    if (cmMax >= obj.i_Diff_epsilon)
                        lineColors{obj.i_Diff_epsilon} = brown;
                        leg{obj.i_Diff_epsilon} = 'Diff';
                    end
                    if (cmMax >= obj.i_Diff_sigma)
                        lineStyles{obj.i_Diff_sigma} = ls_D_m;
                        lineColors{obj.i_Diff_sigma} = brown;
                        leg{obj.i_Diff_sigma} = 'Diff: $$\Sigma$$';
                    end
                    if (cmMax >= obj.i_MT_sigma)
                        lineStyles{obj.i_MT_sigma} = '-.';
                        lineColors{obj.i_MT_sigma} = red;
                        leg{obj.i_MT_sigma} = 'MT: $$\Sigma$$';
                    end
                    for fii = 1:num_fields
                        fn = fieldNames{fii};
                        fnLatex = fieldNames_Latex{fii};
                        figure(fii);
        
                        for cm = 1:cmMax
                            y = K_mu_E_nu_lambdas_AllMethods{cm}(:, fii);
                            lw = lws(cm);
                            ls = lineStyles{cm};
                            lc = lineColors{cm};
                            if (ls == ':')
                                lw = 1.5;
                            end
                            plot(cIs, y, 'LineWidth', lw, 'Color', lc, 'LineStyle', ls);
                            hold on;
                        end
        
                        xh = get(gca, 'XLabel');
                        set(xh, 'String', '$$ c_I $$', 'FontSize', labsz, 'VerticalAlignment','Top', 'Interpreter', 'latex');
                        yh = get(gca, 'YLabel');
                        set(yh, 'String', ['$$ ', fnLatex,' $$'], 'FontSize', labsz, 'VerticalAlignment','Bottom', 'Interpreter', 'latex');
                        legend(leg, 'FontSize', legsz, 'Interpreter', 'latex', 'Location', 'best');
                        legend('boxoff');
                        if (fii ~= 4)
                            ym = min(obj.K_mu_E_nu_lambdas(:, fii));
                            yM = max(obj.K_mu_E_nu_lambdas(:, fii));
                            if ((yM/ym - 1.0) > 0.01)
                                yle = get(gca, 'ylim');
                                ym = max(ym, yle(1));
                                yM = min(yM, yle(2));
                                ylim([ym, yM]);
                            end
                        end
        
                        fne = [rootName, fn, '.png'];
                        print('-dpng', fne);
                        fne = [rootName, fn, '.fig'];
                        savefig(fne);
                    end
                end
            end
        end

        function [cIs, K_mu_E_nu_lambdas_AllMethods, namesOut, objout] = MAIN_Print_Plot_All_Cs_from_config(obj, configName_woExt)
            if (nargin < 2)
                configName_woExt = 'baseTest';
            end
            configName = [configName_woExt, '.txt'];
            fid = fopen(configName, 'r');
            if (fid < 0)
                fprintf(1, 'Cannot open config file: %s', configName);
                pause;
            end
            buf = 'none';
            while (strcmp(buf, 'START') == 0)
                buf = fscanf(fid, '%s', 1);
            end
            buf = fscanf(fid, '%s', 1);
            cIMax = fscanf(fid, '%f', 1);
            buf = fscanf(fid, '%s', 1);
            delcI = fscanf(fid, '%f', 1);
            buf = fscanf(fid, '%s', 1);
            maxNumIteration = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            relTol = fscanf(fid, '%f', 1);            
            buf = fscanf(fid, '%s', 1);
            compute_sigma_versionToo = fscanf(fid, '%d', 1);
            obj = Read_IsoMeanField(obj, fid);
            [cIs, K_mu_E_nu_lambdas_AllMethods, namesOut] = Print_Plot_All_Cs(obj, delcI, cIMax, maxNumIteration, relTol, compute_sigma_versionToo, configName_woExt);
            objout = obj;
        end

        function [K_mu_E_nu_lambdas_AllMethods, objout] = Compute_PropLastcI(obj, ...
                cI, parasM, parasI, mode, input_type)
            if (nargin < 2)
                cI = 0.23905;
            end
            if (nargin < 3)
                parasM = [45.50, 0.18]; % SiC
            end
            if (nargin < 4)
                parasI = [455.0, 0.18]; % SiC
            end
            if (nargin < 5)
                mode = 21;
            end
            if (nargin < 6)
                input_type = 1;
            end
            cIMax = cI;
            delcI = cIMax / 100;
            maxNumIteration = -1; %100;
            relTol = -1; %1e-3;
            compute_sigma_versionToo = 0;
            obj.mode = mode;
            obj.input_type = input_type;
            obj.matrix_pars = parasM;
            obj.num_inc = 1;
            obj.xis = ones(obj.num_inc, 1);
            obj.inclusion_paras = zeros(obj.num_inc, 2);
            for i = 1:obj.num_inc
                obj.inclusion_paras(i, 1) = parasI(i);
            end
            obj = Initialize(obj);
            configName_woExt = 'none';
            obj.verbose = 0;
            [cIs, K_mu_E_nu_lambdas_AllMethodsTmp, namesOut] = Print_Plot_All_Cs(obj, delcI, cIMax, maxNumIteration, relTol, compute_sigma_versionToo, configName_woExt, 0, 0);
            [lastInd, tmpn] = size(K_mu_E_nu_lambdas_AllMethodsTmp{1});
            numModels = min(length(K_mu_E_nu_lambdas_AllMethodsTmp), 9);
            for mi = 1:numModels
                tmpV = K_mu_E_nu_lambdas_AllMethodsTmp{mi}(lastInd,:);
                K_mu_E_nu_lambdas_AllMethods{mi} = tmpV;
            end
            objout = obj;
        end
    end
end
