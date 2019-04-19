classdef VariationalGaussianMixtureModel < handle
    properties (SetAccess = private)
        n_component; % number of mixture components, K
        sample_size;
        ndim; % D
        
        alpha0;
        m0;
        W0;
        nu0;
        beta0;
        component_size;
        alpha;
        beta;
        m;
        W;
        nu;
        r;
        
        ln_lambda;
        ln_pi;
        S;
        Xm;
    end
    
    methods (Access = public)
        function obj = VariationalGaussianMixtureModel(n_component, alpha0)
            if ~exist('n_component', 'var')
                n_component = 10;
            end
            
            if ~exist('alpha0', 'var')
                alpha0 = 1.0;
            end
            
            obj.n_component = n_component;
            obj.alpha0 = alpha0;
        end
        
        function params = get_params(obj)
            params = {obj.alpha, obj.beta, obj.m, obj.W, obj.nu};
        end

        function fit(obj, X, iter_max)
            if ~exist('iter_max', 'var')
                iter_max = 100;
            end
            
            if isempty(obj.r)
                obj.init_params(X);
            end
            
            for i = 1:iter_max
                array = obj.get_params();
                params_old = [];
                for j = 1:length(array)
                    params_old = [params_old; array{j}(:)];
                end
                
                obj.e_like_step(X);
                e_lower_bound = obj.variational_lower_bound();
                obj.m_like_step(X);
                m_lower_bound = obj.variational_lower_bound();
                
                array = obj.get_params();
                params_new = [];
                for j = 1:length(array)
                    params_new = [params_new; array{j}(:)];
                end
                
                rtol=1e-05;
                atol=1e-08;
                if all(abs(params_old - params_new) <= atol + rtol * abs(params_new))
                    break;
                elseif i == iter_max
                    disp('parameters may not have converged');
                else
                    disp(strcat(num2str(i), '-th iteration:',...
                        'lower bound E...', num2str(e_lower_bound),...
                        ' and lower bound M...', num2str(m_lower_bound)));
                end
            end
        end
        
        function label = classify(obj, X)
            [~, label] = max(obj.e_like_step(X), [], 2);
        end

        function [probs, components] = predict_dist(obj, X)
            % PRML (10.81)
            probs = zeros([size(X, 1), 1]);
            components = zeros([size(X, 1), obj.n_component]);
            sum_alpha = sum(obj.alpha);
            for k = 1:obj.n_component
                precision = (obj.nu(k) + 1 - obj.ndim) * obj.beta(k)...
                    * obj.W(:,:,k) / (1 + obj.beta(k));
                tdist = Math.Statistics.StudentsTDistribution(obj.m(:, k),...
                    precision, obj.nu(k) + 1 - obj.ndim);
                components(:, k) = obj.alpha(k) * tdist.pdf(X) / sum_alpha;
                probs = probs + components(:, k);
            end
            %probs = probs / sum(obj.alpha);
        end
    end
    
    methods(Access = private)
        function init_params(obj, X)
            rng(0, 'twister');
            
            [obj.sample_size, obj.ndim] = size(X); % N * D
            obj.alpha0 = ones(obj.n_component, 1) * obj.alpha0; % K * 1
            obj.m0 = zeros(obj.ndim, 1); % D * 1
            obj.W0 = eye(obj.ndim); % D * D
            obj.nu0 = obj.ndim;
            obj.beta0 = 1.0;

            obj.component_size = obj.sample_size / obj.n_component + zeros(obj.n_component, 1); % K * 1
            obj.alpha = obj.alpha0 + obj.component_size; % K * 1
            obj.beta = obj.beta0 + obj.component_size; % K * 1
            indices = randperm(obj.sample_size, obj.n_component);
            obj.m = X(indices, :)'; % D * K
            obj.W = repmat(obj.W0, 1, 1, obj.n_component); % D * D * K
            obj.nu = obj.nu0 + obj.component_size; % K * 1
            
            obj.ln_lambda = zeros(obj.n_component, 1); % K * 1
            obj.ln_pi = zeros(obj.n_component, 1); % K * 1
            obj.S = zeros(obj.ndim, obj.ndim, obj.n_component); % D * D * K
            obj.Xm = zeros(obj.ndim, obj.n_component); % D * K
            
            rng('shuffle');
        end
        
        function r = e_like_step(obj, X)
            sample_size = size(X, 1);
            
            % PRML (10.64)
            expect_value = zeros(sample_size, obj.n_component); % N * K
            for k = 1:obj.n_component
                % These codes are equivalent to the following commented-out
                % codes.
                d = X' - obj.m(:, k); % D * step_size
                expect_value(:, k) = obj.ndim / obj.beta(k) + obj.nu(k) * sum((d' * obj.W(:, :, k)) .* d', 2);
                
                %{
                for n = 1:sample_size
                    d = X(n, :)' - obj.m(:, k); % D * 1
                    expect_value(n, k) = obj.ndim / obj.beta(k) + obj.nu(k) * d' * obj.W(:, :, k) * d;
                end
                %}
            end
            
            % PRML (10.65)
            for k = 1:obj.n_component
                i = 1:obj.ndim;
                obj.ln_lambda(k) = sum(psi((obj.nu(k) + 1 - i) / 2)) + obj.ndim * log(2) + log(det(obj.W(:, :, k)));
            end
            
            % PRML (10.66)
            obj.ln_pi = psi(obj.alpha) - psi(sum(obj.alpha)); % K * 1
            
            % PRML (10.46)0
            ln_rho = obj.ln_pi' + 0.5 * obj.ln_lambda' - 0.5 * obj.ndim * log(2*pi) - 0.5 * expect_value;
            
            % PRML (10.49)
            rho = exp(ln_rho);
            obj.r = rho ./ sum(rho, 2); % N * K
            obj.r(isnan(obj.r) | isinf(obj.r) | obj.r == 0) = sqrt(eps);
            r = obj.r;
        end
        
        function m_like_step(obj, X)
            % PRML (10.51)
            obj.component_size = (sum(obj.r, 1))'; % K * 1
            
            % PRML (10.52)
            obj.Xm = (obj.r' * X ./ obj.component_size)'; % D * K
            
            % PRML (10.53)
            obj.S = zeros(obj.ndim, obj.ndim, obj.n_component); % D * D * K
            for k = 1:obj.n_component
                % These codes are equivalent to the following commented-out
                % codes.
                d = X' - obj.Xm(:, k);
                obj.S(:, :, k) = obj.r(:, k)' .* d * d';
                
                %{
                for n = 1:obj.sample_size
                    d = X(n, :)' - obj.Xm(:, k);
                    obj.S(:, :, k) = obj.S(:, :, k) + obj.r(n, k) * (d * d');
                end
                %}
                                
                obj.S(:, :, k) = obj.S(:, :, k) / obj.component_size(k);
            end
            
            % PRML (10.58)
            obj.alpha = obj.alpha0 + obj.component_size;
            
            % PRML (10.60 - 10.63)
            obj.beta = obj.beta0 + obj.component_size;
            obj.m = (obj.beta0 * obj.m0 + obj.component_size' .* obj.Xm) ./ obj.beta';
            
            for k = 1:obj.n_component
                d = obj.Xm(:, k) - obj.m0;
                obj.W(:, :, k) = inv(inv(obj.W0) + (obj.component_size(k) .* obj.S(:, :, k))...
                    + (obj.beta0 * obj.component_size(k) * (d * d')) / (obj.beta0 + obj.component_size(k)));
            end
            
            obj.nu = obj.nu0 + obj.component_size;
        end
        
        function value = variational_lower_bound(obj)
            % PRML (10.71)
            exp_ln_p_X = 0;
            for k = 1:obj.n_component
                d = obj.Xm(:, k) - obj.m(:, k);
                exp_ln_p_X = exp_ln_p_X...
                    + obj.component_size(k)...
                    * (obj.ln_lambda(k) - obj.ndim / obj.beta(k) - obj.nu(k)...
                    * trace(obj.S(:,:,k) * obj.W(:,:,k)) - obj.nu(k) * d'...
                    * obj.W(:,:,k) * d - obj.ndim * log(2*pi));
            end
            exp_ln_p_X = exp_ln_p_X / 2;
            
            % PRML (10.72)
            exp_ln_p_Z = sum(sum(obj.r .* obj.ln_pi'));
            
            % PRML (10.73)
            ln_C = gammaln(sum(obj.alpha0)) - sum(gammaln(obj.alpha0));
            exp_ln_p_pi = ln_C + obj.alpha0(1) * sum(obj.ln_pi);
            
            % PRML (10.74)
            first_term = 0;
            last_term = 0;
            for k = 1:obj.n_component
                d = obj.m(:, k) - obj.m0;
                first_term = first_term...
                    + obj.n_component * log(obj.beta0 / (2 * pi))...
                    + obj.ln_lambda(k)...
                    - obj.ndim * obj.beta0 / obj.beta(k)...
                    - obj.beta0 * obj.nu(k) * d' * obj.W(:,:,k) * d;
                last_term = last_term...
                    + obj.nu(k) * trace(obj.W0 \ obj.W(:,:,k));
            end
            first_term = first_term / 2;
            last_term = last_term / 2;
            ln_B = (-obj.nu0 / 2) * log(det(obj.W0))...
                - obj.nu0 * obj.ndim * log(2) / 2 ...
                - obj.ndim * (obj.ndim - 1) * log(pi) / 4 ...
                - sum(gammaln((obj.nu0 + 1 - 1:obj.ndim) / 2));
            exp_ln_p_mu_lambda = first_term...
                + obj.n_component * ln_B...
                + 0.5 * (obj.nu0 - obj.ndim - 1) * sum(obj.ln_lambda)...
                - last_term;
            
            % PRML (10.75)
            exp_ln_q_Z = sum(sum(obj.r .* log(obj.r)));
            
            % PRML (10.76)
            ln_C = gammaln(sum(obj.alpha)) - sum(gammaln(obj.alpha));
            exp_ln_q_pi = sum((obj.alpha - 1) .* obj.ln_pi + ln_C);
            
            % PRML (10.77)
            ln_B = zeros(obj.n_component, 1);
            H_lambda = zeros(obj.n_component, 1);
            for k = 1:obj.n_component
                ln_B(k) = (-obj.nu(k) / 2) * log(det(obj.W(:,:,k)))...
                    - obj.nu(k) * obj.ndim * log(2) / 2 ...
                    - obj.ndim * (obj.ndim - 1) * log(pi) / 4 ...
                    - sum(gammaln((obj.nu(k) + 1 - 1:obj.ndim) / 2));
                H_lambda(k) = -ln_B(k) - (obj.nu(k) - obj.ndim - 1) * obj.ln_lambda(k) / 2 ...
                    + obj.nu(k) * obj.ndim / 2;
            end
            exp_ln_q_mu_lambda = sum(obj.ln_lambda + obj.ndim .* log(obj.beta / 2 / pi)...
                - obj.ndim - 2 .* H_lambda) / 2;
            
            value = exp_ln_p_X + exp_ln_p_Z + exp_ln_p_pi + exp_ln_p_mu_lambda...
                - exp_ln_q_Z - exp_ln_q_pi - exp_ln_q_mu_lambda;
        end
    end
end