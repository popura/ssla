classdef StudentsTDistribution < handle
    properties(SetAccess = private)
        mu; % mean (D-dimention vector)
        lambda; % precision matrix (D * D matrix, inverse matrix of covariance matrix)
        nu; % the number of degrees of freedom (positive real number)
        ndim; % the number of dimentions
    end
    
    methods(Access = public)
        function obj = StudentsTDistribution(mu, lambda, nu)
            if isvector(mu)
                if ~iscolumn(mu)
                    obj.mu = mu';
                else
                    obj.mu = mu;
                end
            else
                error('mu should be a scalar or vector');
            end
            
            obj.ndim = size(obj.mu, 1);
            
            if issymmetric(lambda) && size(lambda, 1) == obj.ndim
                obj.lambda = lambda;
            else
                error('lambda should be a D * D symmetric matrix');
            end
            
            if isscalar(nu) && nu > 0
                obj.nu = nu;
            else
                error('nu should be a positive scalar value');
            end
        end
        
        function y = pdf(obj, x)
            % In the case when mu is a vector:
            % x is a N * D matrix, where N is the number of samples.
            % In the case when mu and lambda are scalar values:
            % any matrix or array can be used as x.
            % In this case, this function construes each elements of x as a
            % sample.
            if ~ismatrix(x)
                error('x should be a matrix');
            end
            
            if ~isscalar(obj.lambda) && size(x, 2) ~= obj.ndim
                error('x should be a matrix which includes samples as row vectors');
            end
            
            if obj.ndim == 1
                gamma_ratio = exp(gammaln((obj.nu + 1) / 2) - gammaln(obj.nu / 2));
                sqdelta = obj.lambda * ((x - obj.mu).^2);
                y = gamma_ratio * sqrt(obj.lambda / (obj.nu * pi))...
                    * (1 + sqdelta / obj.nu).^(-(obj.nu + 1) / 2);
            else
                x_t = x';
                gamma_ratio = exp(gammaln((obj.nu + 1) / 2) - gammaln(obj.nu / 2));
                sqdelta = diag((x_t - obj.mu)' * obj.lambda * (x_t - obj.mu));
                y = gamma_ratio * sqrt(det(obj.lambda)) / (obj.nu * pi)^(obj.ndim / 2)...
                    * (1 + sqdelta / obj.nu).^(-(obj.nu + obj.ndim) / 2);
            end
        end
        
        function y = cdf(obj, x)
        end
        
        function y = mean(obj)
            y = obj.mu;
        end
        
        function y = variance(obj)
            y = obj.nu / (obj.nu - 2) * inv(obj.lambda);
        end
        
        function y = mode(obj)
            y = obj.mu;
        end
    end
end