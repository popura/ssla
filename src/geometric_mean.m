function y = geometric_mean(x)
x(x == 0) = sqrt(eps);
y = exp(mean(log(x(:))));
end