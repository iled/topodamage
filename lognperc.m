function Y = lognperc(X, p)
% Y = lognperc(X, p) returns percentiles of the values in a data vector or
% matrix X for the percentages p in the interval [0, 100].
% A check if performed to remove NaN's, Inf's and non-positive values from X.
% The parameters of the underlying normal distribution are estimated using the
% function lognfit from the Statistics and Machine Learning Toolbox.
% Julio Caineta Nov 7 2017

% clean up data from NaN's, Inf's and non-positive values
X = X(~isnan(X) & ~isinf(X) & X > 0);

% estimate parameters of the underlying normal distribution
parmhat = lognfit(X);
mu = parmhat(1);
sigma = parmhat(2);

Y = exp(norminv(p ./ 100, mu, sigma));


