function x = remove_outlier(y,n)
numInputs = nargin;
if numInputs == 1;
    n=3;
end

x = y;
outliers = y - mean(y(:)) > n*std(y(:));
x(outliers) = max(max(x(~outliers)));