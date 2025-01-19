function res = integralL2(f, h)
%% L2 integral

% get stepsize h automatically
if nargin == 1
    n = size(f, 1);
    h = 1 / n;
end

% L2 integral
res = reshape(h * sum(f, 1), [], 1);

end
