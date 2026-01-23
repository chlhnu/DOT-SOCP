function res = integralL2(f, h)
%% L2 integral

if nargin == 1
    n = size(f, 1);
    h = 1 / n;
end

res = reshape(h * sum(f, 1), [], 1);

end
 