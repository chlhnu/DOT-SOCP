function [sigma, factor] = adjust_lagrangianParam(sigma, xi, updateRule, bound)
% Adjust Lagrangian param
% Input
%   sigma: Lagrangian parameter
%   xi: ratio between primal and dual residuals
%   updateRule: "nx2" matrix, 1st column contains xi values and 2nd column contains corresponding factors
%       for example: updateRule = [5,   1.10;
%                                  50,  1.65;
%                                  500, 2.20];
% Output
%   sigma: updated Lagrangian param
%   factor: multiplication factor for updating (factor == 1 means no update)

if nargin == 3
    lowerBound = 1e-3;
    upperBound = 1e3;
else
    lowerBound = bound(1);
    upperBound = bound(2);
end

% If show sigma updating message
if ~exist("printYes", "var")
    printYes = false;
end

% Get sigma-factor accronding to xi
if xi >= 1
    factor = get_factor(xi, updateRule);
elseif xi < 1
    factor = 1 / get_factor(1/xi, updateRule);
end

% Update sigma with safeguard
if (factor ~= 1)
    sigmaOld = sigma;
    sigma    = max(min(sigma * factor, upperBound), lowerBound);
    factor   = sigma / sigmaOld;
end

% Print sigma message
if (printYes) && (factor ~= 1)
    disp("newSigma = sigma * " + num2str(factor) + " = " + num2str(sigma));
end

end

%% calc factor which updating sigma
function factor = get_factor(xi, updateRule)
    factor = 1;

    for i = 1 : size(updateRule, 1)
        if xi >= updateRule(i, 1)
            factor = updateRule(i, 2);
            continue;
        else
            break;
        end
    end
end