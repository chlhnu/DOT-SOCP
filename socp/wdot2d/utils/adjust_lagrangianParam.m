function [sigma, factor] = adjust_lagrangianParam(sigma, xi, updateRule, bound)
%% Adjust Lagrangian parameter

if nargin == 3
    lowerBound = 1e-3;
    upperBound = 1e3;
else
    lowerBound = bound(1);
    upperBound = bound(2);
end

if ~exist("printYes", "var")
    printYes = false;
end

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

if (printYes) && (factor ~= 1)
    disp("newSigma = sigma * " + num2str(factor) + " = " + num2str(sigma));
end

end

%% Factor for updating sigma
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