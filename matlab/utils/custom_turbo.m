function cmap = custom_turbo(n)
if nargin < 1
    n = 256;
end

try
    cmap = turbo(n);
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        cmap = alternative_turbo(n);
    else
        rethrow(ME);
    end
end

end

function cmap = alternative_turbo(n)
    colors = [
        0.18995, 0.07176, 0.23217; 
        0.25107, 0.57053, 0.82644; 
        0.97234, 0.90559, 0.10912; 
        0.81337, 0.11745, 0.13174  
    ];

    x = linspace(0, 1, size(colors, 1));
    xi = linspace(0, 1, n);
    cmap = interp1(x, colors, xi, 'linear');
end
