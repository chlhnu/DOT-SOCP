classdef ModelHandle < handle
    % Handle wrapper to avoid struct copies for model data
    properties
        nx
        ny
        nt
        grad
        c
        normc
        normd
        rho0
        rho1
        cScale
        dScale
        D
        E
    end
    methods
        function obj = ModelHandle(s)
            if nargin > 0
                fns = fieldnames(s);
                for k = 1:numel(fns)
                    fn = fns{k};
                    if ~isprop(obj, fn)
                        addprop(obj, fn);
                    end
                    obj.(fn) = s.(fn);
                end
            end
        end
    end
end
