classdef VarHandle < handle
    % Handle wrapper to avoid struct copies for large iterative variables
    properties
        phi
        q
        z
        alpha
        beta
        cScale
        dScale
        D
        E
        E2
        qInd
        name
        time
    end
    methods
        function obj = VarHandle(s)
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
