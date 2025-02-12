function B = custom_padarray(A, padsize, padval, direction)
    try
        B = padarray(A, padsize, padval, direction);
    catch
        if nargin < 3
            padval = 0;
        end
        if nargin < 4
            direction = 'both';
        end

        sizeA = size(A);
        numDims = numel(sizeA);

        if numel(padsize) < numDims
            padsize = [padsize, zeros(1, numDims - numel(padsize))];
        end

        if strcmp(direction, 'pre')
            padsize = [padsize; zeros(1, numDims)];
        elseif strcmp(direction, 'post')
            padsize = [zeros(1, numDims); padsize];
        elseif strcmp(direction, 'both')
            padsize = [ceil(padsize / 2); floor(padsize / 2)];
        else
            error('Invalid direction. Use ''pre'', ''post'', or ''both''.');
        end

        sizeB = sizeA + sum(padsize, 1);
        B = padval * ones(sizeB, class(A));

        idx = cell(1, numDims);
        for d = 1:numDims
            idx{d} = (1:sizeA(d)) + padsize(1, d);
        end

        B(idx{:}) = A;
    end
end