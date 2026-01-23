function colors = generate_colormap(start_color, end_color, num_colors)
%% Generate a series colors to be used by colormap()
% Input:
%   start_color and end_color both are 1x3 vectors

if exist('num_colors', 'var') && ~isempty(num_colors)
    if ~ (num_colors == round(num_colors) && num_colors >= 2)
        error("Invalid input at position 3.");
    end
end

if isempty(num_colors)
    num_colors = 64;
end

if all((start_color >= 0) .* (start_color <= 255))
    if ~ all((start_color >= 0) .* (start_color <= 1))
        start_color = start_color / 255;
    end
else
    error("Invalid input at position 1.");
end

if all((end_color >= 0) .* (end_color <= 255))
    if ~ all((end_color >= 0) .* (end_color <= 1))
        end_color = end_color / 255;
    end
else
    error("Invalid input at position 2.");
end

colors = [
    linspace(start_color(1), end_color(1), num_colors)', ...
    linspace(start_color(2), end_color(2), num_colors)', ...
    linspace(start_color(3), end_color(3), num_colors)'];

end



