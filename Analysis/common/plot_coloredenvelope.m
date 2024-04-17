
function [patch,h2,h] = plot_coloredenvelope(x,data, ref)

        options.handle     = figure;
        options.color_area = [128 193 219]./255;   
        options.color_line = [ 52 148 186]./255;
        options.alpha      = 0.5;
        options.line_width = 1;

    options.x_axis = x;
     
    % Plotting the result
    h = figure(options.handle);
    x_vector = [options.x_axis, fliplr(options.x_axis)];
    patch = fill(x_vector, [max(data,[],2)',fliplr(min(data,[],2)')], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    h2 = plot(options.x_axis, ref, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    hold off;
    
end