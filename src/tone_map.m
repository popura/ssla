function mapped_images = tone_map(images)
mapped_images = cell(size(images));
for i = 1:length(images)
    luminance = calculate_luminance(images{i});
    white_point = max(luminance(:));    
    mapped_luminance = (luminance ./ (1 + luminance)) .* (1 + (luminance / (white_point^2)));
    mapped_images{i} = replace_colors(images{i}, mapped_luminance, luminance);
end
end