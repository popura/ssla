function output_images = sort_images(input_images)
geo_mean = zeros(size(input_images));
for i = 1:length(input_images(:))
    luminance = calculate_luminance(input_images{i});
    geo_mean(i) = geometric_mean(luminance);
end

[~, index] = sort(geo_mean);
output_images = input_images(index);
end