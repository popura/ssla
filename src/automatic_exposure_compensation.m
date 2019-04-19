function adjusted_images = automatic_exposure_compensation(images)
max_region_number = length(images);
            
sorted_images = sort_images(images);
image_number = length(images);
med = ceil((1 + image_number) / 2);
image = sorted_images{med};

luminance = calculate_luminance(image);
sorted_luminance = sort(luminance(:), 'ascend');

range = sorted_luminance(end) - sorted_luminance(1);
tmp_endpoints = (range / max_region_number) * (0:max_region_number);
endpoints = tmp_endpoints + sorted_luminance(1);

regions = cell(max_region_number, 1);
for i = 1:max_region_number
    j = max_region_number - i + 1;
    if i == 1
        regions{i} = luminance >= endpoints(j);
        continue;
    end
    regions{i} = (luminance >= endpoints(j)) & (luminance < endpoints(j+1));
end

regions{med} = true(size(luminance));

scale_factors = zeros(1, image_number);
adjusted_images = cell(size(images));
for i = 1:image_number
    luminance = calculate_luminance(images{i});
    geo_mean = geometric_mean(luminance(regions{i}));
    scale_factors(i) = 0.18 / geo_mean;
    adjusted_images{i} = scale_factors(i) * images{i};
end

end