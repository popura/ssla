function output_image = replace_colors(input_image, new_luminance, old_luminance)
normalizer = new_luminance ./ old_luminance;
normalizer(old_luminance == 0) = 0;
s = size(input_image);

if length(s) == 3
    output_image(:,:,1) = normalizer .* input_image(:,:,1);
    output_image(:,:,2) = normalizer .* input_image(:,:,2);
    output_image(:,:,3) = normalizer .* input_image(:,:,3);
elseif length(s) == 2
    output_image = normalizer .* input_image;
else
    output_image = [];
end

output_image(isnan(output_image) | isinf(output_image)) = 0;
end