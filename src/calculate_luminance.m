function luminance = calculate_luminance(image)
xyz = rgb2xyz(image, 'ColorSpace', 'linear-rgb', 'WhitePoint', 'D65');
luminance = xyz(:, :, 2);
end