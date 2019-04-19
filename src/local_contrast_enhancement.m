function enhanced_images = local_contrast_enhancement(images)
enhanced_images = cell(size(images));
for i = 1:length(images)
    enhanced_images{i} = retinex_filter(images{i});
end
end