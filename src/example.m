% read multi-exposure images
image_files = dir(fullfile('..', 'data', '*.JPG'));
input_images = cell(1, length(image_files));
for i = 1:length(image_files)
    tmp_image = im2double(imread(fullfile(image_files(i).folder, image_files(i).name)));
    input_images{i} = tmp_image.^2.2; % inverse gamma correction
end

% local contrast enhancement
enhanced_images = local_contrast_enhancement(input_images);

% scene segmentation
extractor = ROIExtractor('approach2');
regions = extractor.extract(enhanced_images);

% luminance scaling
adjuster = BrightnessAdjuster(0.18);
adjusted_images = adjuster.adjust(enhanced_images, regions);

% tone mapping
mapped_images = tone_map(adjusted_images);

output_images = cell(1, length(image_files));
for i = 1:length(image_files)
    tmp_image = mapped_images{i};
    tmp_image(tmp_image < 0) = 0;
    tmp_image(tmp_image > 1) = 1;
    output_images{i} = tmp_image.^(1/2.2); % gamma correction
    figure();
    imshow(output_images{i});
end

% these output images can be fused by any fusion methods.
% here, simple average is used as an example.
fused_image = zeros(size(output_images{1}));
for i = 1:length(image_files)
    fused_image = fused_image + output_images{i};
end
fused_image = fused_image / length(image_files);
figure();
imshow(fused_image);

% if you use MATLAB R2018b or later,
% blendexposure() function can fuse adjusted_images.