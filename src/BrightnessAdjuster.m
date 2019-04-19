classdef BrightnessAdjuster < handle
    properties (SetAccess = protected)
        standard_brightness;
        compensation_mode;
    end
    
    methods (Access = public)
        function obj = BrightnessAdjuster(standard_brightness)            
            if ~exist('standard_brightness', 'var')
                standard_brightness = 0.18;
            end
                                    
            obj.standard_brightness = standard_brightness;
        end

        function adjusted_images = adjust(obj, images, regions)
            adjusted_images = cell(1, length(regions));
            for i = 1:length(regions)
                geometric_means = zeros(size(images));
                min_difference = inf;
                for j = 1:length(images(:))
                    luminance = calculate_luminance(images{j});
                    geometric_means(j) = geometric_mean(luminance(regions{i}));
                    
                    tmp_difference = (obj.standard_brightness - geometric_means(j))^2;
                    if min_difference > tmp_difference
                        min_difference = tmp_difference;
                        nearest_image_index = j;
                    end
                end
                disp(strcat('nearest image:', num2str(nearest_image_index), '-th image'));

                scale_factor = obj.standard_brightness / geometric_means(nearest_image_index);
                lum_old = calculate_luminance(images{nearest_image_index});
                lum_new = scale_factor * lum_old;
                adjusted_images{i} = replace_colors(images{nearest_image_index}, lum_new, lum_old);
            end
            
            adjusted_images = sort_images(adjusted_images);
        end
    end
end