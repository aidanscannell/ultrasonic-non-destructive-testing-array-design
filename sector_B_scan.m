tic
Exercise_5_3

%------------------------------------------------------------------
%CALCULATE DELAYS
%------------------------------------------------------------------   
max_time = 2*distance_to_wall/velocity;
r = sqrt(mx.^2 + mz.^2);
r_4D = repmat(r, [1,1,1,transducer_elements]);
angle_sector = atan(mx ./ mz);
angle_sector = repmat(angle_sector, [1,1,1,transducer_elements]);
ms_i = repmat(ms, [1,1,1,transducer_elements]);
ms_j = permute(ms_i, [1 2 4 3]);
x_tx_sin = abs(ms_i .* sin(angle_sector));
x_rx_sin = abs(ms_j .* sin(angle_ sector));

d_txrx_new = (2 .* r_4D + x_tx_sin + x_rx_sin) ./ velocity;

%------------------------------------------------------------------
%FIND TIME INDEXES
%------------------------------------------------------------------ 
% index = round((d_txrx_new ./ max_time) .* fft_points); 
% index(index==0) = 1;
% index(index>fft_points) = fft_points;
% index = round(index);

index = round((d_txrx_new ./ max_time) .* length(time)); 
index(index==0) = 1;
index(index>length(time)) = length(time);
%------------------------------------------------------------------
%APERTURE
%------------------------------------------------------------------ 
D = transducer_pitch * 5;

%------------------------------------------------------------------
%INITIALISE INTENSITY MATRIX
%------------------------------------------------------------------ 
I = zeros(length(x),length(z),length(source_x_positions),length(source_x_positions));

%------------------------------------------------------------------
%LOOP FOR INTENSITY
%------------------------------------------------------------------ 
for ii = 1 : length(x)
    for jj = 1 : length(z)
        for kk = 1 : length(source_x_positions)
            for ll = 1 : length(source_x_positions)
                I(ii,jj,kk,ll) = time_domain(kk,ll,index(ii,jj,kk,ll));
            end
        end
    end
end

I = sum(I,3);
I = sum(I,4);
%I = 20 * log10(I/max(max(I)));
%------------------------------------------------------------------
%PLOT INTENSITY
%------------------------------------------------------------------ 
hold off
imagesc(abs(I))
hold on
scatter(focal_points_index(:,1),focal_points_index(:,2), 20, 'filled', 'r');
scatter(100,back_wall_index,20, 'filled', 'r');
title('TFM Intensity Plot')
xticklabels = -10:1:10;
xticks = linspace(1, size(I, 2), numel(xticklabels));
yticklabels = 0:1:22;
yticks = linspace(1, size(I, 2), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('X Position (mm)'); 
ylabel('Z Position (mm)');
toc