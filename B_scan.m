tic
Exercise_5_3

%------------------------------------------------------------------
%CALCULATE B SCAN TIME
%------------------------------------------------------------------         
time_B_scan = (2.*mz)./velocity;
%index = round((time_B_scan ./ max_time) .* fft_points); 

d_txrx_repmat = repmat(time_B_scan, [1 1 1 transducer_elements]);

%------------------------------------------------------------------
%FIND TIME INDEXES
%------------------------------------------------------------------ 
index = round((d_txrx_repmat ./ max_time) .* fft_points); 
index(index==0) = 1;
index(index>fft_points) = fft_points;

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
                allow = abs(source_x_positions(kk) - mx(1,jj,1));
                if D/2 >= allow
                    I(ii,jj,kk,ll) = time_domain(kk,ll,index(ii,jj,kk,ll));
                end
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
title('TFM Intensity Plot')
xticklabels = -10:2:10;
xticks = linspace(1, size(I, 2), numel(xticklabels));
yticklabels = 0:2:grid_size*1000;
yticks = linspace(1, size(I, 2), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('X Position (mm)'); 
ylabel('Z Position (mm)');
toc