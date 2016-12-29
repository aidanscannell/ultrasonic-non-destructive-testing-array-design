tic
Exercise_5_3

%------------------------------------------------------------------
%CALCULATE B SCAN TIME
%------------------------------------------------------------------         
time_B_scan = (2.*mz)./velocity;
d_txrx_repmat = repmat(time_B_scan, [1 1 1 transducer_elements]);

%------------------------------------------------------------------
%FIND TIME INDEXES
%------------------------------------------------------------------ 
index = round((d_txrx_repmat ./ max_time) .* fft_points); 
index(index==0) = 1;
index(index>fft_points) = fft_points;

%------------------------------------------------------------------
%APERTURE
%------------------------------------------------------------------ 
D = transducer_pitch * 5 - spacing;

%------------------------------------------------------------------
%INITIALISE INTENSITY MATRIX
%------------------------------------------------------------------ 
I = zeros(length(x),length(z),length(source_x_positions),length(source_x_positions));

%------------------------------------------------------------------
%LOOP FOR INTENSITY
%------------------------------------------------------------------ 
aa=0;
for ii = 1 : length(x)
    for jj = 1 : length(z)
        for kk = 1 : length(source_x_positions)
            for ll = 1 : length(source_x_positions)
                allow = abs(source_x_positions(kk) - mx(1,jj,1));
                if D/2 >= allow
                    I(ii,jj,kk,ll) = time_domain(kk,ll,index(ii,jj,kk,ll));
                    aa = aa + 1;
                end
            end
        end
    end
end
toc
I = sum(I,3);
I = sum(I,4);

%------------------------------------------------------------------
%DETERMINE API
%------------------------------------------------------------------
I_API = abs(I/max(max(I)));

API_cutoff = I_API;
API_cutoff(API_cutoff <= 0.1) = 0;

maxima = imregionalmax(API_cutoff(1:grid_pts-20,10:grid_pts-10));
[row col] = find(maxima);

figure(length(row)+1)
hold off
imagesc(maxima)
hold on
title('Focused B Scan API Maxima')
xlabel('x Position')
ylabel('z Position')
hold off
tol = 3;
interp_val = 3;
pixel_size = (grid_size/grid_pts)^2 ./ (interp_val * 2^(interp_val-1) - 1)^2;

for ii = 1 : length(row)
    API = I_API(row(ii)-tol:row(ii)+tol,col(ii)-tol:col(ii)+tol);
    interp_I = interp2(API,interp_val);
    interp_I = interp_I/max(max(interp_I));
    interp_I(interp_I < 0.5) = 0;
    %imagesc(abs(interp_I))
    %API(API>=0.5) = 1;
    area_data = regionprops('struct',im2bw(interp_I),'Centroid','Area');
    data(ii,1) = row(ii);
    data(ii,2) = col(ii);
    data(ii,3) = area_data(1).Area * pixel_size;
    data(ii,4) = data(ii,3) / wavelength^2;
    
    figure(ii)
    imagesc(abs(interp_I))
end


I_2 = interp2(I,3);
%------------------------------------------------------------------
%PLOT INTENSITY
%------------------------------------------------------------------ 
hold off
figure(length(row)+2)
I_3 = -20*log10(I_2/max(max(I_2)));
I_3(isinf(I_3))=0;
imagesc(abs(I_3))
hold on
%scatter(focal_points_index(:,1),focal_points_index(:,2), 20, 'filled', 'r');
title('Plane B Scan Intensity Plot')
xticklabels = -10:2:10;
xticks = linspace(1, size(I_3, 2), numel(xticklabels));
yticklabels = 0:2:grid_size*1000;
yticks = linspace(1, size(I_3, 2), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('X Position (mm)'); 
ylabel('Z Position (mm)');
