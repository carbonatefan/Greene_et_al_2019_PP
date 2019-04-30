%This script is published alongside Greene et al.,2019 'Early Cenozoic 
%Decoupling of Climate and Carbonate Compensation Depth Trends'. Please 
%cite usage accordingly.

% plot_CCDcontour
%
%   ***********************************************************************
%   *** PLOT P-E CCD DATA *************************************************
%   ***********************************************************************
%
%   plot_peccd(PFILE) plots PE CCD data and takes 1 arguments:
%
%   PFILE [STRING] (e.g. 'MyData.dat')
%   --> the data file to be plotted (4 columns) where column 1 is
%   (paleo)longitude, column 2 is (paleo)latitude, column 3 is
%   (paleo)depth, column 4 is wt% CaCO3.
%
%   ***********************************************************************

%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/03/03: CREATED
%   14/03/04: added further plot output
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES, LOAD DATA ********************** %
% *********************************************************************** %
%
clear all;
close all;
% set CCD data filename
filename = 'file.txt'; %Enter your filename here
% set min,max depth range
D_max = 6000;
D_min = 0;
C_min=0;
C_max=100;
% set window width (depth space) and depth increment (for sliding window)
D_width = 1000;
D_dD = 100;
C_width=10;
C_dC=10;
% set file and date strings for file-saving
str_window = num2str(D_width);
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% define wt% bin mid-points
%%%bin_wtpct = [5 20 40 60 80 95];
%%%bin_dwtpct = [10 20 20 20 20 10];
bin_wtpct = [C_min:C_dC:C_max];
bin_dwtpct = 0.0*bin_wtpct + C_dC;
n_max = length(bin_wtpct);
plot_wtpct = [0:2:100];
plot_wtpct_n = length(plot_wtpct);
% define depth mid-points
bin_D = [(D_min + D_width/2.0):D_dD:(D_max - D_width/2.0)];
bin_D_n=length(bin_D);
plot_D = [D_min:200:D_max];
plot_D_n = length(plot_D);
% set up basic grids
[X Y] = meshgrid(bin_wtpct,bin_D);
Z = NaN*zeros(size(X));
Z_n = Z;
Z_w1 = Z;
% set up interpolated grid
[XI YI] = meshgrid(plot_wtpct,plot_D);
% load CCD data
data = load(filename);
%
% *********************************************************************** %

% *********************************************************************** %
% *** PROCESS DATA ****************************************************** %
% *********************************************************************** %
%
% bin in 2D
for n = 1:n_max,
    % initialize inner loop
    data_Dbin = data(find(abs(data(:,4)-bin_wtpct(n)) <= (bin_dwtpct(n)/2.0)),:); 
    m = 1;
    D_mid = D_min + D_width/2.0;
    while D_mid <= (D_max - D_width/2.0),
        loc_indx = find( abs(data_Dbin(:,3)-D_mid) <= D_width/2.0 );
        Z_n(m,n) = length(loc_indx); 
        Z_w1(m,n) = sum( 1.0 - abs(data_Dbin(loc_indx,3)-D_mid)/(D_width/2.0) ); 
        % update inner loop counters
        m = m+1;
        D_mid = D_mid + D_dD;
    end
end
% normalize binned data in wt% space
% NOTE: assumes triangular weighting function
m = 1;
D_mid = D_min + D_width/2.0;
while D_mid <= (D_max - D_width/2.0),
    loc_count = sum(Z_n(m,:));
    if (loc_count > 0), 
        Z_n(m,:) = Z_n(m,:)/loc_count; 
        Z_w2(m,:) = Z_w1(m,:)/sum(Z_w1(m,:)); 
    else
        Z_n(m,:) = NaN;
        Z_w1(m,:) = NaN;
    end
    m = m+1;
    D_mid = D_mid + D_dD;
end
% interpolate data
for i=1:size(Z_w2,1)
    wtpct_intp(i,:)=interp1(bin_wtpct,Z_w2(i,:),plot_wtpct);
end
ZI_n = interp2(X,Y,Z_n,XI,YI,'linear');
ZI_w1 = interp2(X,Y,Z_w2,XI,YI,'linear');
ZI_w2 = interp2(X,Y,Z_w2,XI,YI,'cubic');
ZI_w3 = interp2(X,Y,Z_w2,XI,YI,'linear');
% find 50 wt% of cumulative interpolated CaCO3 at each interpolated depth
ZI_w3(find(isnan(ZI_w2))) = 0.0;
wtpct_divide = zeros(1,plot_D_n);
for m = 1:plot_D_n,
    loc_Z = ZI_w3(m,:)/sum(ZI_w3(m,:));
    n = 1;
    loc_sum = 0.0;
    while loc_sum <=0.5,
        loc_sum = sum(loc_Z(1:n));
        wtpct_divide(m) = plot_wtpct(n);
        n = n + 1;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT FIGURES ****************************************************** %
% *********************************************************************** %
%
% *** FIGURE 1--Raw Binned Data ********************************************************** %
%
figure;
colormap parula
hold on;
caxis([0.0 1.0]);
contourf(X,Y,Z_n,20,'LineWidth',0.1);
%scatter(data(:,4),data(:,3),'o','Filled','Sizedata',40,'MarkerEdgeColor','k','MarkerFaceColor','w');
% set axes and labels
axis([0 100 D_min D_max]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
xlabel('wt% CaCO_3');
ylabel('paleo-depth (m)');
title('Raw binned data');
% save
%print('-dpsc2', [filename '.' str_window '.raw.' str_date '.ps']);
%
% *** FIGURE 2--Weighted and Normalized Binned Data (no interpolation!) Use Me!! :)********************************************************** %
%
figure;
colormap parula
hold on;
caxis([0.0 1.0]);
contourf(bin_wtpct,bin_D,Z_w2,20,'LineWidth',0.1);
%scatter(data(:,4),data(:,3),'o','Filled','Sizedata',40,'MarkerEdgeColor','k','MarkerFaceColor','w');
plot(wtpct_divide(:),plot_D(:),'w--','linewidth',1.5);
% set axes and labels
axis([0 100 D_min D_max]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
xlabel('wt% CaCO_3');
ylabel('paleo-depth (m)');
title('Weighted and normalized binned data');
% save
%print('-dpsc2', [filename '.' str_window '.weighted_normalized.' str_date '.ps']);
%
% *** FIGURE 3--Weighted and Normalized Binned Data (interpolated only in wtpct space) ********************************************************** %
%
figure;
colormap parula
hold on;
caxis([0.0 1.0]);
contourf(plot_wtpct,bin_D,wtpct_intp,20,'LineWidth',0.1);
%scatter(data(:,4),data(:,3),'o','Filled','Sizedata',40,'MarkerEdgeColor','k','MarkerFaceColor','w');
plot(wtpct_divide(:),plot_D(:),'w--','linewidth',1.5);
% set axes and labels
axis([0 100 D_min D_max]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
xlabel('wt% CaCO_3');
ylabel('paleo-depth (m)');
title('Weighted and normalized binned data -- interpolated (only in wtpct space)');
% save
%print('-dpsc2', [filename '.' str_window '.interp_raw.' str_date '.ps']);
%
% *** FIGURE 4 --Weighted and Normalized Binned Data (interpolated bilinearly)********************************************************** %
%
figure;
colormap parula
hold on;
%%%caxis([0.0 1.0]);
contourf(XI,YI,ZI_w1,20,'LineWidth',0.1);
plot(wtpct_divide(:),plot_D(:),'w--','linewidth',1.5);
%scatter(data(:,4),data(:,3),'o','Filled','Sizedata',40,'MarkerEdgeColor','k','MarkerFaceColor','w');
% set axes and labels
axis([0 100 D_min D_max]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
xlabel('wt% CaCO_3');
ylabel('paleo-depth (m)');
title('Weighted and normalized binned data -- interpolated (linear)');
% save
%print('-dpsc2', [filename '.' str_window '.interp_weighted_raw.' str_date '.ps']);
%
% *** FIGURE 5 --Weighted and Normalized Binned Data (cubic interpolation)********************************************************** %
%
figure;
colormap parula
hold on;
caxis([0.0 1.0]);
contourf(XI,YI,ZI_w2,20,'LineWidth',0.1);
%scatter(data(:,4),data(:,3),'o','Filled','Sizedata',40,'MarkerEdgeColor','k','MarkerFaceColor','w');
plot(wtpct_divide(:),plot_D(:),'w--','linewidth',1.5);
% set axes and labels
axis([0 100 D_min D_max]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
xlabel('wt% CaCO_3');
ylabel('paleo-depth (m)');
title('Weighted and normalized binned data -- interpolated (cubic)');
colorbar
% save
%print('-dpsc2', [filename '.' str_window '.interp_weighted_normalized.' str_date '.ps']);
%
% *********************************************************************** %
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
% *********************************************************************** %
