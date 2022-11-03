function [fig,para4] = striped_interp(x,y,z,detail,intertype,extraptype)
% Creates a filled contour plot using high res vertical interpolation then
% low res horizontal interpolation
%   This function takes an x, y, and z variable of the same size and makes
%   a contourf plot with points where sampled. x, y, and z all need to be
%   the same size
%   detail corresponds to the vertical resolution. A single number
%   intertype corresponds to the type of interpolation in interp1
%   extraptype tells the function if it extrapolates ('extrap') or not
%   ('NaN')

chem_all=z; % change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=y; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);

if exist('detail','var')==0
    detail=1;
end

if exist('intertype','var')==0
    intertype='linear';
end

if exist('extraptype','var')==0
    extraptype=NaN;
end

filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
xp=x(filt); %filters the x-axis to keep the same dimensions
para=chem_all(filt); %filters the heat map
XIc=linspace(min(xp),max(xp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(para));
x1=xp(~isnan(para));
para=para(~isnan(para));

depth_range=[min(y_var1):detail:max(y_var1)];

%determine the number of stations and the positions for them
stn_coords=unique(round(x1,2,'significant'));

for counter_var=1:length(stn_coords) % for each station
    ind=find(round(x1,2,'significant')==stn_coords(counter_var)); % identify the index of each station in the broader data set
    y_var2=y_var1(ind);
    para2=para(ind);
    
    [~, uniques] = unique(y_var2); % ind = index of first occurrence of a repeated value
    y_var3=y_var2(uniques);
    para3=para2(uniques);
    
    [~, uniques] = unique(y_var3); % ind = index of first occurrence of a repeated value
    if length(uniques)>1
        para4(:,counter_var)=interp1(y_var3(uniques),para3(uniques),depth_range,intertype,extraptype); % interpolate data onto pre-definied grid, linearly, exptrapolated values are set to NaN
    end
    if length(uniques)<=1
        para4=zeros(length(depth_range),length(stn_coords));
    end
end

fig=contourf(stn_coords,depth_range,para4,[min(para4,[],'all'):(max(para4,[],'all')-min(para4,[],'all'))/100:max(para4,[],'all')],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(x1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
set(gca,'box','on')
hold off

end

