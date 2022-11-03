%% 

close all

%% Load in data
load('ETNP_df_pyompasoln') % df is the mean_endmember_skeleton, std_df is the std_endmember_skeleton

% for mean_endmembers
refmag=16^2+1^2; % reference stoichiometry
ref=[16,1];
aero_angle=[df.NO3_to_aerobic_remin_ratio,df.PO4_to_aerobic_remin_ratio]; % define the vector direction
df.aero_proj=df.aerobic_remin.*aero_angle*ref'./refmag; % decouple nitrite reoxidation from our remineralization

refmag=94.4^2+1^2; % reference stoichiometry
ref=[-94.4,1];
anaero_angle=[df.NO3_to_anaerobic_remin_ratio,df.PO4_to_anaerobic_remin_ratio]; % define the vector direction
df.anaero_proj=df.anaerobic_remin.*anaero_angle*ref'./refmag; % decouple nitrite reoxidation from our remineralization

ref=-94.4;
df.no2_reox=100*(ref-df.NO3_to_anaerobic_remin_ratio)./ref;

% for std_endmembers
refmag=16^2+1^2; % reference stoichiometry
ref=[16,1];
aero_angle=[std_df.NO3_to_aerobic_remin_ratio,std_df.PO4_to_aerobic_remin_ratio]; % define the vector direction
std_df.aero_proj=std_df.aerobic_remin.*aero_angle*ref'./refmag; % decouple nitrite reoxidation from our remineralization

refmag=94.4^2+1^2; % reference stoichiometry
ref=[-94.4,1];
anaero_angle=[std_df.NO3_to_anaerobic_remin_ratio,std_df.PO4_to_anaerobic_remin_ratio]; % define the vector direction
std_df.anaero_proj=std_df.anaerobic_remin.*anaero_angle*ref'./refmag; % decouple nitrite reoxidation from our remineralization

ref=-94.4;
std_df.no2_reox=100*(ref-std_df.NO3_to_anaerobic_remin_ratio)./ref;

res=sqrt(df.Consvtemp_resid.^2+df.Abssal_resid.^2+df.PO4_resid.^2+df.NO3_resid.^2);


%% plot water masses and water mass uncs

wms={'13CW','NEPIW','AAIW'};

for i=1:3
    figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
    chem_all=100*table2array(df(:,28+i)); %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
    y_var_all=df.pdens; %selects press or pdens to plot, default is press
    min_depth=min(y_var_all); %in pressure or potential density
    max_depth=max(y_var_all);
    detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
    intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
    depth_var='Potential density/kg m^{-3}';
    filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
    y_var=y_var_all(filt); %filters the y-axis
    longp=df.long(filt); %filters the x-axis to keep the same dimensions
    chem=chem_all(filt); %filters the heat map
    XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
    YIc=linspace(min(y_var),max(y_var),detail);
    y_var1=y_var(~isnan(chem));
    long1=longp(~isnan(chem));
    chem1=chem(~isnan(chem));
    chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
    h=contourf(XIc,YIc,chem2,[floor(min(chem1)):1:floor(max(chem1))],'LineStyle','None');
    hold on %allows you to plot the sample points over the transect
    plot(long1,y_var1,'k.','markersize',6) %plots the sample points
    set(gca,'YDir','reverse'); %flips the depth axis
    caxis([0 100]) %sets the range of colorbar. Comment out for an auto scale
    h=colorbar; %displays the colorbar
    cmap=cmocean('tempo');
    colormap(cmap)
    h.Ruler.TickLabelFormat='%g%%';
    set(gca,'position',[.1 .1 .69 .6])
    set(get(h,'label'),'string','Water mass content','FontSize',12);
    title(wms{i});
    % xlabel(['Longitude/' char(176) 'E'])
    set(gca,'XTickLabels',[])
    
    ylabel(depth_var)
    ylim([26.4 27.2])
    set(gca,'box','on')
    hold off
    
    saveas(gcf,[wms{i} '.svg']);
    
%     std_endmember
    figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
    chem_all=100*table2array(std_df(:,28+i))/3; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
    y_var_all=std_df.pdens; %selects press or pdens to plot, default is press
    min_depth=min(y_var_all); %in pressure or potential density
    max_depth=max(y_var_all);
    detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
    intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
    depth_var='Potential density/kg m^{-3}';
    filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
    y_var=y_var_all(filt); %filters the y-axis
    longp=df.long(filt); %filters the x-axis to keep the same dimensions
    chem=chem_all(filt); %filters the heat map
    XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
    YIc=linspace(min(y_var),max(y_var),detail);
    y_var1=y_var(~isnan(chem));
    long1=longp(~isnan(chem));
    chem1=chem(~isnan(chem));
    chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
    h=contourf(XIc,YIc,chem2,[floor(min(chem1)):0.01:floor(max(chem1))],'LineStyle','None');
    hold on %allows you to plot the sample points over the transect
    plot(long1,y_var1,'k.','markersize',6) %plots the sample points
    set(gca,'YDir','reverse'); %flips the depth axis
    caxis([0 0.25]) %sets the range of colorbar. Comment out for an auto scale
    h=colorbar; %displays the colorbar
    colormap(cmap)
    h.Ruler.TickLabelFormat='%g%%';
    set(gca,'position',[.1 .1 .69 .6])
%     set(get(h,'label'),'string',['Uncertainty in ' wms{i} ' content'],'FontSize',12);
    set(get(h,'label'),'string','Uncertainty in water mass content','FontSize',12);

%     title(wms{i});
    xlabel(['Longitude/' char(176) 'E'])
%     set(gca,'XTickLabels',[])
    
    ylabel(depth_var)
    ylim([26.4 27.2])
    set(gca,'box','on')
    hold off
    
    saveas(gcf,[wms{i} ' unc.svg']);
end


%% aerobic remin

figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
chem_all=df.aero_proj; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
longp=df.long(filt); %filters the x-axis to keep the same dimensions
chem=chem_all(filt); %filters the heat map
XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(chem));
long1=longp(~isnan(chem));
chem1=chem(~isnan(chem));
chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
h=contourf(XIc,YIc,chem2,[min(chem1):0.05:max(chem1)],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
% caxis([-0.5 2]) %sets the range of colorbar. Comment out for an auto scale
caxis([0 2]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('amp');
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Aerobic remineralization/{\mu}mol kg^{-1}','PO_4^{3-} equivalents'},'FontSize',12);
% xlabel(['Longitude/' char(176) 'E'])
set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

saveas(gcf,'aero.svg');

%     std_endmember
figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
chem_all=std_df.aero_proj/3; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=std_df.pdens; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
longp=df.long(filt); %filters the x-axis to keep the same dimensions
chem=chem_all(filt); %filters the heat map
XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(chem));
long1=longp(~isnan(chem));
chem1=chem(~isnan(chem));
chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
h=contourf(XIc,YIc,chem2,[min(chem1):0.00001:max(chem1)],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
caxis([0 0.001]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Unc in aerobic remineralization/','{\mu}mol kg^{-1} PO_4^{3-} equivalents'},'FontSize',12);
% title(wms{i});
xlabel(['Longitude/' char(176) 'E'])
% set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

saveas(gcf,'aero_unc.svg');

%% anaerobic remin

figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
chem_all=df.anaero_proj; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
longp=df.long(filt); %filters the x-axis to keep the same dimensions
chem=chem_all(filt); %filters the heat map
XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(chem));
long1=longp(~isnan(chem));
chem1=chem(~isnan(chem));
chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
h=contourf(XIc,YIc,chem2,[min(chem1):0.01:max(chem1)],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
% caxis([-0.5 2]) %sets the range of colorbar. Comment out for an auto scale
caxis([0.05 0.4]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('algae');
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Anaerobic remineralization/{\mu}mol kg^{-1}','PO_4^{3-} equivalents'},'FontSize',12);
% xlabel(['Longitude/' char(176) 'E'])
set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

saveas(gcf,'anaero.svg');

%     std_endmember
figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
chem_all=abs(std_df.anaero_proj)/3; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=std_df.pdens; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
longp=df.long(filt); %filters the x-axis to keep the same dimensions
chem=chem_all(filt); %filters the heat map
XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(chem));
long1=longp(~isnan(chem));
chem1=chem(~isnan(chem));
chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
h=contourf(XIc,YIc,chem2,[min(chem1):0.00001:max(chem1)],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
caxis([0 0.003]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Unc in anaerobic remineralization/','{\mu}mol kg^{-1} PO_4^{3-} equivalents'},'FontSize',12);
% title(wms{i});
xlabel(['Longitude/' char(176) 'E'])
% set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

saveas(gcf,'anaero_unc.svg');

%% residual

figure(); %declare the figure. Running this script again will overwrite the figure unless you change the number
chem_all=res; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(y_var_all); %in pressure or potential density
max_depth=max(y_var_all);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
filt=(y_var_all < max_depth) & (y_var_all > min_depth); %logical indexing for filt as a logical vector
y_var=y_var_all(filt); %filters the y-axis
longp=df.long(filt); %filters the x-axis to keep the same dimensions
chem=chem_all(filt); %filters the heat map
XIc=linspace(min(longp),max(longp),detail)'; % each is a matrix with regular spacing in X and Y dimensions
YIc=linspace(min(y_var),max(y_var),detail);
y_var1=y_var(~isnan(chem));
long1=longp(~isnan(chem));
chem1=chem(~isnan(chem));
chem2=griddata(long1,y_var1,chem1,XIc,YIc,intertype);
h=contourf(XIc,YIc,chem2,[min(chem1):0.001:max(chem1)],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
caxis([0 0.06]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('tempo');
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string','Residuals of fit','FontSize',12);
xlabel(['Longitude/' char(176) 'E'])
% set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

saveas(gcf,'resid.tif');


