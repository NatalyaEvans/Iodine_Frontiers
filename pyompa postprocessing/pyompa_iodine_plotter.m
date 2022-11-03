%% 

%% Load in data
load('ETNP_df_pyompasoln')

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

res=sqrt(df.Consvtemp_resid.^2+df.Abssal_resid.^2+df.PO4_resid.^2+df.NO3_resid.^2);

%% linear regression instead of SVD

df2=df;
df2(isnan(df2.aero_proj)==1,:)=[];
df2(isnan(df2.anaero_proj)==1,:)=[];
% df2.no2_reox_proj=df2.anaero_proj.*df2.no2_reox;
% lme2 = fitlme(df2,'Iodate ~ -1 + CW_frac_total + NEPIW_frac_total + AAIW_frac_total + aero_proj + no2_reox_proj')
lme2 = fitlme(df2,'Iodate ~ -1 + CW_frac_total + NEPIW_frac_total + AAIW_frac_total + aero_proj + anaero_proj')

lme2.Rsquared

a=lme2.Coefficients.Estimate;
sa=lme2.Coefficients.SE;

writetable(dataset2table(lme2.Coefficients),'io3_deconvolution.csv')

sim_I=[df2.CW_frac_total, df2.NEPIW_frac_total, df2.AAIW_frac_total, df2.aero_proj, df2.anaero_proj]*a; % calculated the simulate iodate
diff_I=sim_I'-df2.Iodate; % calculate the residual iodate

%% plot regression results as a scatter plot

load('transMap');
transMap=transMap.transMap;

y=df2.Iodate;
figure(1)
plot(y,sim_I,'ko','MarkerFaceColor','k')
hold on
[m,b,r,sm,sb] = lsqfitma(y,sim_I);
plot(sort(y),m*sort(y)+b,'k--')
plot([min(y):1:max(sim_I)],[min(y):1:max(sim_I)],'m-.')
{'m','sm','b','sb'}
[m,sm,b,sb]
% ylim([0 600])
% title('{\sigma}_{\theta} = 26.4-27.2 kg m^{-3}')
xlabel('Iodate/nM')
ylabel('Simulated iodate/nM')
legend('Data','Best fit line','1:1 line','Location','Southeast')
box on
hold off

%% plot SVD results as a section plot

chem_all=df.Iodate; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(df.pdens); %in pressure or potential density
max_depth=max(df.pdens);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}';
%cmax=800; %Sets the max color for interpolation. For iodine use 504 from Rin's paper. NO2 looks good with 2.2
figure(2); %declare the figure. Running this script again will overwrite the figure unless you change the number
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
%caxis([0 cmax]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('haline');
colormap(cmap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string','Iodate/nM','FontSize',12);
% xlabel(['Longitude/' char(176) 'E'])
set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

% saveas(gcf,'io3_pdens_section.svg');

%% plot regression results as a section plot

chem_all=diff_I; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(df.pdens); %in pressure or potential density
max_depth=max(df.pdens);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='natural'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}'; 
%cmax=800; %Sets the max color for interpolation. For iodine use 504 from Rin's paper. NO2 looks good with 2.2
figure(3); %declare the figure. Running this script again will overwrite the figure unless you change the number
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
%caxis([0 cmax]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
colormap(transMap)
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string','Simulated - measured iodate/nM','FontSize',12);
% xlabel(['Longitude/' char(176) 'E'])
set(gca,'XTickLabels',[])

ylabel(depth_var)
ylim([26.4 27.2])
set(gca,'box','on')
hold off

% saveas(gcf,'sim_meas_io3_section.svg');


%% Plot iodate consumed by anaerobic remineralization

chem_all=a(5)*df.anaero_proj; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(df.pdens); %in pressure or potential density
max_depth=max(df.pdens);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='natural'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}'; 
%cmax=800; %Sets the max color for interpolation. For iodine use 504 from Rin's paper. NO2 looks good with 2.2
figure(4); %declare the figure. Running this script again will overwrite the figure unless you change the number
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
h=contourf(XIc,YIc,chem2); %does a better interpolated color fit than contourf
contourf(XIc,YIc,chem2,[floor(min(chem1)):10:floor(max(chem1))],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
caxis([-400 -100]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('-dense');
colormap(cmap);
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Iodate removed by', 'anaerobic remineralization/nM'},'FontSize',12);
xlabel(['Longitude/' char(176) 'E'])
ylabel(depth_var)
set(gca,'box','on')
ylim([26.4 27.2])
hold off

% saveas(gcf,'sim_io3_anaero.svg');


%% Plot iodate removed by the contribution of 13CW to the water

chem_all=(a(1)-a(2))*df.CW_frac_total; %change this assignment to alter what you are plotting. RHS is the name in the workspace defined in import_data
y_var_all=df.pdens; %selects press or pdens to plot, default is press
min_depth=min(df.pdens); %in pressure or potential density
max_depth=max(df.pdens);
detail=100; %how thorough the chem interpolation is. The smaller the number, the smaller the area of interpolation will be and the less defined. Higher values risk creating hot spots and an artifically lower blanket appearance
intertype='linear'; %type of interpolation. Default is linear, but cubic, natural, neighbor, and v4 are available. More information can be found at https://www.mathworks.com/help/matlab/ref/griddata.html#bvkwume-method
depth_var='Potential density/kg m^{-3}'; 
%cmax=800; %Sets the max color for interpolation. For iodine use 504 from Rin's paper. NO2 looks good with 2.2
figure(5); %declare the figure. Running this script again will overwrite the figure unless you change the number
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
h=contourf(XIc,YIc,chem2); %does a better interpolated color fit than contourf
contourf(XIc,YIc,chem2,[floor(min(chem1)):1:floor(max(chem1))],'LineStyle','None');
hold on %allows you to plot the sample points over the transect
plot(long1,y_var1,'k.','markersize',6) %plots the sample points
set(gca,'YDir','reverse'); %flips the depth axis
caxis([-350 0]) %sets the range of colorbar. Comment out for an auto scale
h=colorbar; %displays the colorbar
cmap=cmocean('-dense');
colormap(cmap);
set(gca,'position',[.1 .1 .69 .6])
set(get(h,'label'),'string',{'Iodate removed by 13CW content/nM'},'FontSize',12);
xlabel(['Longitude/' char(176) 'E'])
ylabel(depth_var)
set(gca,'box','on')
ylim([26.4 27.2])
hold off

% saveas(gcf,'sim_io3_13CW.svg');

%% Estimate remin rate in the ODZ

df.Depthm=-gsw_z_from_p(df.press,df.lat);

% anaerobic remin rates for C using Devol and Hartnett measurements
depths=[100 145 190 310 322 500 620 800 1020];
cox=[5.52 5.36 5.57 3.62 3.08 2.15 3.55 2.81 1.44]; % mmol C/m2 day
cox_err=[0.71 0.68 0.76 0.45 0.39 0.35 0.44 0.38 0.19];

pin=[7.4,0.36];
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(depths',cox',pin,'modfunc');
out = p(1).*(depths/100).^p(2);

figure(8)
errorbar(cox,depths,[],[],cox_err,cox_err,'ko-','MarkerFaceColor','k')
hold on
plot(out,depths,'m--')
axis ij
xlabel('Carbon oxidation rate/mmol m^{-2} day^{-1}')
ylabel('Depth/m')
legend('Data','Best fit curve','Location','Southeast')
ylim([0 1100])
box on
grid on
hold off

p(1); % this is the remin rate of C at 100 m in mmol/m^3 day, or uM/day

% rate of iodate reduction

% a in nM IO3/uM PO4
% p in mmol/m^3 day
iod_red=[-lme2.Coefficients.Estimate(5), lme2.Coefficients.SE(5)]*(p(1)/106)
iod_red=[-lme2.Coefficients.Estimate(5), lme2.Coefficients.SE(5)]*(out(4)/106) % using the 310 m point because it is anoxic and a good curve fit





