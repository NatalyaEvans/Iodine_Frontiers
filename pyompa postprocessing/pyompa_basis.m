%% 

% for WOD data with nitrite reox

%% Load in data
load('ETNP_df_pyompasoln')

% res=sqrt((base.CT_resid.^2)/12+(base.SA_resid.^2)/8+(base.Phosphate_resid.^2)/6+(base.Nitrate_resid.^2)/4+(base.Silicate_resid.^2)/2+(base.tCO2_resid.^2)/4);
% base(res>4,:)=[]; % removes large residuals

load('WOD_data_C_filter'); % loads in water mass evolution data
col=[[0 0.4470 0.7410]',[0.5500 0.4250 0.0980]',[0.9290 0.6940 0.1250]',[0.4940 0.1840 0.5560]',[0.4660 0.6740 0.1880]',[0.3010 0.7450 0.9330]',[0.6350 0.0780 0.1840]',[1 0 0]']; %colors for iterating on loops

wm(1,:)=[2.35500000000000;2.65000000000000;3.13000000000000]; % po4
wm(2,:)=[31.2000000000000;35.6000000000000;43.21000000000000]; % no3
wm_names={'13CW','NEPIW','AAIW'};
aero=[[1,14]',[1,17]'];
anaero=[[1,-65]',[1,-40]'];

xshift=0.075;
yshift=0.075*(1/3);

clear box
box(1,:)=[-xshift+yshift,-xshift-yshift,xshift-yshift,xshift+yshift,-xshift+yshift];
xshift=xshift*15.5;
yshift=yshift*(52.5);
box(2,:)=[-xshift-yshift,-xshift+yshift,xshift+yshift,xshift-yshift,-xshift-yshift];

%% N-P

ind=base.pdens>=26.4 & base.pdens<=27.2;
no3=base.Nitrate(ind);
po4=base.Phosphate(ind);

figure(1)
p4=plot(po4,no3,'ko','MarkerSize',3,'MarkerFaceColor','k'); % 
hold on
i=3;

% plot WOA data
p1=plot(AAIW_Cfilt.Phosphatemmolkg, AAIW_Cfilt.Nitratemmolkg,'Color',col(:,i),'Marker','x','LineStyle','none','MarkerSize',5);
hold on
i=i-1;

p2=plot(EqPIW_Cfilt.Phosphatemmolkg, EqPIW_Cfilt.Nitratemmolkg,'Color',col(:,i),'Marker','x','LineStyle','none','MarkerSize',5);
i=i-1;

p3=plot(x13CW_Cfilt.Phosphatemmolkg, x13CW_Cfilt.Nitratemmolkg,'Color',col(:,i),'Marker','x','LineStyle','none','MarkerSize',5);
i=i-1;

p5=plot(wm(1,:),wm(2,:),'kd','MarkerFaceColor','r','MarkerSize',8);
text(wm(1,:)-[0.1,0.2,0.2],wm(2,:)+[6.5,4.5,2.5],wm_names)

shape=[wm(:,1),wm(:,1)+aero(:,1),wm(:,1)+aero(:,2)];
p6=fill(shape(1,:),shape(2,:),'magenta','FaceAlpha',0.3,'EdgeColor','None');

shape=[wm(:,1),wm(:,1)+1.1*aero(:,1),wm(:,1)+1.1*aero(:,2)];
arrowX=[shape(1,2)-0.06, shape(1,2), shape(1,2)-0.02];
arrowY=[shape(2,2), shape(2,2), shape(2,2)-1.25];
plot(shape(1,1:2),shape(2,1:2),'k')
plot(arrowX,arrowY,'k')

arrowX=[shape(1,3)-0.06, shape(1,3), shape(1,3)-0.01];
arrowY=[shape(2,3), shape(2,3), shape(2,3)-1.5];
plot(shape(1,[1,3]),shape(2,[1,3]),'k')
plot(arrowX,arrowY,'k')

shape=[wm(:,2),wm(:,2)+0.5*anaero(:,1),wm(:,2)+0.5*anaero(:,2)];
p7=fill(shape(1,:),shape(2,:),'green','FaceAlpha',0.3,'EdgeColor','None');

arrowX=[shape(1,2)-0.1, shape(1,2), shape(1,2)+0.06];
arrowY=[shape(2,2)+1, shape(2,2), shape(2,2)+2.2];
plot(shape(1,1:2),shape(2,1:2),'k')
plot(arrowX,arrowY,'k')

arrowX=[shape(1,3)-0.07, shape(1,3), shape(1,3)+0.02];
arrowY=[shape(2,3)+0.75, shape(2,3), shape(2,3)+1.9];
plot(shape(1,[1,3]),shape(2,[1,3]),'k')
plot(arrowX,arrowY,'k')

for i=1:3
%     display(i)
%     wm(1,i)+box(1,:)
%     wm(2,i)+box(2,:)
    p8=plot(wm(1,i)+box(1,:),wm(2,i)+box(2,:),'r','LineWidth',1.5);
end

xlim([1 3.5])
xlabel('PO_4^{3-}/{\mu}mol kg^{-1}')
ylabel('NO_3^{-}/{\mu}mol kg^{-1}')
[h,icons] = legend([p3, p2, p1, p4, p5, p6, p7],'WOD18 13CW','WOD18 NEPIW','WOD18 AAIW','FK180624 data','Endmembers','Aerobic remin.','Anaerobic remin.','Location','South','NumColumns',2);
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons(8:10),'MarkerSize',10);

hold off

