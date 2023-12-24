warning off                                                                                                                              
                                                                                                                                         
rootpath='D:\Droughts_scales';                                                                                                           
inputpath=[rootpath,'\output'];                                                                                                          
outpath=[rootpath,'\figures'];                                                                                                           
if ~exist(outpath)                                                                                                                       
    mkdir(outpath)  % if not exist, creat a new folder                                                                                   
end                                                                                                                                      
[mask_land,r]=geotiffread('D:\Droughts_scales\landmask_0.5degree.tif');                                                 
load('D:\Droughts_scales\mask_ice_bare.mat');                                                               
id1=find(mask_ice_bare~=1);                                                                                 
id2=find(mask_land==1);                                                                                     
mask=intersect(id1,id2);    
load('D:\Droughts_scales\Areas_05.mat');%the pixel areas at 0.5X0.5 degree
load('D:\Droughts_scales\continental_raction.mat');%the proportion of land areas at 0.5X0.5 degree
Areas=A.*fr;
%Areas=Areas(21:300,:); % mask the areas in polar areas
%%%%%%%%%%%
[mask_land,~]=geotiffread('D:\Droughts_scales\landmask_0.5degree.tif');                                                                                  
load('D:\Droughts_scales\mask_ice_bare.mat');                                                                                                                   
id1=find(mask_ice_bare~=1);                                                                                                                  
id2=find(mask_land==1);                                                                                                                      
mask=intersect(id1,id2);                                                                                                                     
%%Mask data  
load('D:\Droughts_scales\Areas_05.mat')
load('D:\Droughts_scales\continental_raction.mat');
Areas=A.*fr;                      
%%Mask data                                                                                                                                  
[lat,lon] = cdtgrid(0.5);                                                                                                                    
txt={'Arid';'Semi-Arid';'Sub-Humid';'Humid'};                                                                                                
%isin=inpolygon(lon,lat,content.X,content.Y);                                                                                                
[img, ~] = geotiffread('D:\Droughts_scales\ai05.tif');                                                                                                   
%info=geotiffinfo('');                                                                                                                       
img=single(img);                                                                                                                             
img(img==img(1,1))=nan;                                                                                                                      
AI=img;                                                                                                                       
%%%%%%%%%%%%%%%%%%%
[lucc,~]=geotiffread('D:\Droughts_scales\MODIS\LCType05.tif');                                                                              
%imagesc(lucc)                                                                                                                          
lucc=single(lucc);                                                                                                                      
id=find(lucc==0);%0 is ocean                                                                                                                  
lucc(id)=nan; 
lucc(lucc<=10)=1;
lucc(lucc>10)=nan;                                                        
%%%%%%%%%%%%%%%%%%         
%%%%%%%%%%%%%%%%%%%
[ghm,~]=geotiffread('D:\Droughts_scales\input\gHM_84_05_0-1_mean.tif');                                                                                                                                                                                                 
ghm=single(ghm);                                                                                                                      
ghm(ghm==ghm(1,1))=nan;
ghm(ghm>0.6)=nan;                                                                                                                    
ghm(ghm<=0.6)=1;                                                  
%%%%%%%%%%%%%%%%%%    
load('D:\Droughts_scales\input\Vegetation_id.mat');
idv=nan(360,720);
idv(mask)=faster_id;
%%%%%%%%%%%%%%%%%%   
window=17;
bins=34-window+1
%%%%%%%%%%%%%%%%%%%%%%%%%    


ndvi=ndvi/10000; 
[rows,cols,nts]=size(ndvi);
ndvi=reshape(ndvi,rows*cols,12,nts/12);%360*720    12    18
ndvi=ndvi(mask,:,:); 

for n=1:bins
disp(['Round----',num2str(n)])
%calculate monthly Pearson coefficent; 
datestr(now)
Pearson_Rval_ndvi_Spei=[];
Pearson_Pval_ndvi_Spei=[];
%%%%%%%%%%%%%%%%%%%%%%%
for scale=1:24 %%%1-24 month scales
  load(['spei198201-201512_',num2str(scale,'%02d'),'.mat']); %spei01-24
  spei=reshape(spei,rows*cols,12,nts/12);
  spei=spei(mask,:,:);
  
  for m=1:12 %1-12month
       %%%%calculate the m th month£¬r(NDVI,SPEI-j)
       disp([['Round',num2str(n),':--',num2str((scale-1)*12+m),': Month',num2str(m),'---','Scale',num2str(scale)]])
       ndvi_month=squeeze(ndvi(:,m,n:n+window-1));
       ndvi_month_z=zscore(ndvi_month,0,2);
       spei_month=squeeze(spei(:,m,n:n+window-1));
       spei_month_z=zscore(spei_month,0,2);
        tic
        for i=1:length(ndvi)
            x=ndvi_month_z(i,:);
            y=spei_month_z(i,:);
            [r p]=corrcoef(x,y);
            Pearson_Rval_ndvi_Spei(scale).scale(m,i)=r(2,1);
            Pearson_Pval_ndvi_Spei(scale).scale(m,i)=p(2,1);            
        end
        toc
    end %end month
end %end scale

save(['D:\Droughts_scales\output\GIMMS_NDVI_Pearson_Rval_',num2str(n,'%02d'),'_',num2str(window,'%02d'),'.mat'],'Pearson_Rval_ndvi_Spei','-v7.3');
save(['D:\Droughts_scales\output\GIMMS_NDVI_Pearson_Pval_',num2str(n,'%02d'),'_',num2str(window,'%02d'),'.mat'],'Pearson_Pval_ndvi_Spei','-v7.3');

end %end n

disp('finished Pearson')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Calculate maximum Pearson Coefficient and Optimal timescales.
%%%
warning off
cd('D:\Droughts_scales\output')

window=17;
bins=34-window+1

Coef_max_pearson_t=nan(360,720,bins,'single');
Scale_max_pearson_t=nan(360,720,bins,'single');
Month_max_pearson_t=nan(360,720,bins,'single');
Pval_max_pearson_t=nan(360,720,bins,'single');


for n=1:bins
n
Pearson_Rval_ndvi_Spei=importdata(['D:\Droughts_scales\output\GIMMS_NDVI_Pearson_Rval_',num2str(n,'%02d'),'_',num2str(window,'%02d'),'.mat']);
Pearson_Pval_ndvi_Spei=importdata(['D:\Droughts_scales\output\GIMMS_NDVI_Pearson_Pval_',num2str(n,'%02d'),'_',num2str(window,'%02d'),'.mat']);


%%%find optimal timescale with maximum Pearson Coefficient
for j=1:24  
    Pearson_sum((j-1)*12+1:j*12,:)=Pearson_Rval_ndvi_Spei(j).scale;
    Pearson_Pval_sum((j-1)*12+1:j*12,:)=Pearson_Pval_ndvi_Spei(j).scale;
end

[max_a,index]=max(Pearson_sum);
%Coef_max_pearson=max_a;
in=sub2ind(size(Pearson_sum),index,(1:length(index))); % max by col 
%[max_a,index]=max(Pearson_sum,[],2);in=sub2ind(size(Pearson_sum),(1:length(index))',index);% max by rows
Coef_max_pearson=Pearson_sum(in);
Coef_max_pearson(find(isnan(max_a)))=nan;
%%%%
Scale_max_pearson=ceil(index/12);
Scale_max_pearson(find(isnan(max_a)))=nan;
%%%%%Pval
Pval_max_pearson=Pearson_Pval_sum(in);
Pval_max_pearson(find(isnan(max_a)))=nan;
%%%%%
da=nan(360,720,'single');
da(mask)=Coef_max_pearson;
Coef_max_pearson=da;
%
da=nan(360,720,'single');
da(mask)=Scale_max_pearson;
Scale_max_pearson=da;
%
da=nan(360,720,'single');
da(mask)=Pval_max_pearson;
Pval_max_pearson=da;
%%%%%
Coef_max_pearson_t(:,:,n)=Coef_max_pearson;
Scale_max_pearson_t(:,:,n)=Scale_max_pearson;
Pval_max_pearson_t(:,:,n)=Pval_max_pearson;
end %end n

%%%%%
save(['Coef_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat'],'Coef_max_pearson_t','-v7.3');
save(['Scale_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat'],'Scale_max_pearson_t','-v7.3');
save(['Pval_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat'],'Pval_max_pearson_t','-v7.3');
disp('finished Time-scales & Maxmium R')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('load data')
warning off
cd('D:\Droughts_scales\output')
load(['Coef_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Scale_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Pval_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);

r=Coef_max_pearson_t;
s=Scale_max_pearson_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%slope%%%%%%%%%%%%%%
%calculate%%TSd%slope%%%%%%%%%%%%%
disp('Start slope  Time-scales')
data=reshape(Scale_max_pearson_t,360*720,bins);
x=(1:size(data,2))';
%%%%%%%
slopes=nan(length(data),1);
sigs=nan(length(data),1);
tic
for i=1:length(data)
if all(~isnan(data(i,:)'))
[~, ~, ~, sig, ~, ~, ~,slope, ~, ~, ~, ~, ~, ~ ,~ ,~] = ktaub([x data(i,:)'],0.05);   
slopes(i)=slope;
sigs(i) = sig;
end
end
toc
%%
%%%%%%%%%%%%%%%%%%%%%%%
scale_slope =reshape(slopes,360,720);
scale_slope(lucc~=1)=nan;
scale_slope(idv~=1)=nan;         
scale_slope(ghm~=1)=nan;      
scale_sig =reshape(sigs,360,720);
scale_sig(lucc~=1)=nan;
scale_sig(idv~=1)=nan;         
scale_sig(ghm~=1)=nan;     
%%%%%%%
cd('D:\Droughts_scales\output')
save(['TSd_max_Pearson_GIMMS_NDVI_Slope_',num2str(window,'%02d'),'.mat'],'scale_slope'); 
save(['TSd_max_Pearson_GIMMS_NDVI_Sig_',num2str(window,'%02d'),'.mat'],'scale_sig'); 
disp('finished slope Time-scales')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate%%R%slope%%%%%%%%%%%%%
disp('Start slope R')
data=reshape(Coef_max_pearson_t,360*720,bins);
x=(1:size(data,2))';
%%%%%%%
slopes=nan(length(data),1);
sigs=nan(length(data),1);
tic
for i=1:length(data)
if all(~isnan(data(i,:)'))
[~, ~, ~, sig, ~, ~, ~,slope, ~, ~, ~, ~, ~, ~ ,~ ,~] = ktaub([x data(i,:)'],0.05);  
slopes(i)=slope;
sigs(i) = sig;
end
end
toc
%%
%%%%%%%%%%%%%%%%%%%%%%%
r_slope =reshape(slopes,360,720);
r_slope(lucc~=1)=nan;
r_slope(idv~=1)=nan;         
r_slope(ghm~=1)=nan;      
r_sig =reshape(sigs,360,720);
r_sig(lucc~=1)=nan;
r_sig(idv~=1)=nan;         
r_sig(ghm~=1)=nan;    
%%%%%%%
cd('D:\Droughts_scales\output')
save(['Rd_max_Pearson_GIMMS_NDVI_Slope_',num2str(window,'%02d'),'.mat'],'r_slope'); 
save(['Rd_max_Pearson_GIMMS_NDVI_Sig_',num2str(window,'%02d'),'.mat'],'r_sig'); 
disp('finished slope R')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start corr(Time-Scale,R)')
warning off
cd('D:\Droughts_scales\output')
load(['Coef_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Scale_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Pval_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);

r=Coef_max_pearson_t;
s=Scale_max_pearson_t;


s=reshape(s,360*720,bins);
r=reshape(r,360*720,bins);
rva=r(mask,:);
sca=s(mask,:);
R=nan(size(rva,1),1);
P=nan(size(rva,1),1);
tic
for n=1:length(rva)
            x=rva(n,:);
            y=sca(n,:);
            [r,p]=corr(x',y','type','spearman'); 
            R(n)=r;
            P(n)=p;  
end
toc

da=nan(360,720);
da(mask)=R;
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
R=da;

da=nan(360,720);
da(mask)=P;
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
P=da;

save(['Rval_ScalevsR_max_Pearson_GIMMS_NDVI_R_',num2str(window,'%02d'),'.mat'],'R'); 
save(['Pval_ScalevsR_max_Pearson_GIMMS_NDVI_Sig_',num2str(window,'%02d'),'.mat'],'P'); 
disp('finished corr(Time-Scale,R)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start spatial corr(TSd,Rd)')
warning off
cd('D:\Droughts_scales\output')
load(['Coef_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Scale_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);
load(['Pval_max_Pearson_Rval_GIMMS_NDVI_Spei_t_',num2str(window,'%02d'),'.mat']);

r=Coef_max_pearson_t;
s=Scale_max_pearson_t;

R=nan(size(s,3),1);
P=nan(size(s,3),1);
for t=1:size(s,3)
t
da=r(:,:,t);
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
x=da;
da=s(:,:,t);
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
y=da;
z=x.*y;

id=find(~isnan(z));
x=x(id);
y=y(id);
[r0,p0]=corr(x,y,'type','spearman');
R(t)=r0;
P(t)=p0;  
end
cd('D:\Droughts_scales\output')
save(['Rt_TSd-Rd_max_Pearson_GIMMS_NDVI_R_',num2str(window,'%02d'),'.mat'],'R'); 
save(['Pt_TSd-Rd_max_Pearson_GIMMS_NDVI_Sig_',num2str(window,'%02d'),'.mat'],'P'); 
disp('finished corr(spatial Time-Scale,R)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%calculate global area-weighted mean of TSd and Rd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
y=[];
ci95=[];
sen=Scale_max_pearson_t;
num=size(sen,3);
sen(isinf(sen))=nan;

quantile(sen(:),0.9999),quantile(sen(:),0.0001) 
 %%%%%%%%%%      
for t=1:size(sen,3)
da=sen(:,:,t);
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
%da=da(21:300,:,:);

id=find(~isnan(da));da=da(id);a=Areas(id);a=a./nanmax(a(:));b=da.*a; %area-weighted
y(t)=nanmean(b(:)) ;
%calculated 95% values of the t-distributio
data=b/nansum(a(:));
N = size(data,1);                 
yMean = nanmean(data);                 
ySEM = nanstd(data)/sqrt(N);            
CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:)); 
ci95(t)=yCI95(2);
end
cd('D:\Droughts_scales\output')
save(['Global_Area-weighted-mean_TS_max_Pearson_GIMMS_NDVI_y_',num2str(window,'%02d'),'.mat'],'y'); 
save(['Global_Area-weighted-mean_TS_max_Pearson_GIMMS_NDVI_ci95_',num2str(window,'%02d'),'.mat'],'ci95'); 
%%%%%%%%%%

%%%%%%%%%%
sen=Coef_max_pearson_t;
y=[];
ci95=[];                                                                                                 
for t=1:size(sen,3)
da=sen(:,:,t);
da(lucc~=1)=nan;
da(idv~=1)=nan;         
da(ghm~=1)=nan;    
%da=da(21:300,:,:);
id=find(~isnan(da));da=da(id);a=Areas(id);a=a./nanmax(a(:));b=da.*a; %area-weighted
y(t)=nanmean(b(:)) ;
%calculated 95% values of the t-distributio
data=b;
N = size(data,1);                    
yMean = nanmean(data);                 
ySEM = nanstd(data)/sqrt(N);         
CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:)); 
ci95(t)=yCI95(2);
end
save(['Global_Area-weighted-mean_R_max_Pearson_GIMMS_NDVI_y_',num2str(window,'%02d'),'.mat'],'y'); 
save(['Global_Area-weighted-mean_R_max_Pearson_GIMMS_NDVI_ci95_',num2str(window,'%02d'),'.mat'],'ci95'); 
%%%%%%%%%%
close;[l,p] = boundedline((1:length(y)), y,ci95, '-b','alpha');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
