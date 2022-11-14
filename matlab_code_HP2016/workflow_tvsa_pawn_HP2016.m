% This workflow scripts implements the time-varying GSA
% of parameter versus data uncertainty reported in:
%
% Pianosi, F. and Wagener, T. (2016), Understanding the time-varying 
% importance of different uncertainty sources in hydrological modelling 
% using global sensitivity analysis, Hydrological Processes, 
% doi:10.1002/hyp.10968
% http://onlinelibrary.wiley.com/doi/10.1002/hyp.10968/abstract
%
% The code works together with the SAFE Toolbox. More information and
% dowload at: www.bris.ac.uk/cabot/resources/safe-toolbox/
%
% This code developed by Francesca Pianosi, University of Bristol
% email: francesca.pianosi@bristol.ac.uk

clear all

safe_dir = '/Users/francesca/SAFE_code/safe_R1.1';

% Addpath of SAFE:
addpath(genpath(safe_dir))

% Load and prepare data:
%filename = '03451500.txt'; River='French'; snow = 0 ;% French Broad River at Ashville
filename = '08167500.txt'; River='Guadalupe'; snow = 0 ;% Guadalupe Rv near Spring Branch
%filename = '05455500.txt'; River='English'; snow = 1 ;% English River at Kalona
ini_date =[ 1948 10 1  ];
fin_date =[ 1951 09 30 ];
[prec,ept,temp,flow,tt,dy,ly] = load_mopex_data(filename,ini_date,fin_date);

%% STEP 1: Define parameter ranges

% a) Sample parameters:
Npar = 10000;
Xpar_lab = {'BETA','LP','FC','PERC','K0','K1','K2','UZL','MAXBAS'};
xmin =       [ 0   0.3    1   0     0.05 0.01 0.05   0 1] ; 
xmax =       [ 7   1   2000 100     1    1    0.1  100 6] ;
M = length(xmin); 
DistrFun = cell(1,M); for i=1:M-1; DistrFun{i} = 'unif' ; end; DistrFun{end} = 'unid' ;
DistrPar = cell(M,1); for i=1:M-1; DistrPar{i} = [ xmin(i) xmax(i) ] ; end; DistrPar{end} = xmax(end) ;
Xpar = AAT_sampling('lhs',M,DistrFun,DistrPar,Npar);
% snow parameters (will be actually used in snow-affected catchments only):
Xpar_lab_snow = {'TS','CFMAX','CFR','CWH'};
xmin_snow =[ -3  0 0 0.00 ] ; 
xmax_snow =[  3 20 1 0.8  ] ;
M = length(xmin_snow); 
DistrPar = cell(M,1); for i=1:M; DistrPar{i} = [ xmin_snow(i) xmax_snow(i) ] ; end
Xpar_snow = AAT_sampling('lhs',M,'unif',DistrPar,Npar);

% b) Run simulations:
Case   =  1; % Case=1: interflow is dominant / Case = 2: percolation is dominant
warmup =  0; % days for warmup
T      = length(prec) ;
Q_sim  = nan(T,Npar)  ;
NobjFun= 3 ; % number of objective functions
ObjFun = nan(Npar,NobjFun) ;
for i=1:Npar    
    if snow == 1
        P = snow_routine(temp,prec,Xpar_snow(i,:)) ; % snow-affected
    else
        P = prec ; 
    end
    %
    Qs = hbv_sim(P,ept,Xpar(i,:),Case);
    Q_sim(:,i)  = Qs ;
    ObjFun(i,1) = sqrt(mean((Qs(warmup+1:end) - flow(warmup+1:end)).^2)); %RMSE
    ObjFun(i,2) = mean(abs(Qs(warmup+1:end) - flow(warmup+1:end)))      ; %AME
    ObjFun(i,3) = abs(mean(Qs(warmup+1:end) - flow(warmup+1:end)))      ; %BIAS
end

% c) Set and evaluate thresholds for objective functions:
%thresholds = [ Inf Inf Inf ] ; % Use this if you dont want any threshold
%thresholds = [ 1 0.8 0.1 ] ; % reasonable thresholds for French
thresholds = [ 0.12 0.1 0.1 ] ; % reasonable thresholds for Guadalupe
%thresholds = [ 1.1 0.43 0.5 ] ; % reasonable thresholds for English
% Find parameter sets that satisfy all the thresholds at same time:
idx = sum( ObjFun <= repmat(thresholds,Npar,1 ) , 2 ) == NobjFun ;
% Plot results (scatter plot):
figure
k=1; M = size(Xpar,2) ;
for j=1:NobjFun
    for i=1:M
        subaxis(NobjFun,M,k,'SpacingVert', 0.05,'SpacingHoriz',0.01)
        plot(Xpar(idx,i),ObjFun(idx,j),'.')
        if i~=1; set(gca,'YTickLabel',{}); end 
        if j~=NobjFun; set(gca,'XTickLabel',{}); else xlabel(Xpar_lab{i}); end
        k=k+1;
    end
end
sum(idx)
% Plot results (PCP):
figure; parcoor(Xpar,Xpar_lab,[],idx)
% Plot results (hydrograph):
figure
plot(Q_sim(:,idx),'k'); hold on
plot(flow,'r')

% PCP of snow parameters:
figure; parcoor(Xpar_snow,Xpar_lab_snow,[],idx)

% d) Based on previous analysis, define refined parameter ranges:
[ min(Xpar(idx,:)) ; max(Xpar(idx,:)) ]
% Rounding for French:
%xmin_ref = [ 0.4 0.3 125  10 0.15 0.10 0.05 28 1 ] ;
%xmax_ref = [ 2.6 1.0 500 100 1.50 0.50 0.10 90 4 ] ;
% Rounding for Guadalupe:
xmin_ref = [ 1.5 0.3  100   3 0.05 0.05 0.05   0 1 ] ;
xmax_ref = [ 6.0 1.0 1000 100 1.00 1.00 0.10 100 4 ] ;
% Rounding for English:
%xmin_ref = [ 1.0 0.3  50   0 0.05 0.10 0.05   0 1 ] ;
%xmax_ref = [ 7.0 1.0 600 100 2.00 1.00 0.10 100 4 ] ;
%[ min(Xpar_snow(idx,:)) ; max(Xpar_snow(idx,:)) ]
%xmin_snow_ref = [ -3  1 0 0    ] ;
%xmax_snow_ref = [  3 20 1 0.8  ] ;
  
%% STEP 2: CREATE LIST OF INPUT FACTORS SAMPLES

% a) create list of soil parameters:
Npar = 5000 ;
DistrPar = cell(3,1); for i=1:3; DistrPar{i} = [ xmin_ref(i) xmax_ref(i) ] ; end
Xparsoil = AAT_sampling('lhs',3,'unif',DistrPar,Npar);

% b) create list of routing parameters:
DistrFun = cell(1,6); for i=1:5; DistrFun{i} = 'unif' ; end; DistrFun{end} = 'unid' ;
DistrPar = cell(6,1); for i=1:5; DistrPar{i} = [ xmin_ref(3+i) xmax_ref(3+i) ] ; end; DistrPar{end} = xmax_ref(end) ;
Xparroute = AAT_sampling('lhs',6,DistrFun,DistrPar,Npar);

% c) create list of snow parameters ('dummy' list if not snow-affected):
if snow == 1
DistrPar = cell(4,1); for i=1:4; DistrPar{i} = [ xmin_snow_ref(i) xmax_snow_ref(i) ] ; end
    Xparsnow = AAT_sampling('lhs',4,'unif',DistrPar,Npar);
else
    Xparsnow = repmat([-Inf,0,0,0],Npar,1); % not snow-affected
end

% d) Create list of perturbed rainfall data:
rM = 0.4  ; % allowed ratio of variability (e.g. rM=0.4 means +/-40%)
Nr = 5000 ; % number of perturbed time series to be generated
Xprec = create_prec_list(prec,rM,Nr) ;
%figure; plot(Xprec)

% e) Create list of potential evaporation data:
eM = 0.2  ; % allowed ratio of variability
Ne = 5000 ; % number of perturbed time series to be generated
Xept = create_ept_list(ept,eM,Ne) ;
%figure; plot(Xept)

% f) Create list of perturbed flow data:
eQ = 0.4  ; % allowed ratio of variability
Nf = 5000 ; % number of perturbed time series to be generated
alfa = 0.8; 
Xflow = create_flow_list(flow,eQ,alfa,Nf);
%figure; plot(Xflow)

%% STEP 3: Apply PAWN

% a) define some additional tuning parameters:
W = 15; % semi-length of window to compute time-varying objective fun
X_labels = {'Prec','PET','Snow','Soil','Route','Flow'} ;
M = length(X_labels) ;

% b) Unconditional sampling and model evaluation:
NU = 2000 ;
XU = AAT_sampling('rsu',M,'unid',{Nr,Ne,Npar,Npar,Npar,Nf},NU);% matrix (NU,M)
[TSU,TSobsU] = hbv_snow_sim_1(XU,Xprec,temp,Xept,Xflow,Xparsnow,Xparsoil,Xparroute,Case);
RMSEtU       = tv_objfun(TSU(:,warmup+1:end),TSobsU(:,warmup+1:end),'RMSE',W) ;

% c) Conditional sampling and model evaluation:
n  = 20   ;
NC = 1000 ;
[ XC, xc ] = pawn_sampling('rsu',M,'unid',{Nr,Ne,Npar,Npar,Npar,Nf},n,NC);
RMSEtC = cell(M,n);
TSC    = cell(M,n);
tic
for i=1:M
    for k=1:n
        [TSsim,TSobs] = hbv_snow_sim_1(XC{i,k},Xprec,temp,Xept,Xflow,Xparsnow,Xparsoil,Xparroute,Case);
        TSC{i,k}      = TSsim ;
        RMSEtC{i,k}   = tv_objfun(TSsim(:,warmup+1:end),TSobs(:,warmup+1:end),'RMSE',W) ;
    end
    fprintf('Input %d/%d completed\n',i,M);
end; comp_time=toc;

% d) Create support for estimation of CDFs:
meanflow = nan(T,1) ;
for t=1:T
    meanflow(t)=mean(flow(max(1,t-W):min(T,t+W))) ;
end
threshold = 2*meanflow ;
[YF, Y_min,Y_max] = cdf_support(RMSEtU,RMSEtC,threshold,100) ;

% e) Compute KS:
[KS_max_RMSE, KS] = pawn_ks_timevarying(RMSEtU,RMSEtC,YF,threshold)   ;
% Plot:
fs = 26;
hfig=figure;
clrs = gray ;
clrs = clrs(end:-1:1,:);
colormap(clrs)
imagesc4pdf([KS_max_RMSE ones(M,1) zeros(M,1)])
%colorbar('location','WestOutside')
colorbar
hold on
plot(M+0.5-flow/max(flow)*M,'r','LineWidth',2)
set(gca,'YTick',1:M,'YTickLabel',X_labels,'XTick',dy,'XTickLabel',ly,'XLim',[365 T],'FontSize',fs)
set(hfig, 'Position', [0 0 900 400])
box on
ylabel(River,'FontSize',fs)
%file_name = [ num2str(filename(1:end-4)) '_tv'] ;
%saveas(hfig,file_name,'fig');

% f) Compare with non-TV GSA:
%RMSEcU = RMSE(TSU,flow');
RMSEcU = RMSE(TSU(:,366:end),Xflow(366:end,XU(:,6))');
RMSEcC = cell(M,n);
for i=1:M
    for k=1:n
        %RMSEcC{i,k} = RMSE(TSC{i,k},flow');
        TSsim = TSC{i,k};
        XCik  = XC{i,k} ;
        TSobs = Xflow(:,XCik(:,6))' ;
        RMSEcC{i,k} = RMSE(TSsim(:,366:end),TSobs(:,366:end));
    end
end
Nboot = 100 ;
KS_max_RMSEc = pawn_indices( RMSEcU, RMSEcC, 'max',[], Nboot,[],'below',2*mean(flow) ) ;
% Plot:
hfig=figure;
KS_max_RMSEc_ = KS_max_RMSEc(end:-1:1);
X_labels_ = X_labels(end:-1:1) ;
barh(KS_max_RMSEc_,'FaceColor',[0.5 .5 .5]);set(gca,'YTickLabel',X_labels_)
set(gca,'FontSize',fs,'XLim',[0,max(KS_max_RMSEc)+0.01],'YLim',[0.5,M+0.5])
set(hfig, 'Position', [0 0 300 400])
%file_name = [ num2str(filename(1:end-4)) '_con'] ;
%saveas(hfig,file_name,'fig');

%% Analysis changing window size:
W2 = 1 ;
warmup = 0 ;
RMSEtU2 = tv_objfun(TSU(:,warmup+1:end),TSobsU(:,warmup+1:end),'RMSE',W2) ;
RMSEtC2 = cell(M,n);
for i=1:M
    for k=1:n
        TSsim = TSC{i,k};
        XCik  = XC{i,k} ;
        TSobs = Xflow(:,XCik(:,6))' ;
        RMSEtC2{i,k} = tv_objfun(TSsim(:,warmup+1:end),TSobs(:,warmup+1:end),'RMSE',W2) ;
    end
    fprintf('Input %d/%d completed\n',i,M);
end
% Create support for estimation of CDFs:
YF2 = cdf_support(RMSEtU2,RMSEtC2,threshold,100) ;
% Compute KS:
[KS_max_RMSE2,KS2] = pawn_ks_timevarying(RMSEtU2,RMSEtC2,YF2,threshold) ;
% Plot:
fs = 26;
hfig=figure;
clrs = gray ;
clrs = clrs(end:-1:1,:);
colormap(clrs)
imagesc4pdf([KS_max_RMSE2 ones(M,1) zeros(M,1)])
%colorbar('location','WestOutside')
colorbar
hold on
plot(M+0.5-flow/max(flow)*M,'r','LineWidth',2)
set(gca,'YTick',1:M,'YTickLabel',X_labels,'XTick',dy,'XTickLabel',ly,'XLim',[365 T],'FontSize',fs)
set(hfig, 'Position', [0 0 900 400])
box on
ylabel(River,'FontSize',fs)
%file_name = [ num2str(filename(1:end-4)) '_tv_w2'] ;
%saveas(hfig,file_name,'fig');

%% Apply DYNIA

perc  = 10 ; % Percent of the input/output sample size
nbins = 5 ; 
X_tmp = Xparsoil ; ii=4 ; X_labels_tmp = {'BETA','LP','FC'};
%X_tmp = Xparroute ; ii=5 ; X_labels_tmp = {'PERC','K0','K1','K2','UZL','MAXBAS'};
for i=1:length(X_labels_tmp)
    % Compute posterior probability:
    [ xi, fi ] = dynia_histograms(X_tmp(XU(:,ii),i),RMSEtU,perc,nbins) ;    
    % Plot results:
    dynia_plot(xi,fi,X_labels_tmp{i},flow) ;
    set(gca,'XLim',[366,T])
    hold on
    plot(0.5+prec/max(prec)*nbins,'b','LineWidth',2)
end

%% Save results
%file_name = [ num2str(filename(1:end-4)) '_res'];
%save(file_name);

%% More detailed Analysis

% Detailed analysis for Guadalupe
thresholds = [ 0.12 0.1 0.1 ] ;
idx = sum( ObjFun <= repmat(thresholds,size(Xpar,1),1 ) , 2 ) == NobjFun ;
fs = 24;
hfig=figure; hold on 
plot(Q_sim(:,idx),'Color',[126 126 126]/256)
%plot(Xflow(:,xc{end}),'b')
plot(flow,'r','LineWidth',1.5)
plot(1.1-prec/max(prec)*max(flow)/10,'k','LineWidth',2.5)
xlabel('time (day)','FontSize',fs)
ylabel('Flow (mm/day)','FontSize',fs)
axis([365*2,T,0,1.1])
box on
set(gca,'FontSize',fs)
set(hfig, 'Position', [0 0 630 400])
%
X_tmp = Xparsoil ; ii = 4 ; % group: soil parameters
X_labels_tmp = {'BETA','LP','FC'};
i  = 1 ; % parameter in group 'soil'
perc  = 5 ; nbins = 5 ; 
[ xi, fi ] = dynia_histograms(X_tmp(XU(:,ii),i),RMSEtU,perc,nbins) ;    
hfig = dynia_plot(xi,fi,X_labels_tmp{i},flow) ;
set(hfig, 'Position', [0 0 700 400])
set(gca,'XLim',[365*2,T])
xlabel('time (day)','FontSize',fs)
plot(nbins+0.5-flow/max(flow(365*2:T))*(nbins-1),'k','LineWidth',2)
%plot(prec/max(prec(365*2:T))*(nbins-1)-0.5,'k')

figure
plot(reshape(flow(1:365*3),365,3))
% effect of 'BETA':
fs = 24 ;
FC = max(X_tmp(:,3));
sm = [0:FC] ;
hfig=figure; hold on
BETA1 = 1.96 ;
BETA2 = 5.55 ;
plot((sm/FC).^BETA1,'k','LineWidth',2)
plot((sm/FC).^BETA2,':k','LineWidth',2)
legend(['BETA = ' num2str(BETA1) ],['BETA = ' num2str(BETA2) ])
xlabel('soil moisture (mm)','FontSize',fs); ylabel('Runoff/Prec','FontSize',fs)
set(gca,'FontSize',fs)
set(hfig, 'Position', [0 0 500 400])
box on
% ... increasing BETA reduces runoff...

% Detailed analysis for English:
t = 500 ; % time step
i = 6   ; % uncertainty source
% Plot flow perturbation:
pert = Xflow(t-W:t+W,xc{i}) - repmat(flow(t-W:t+W),1,n) ;
[tmp,idx1] = max(mean(pert)) ;
[tmp,idx2] = min(mean(pert)) ;
hfig=figure;
hold on
plot(t-W:t+W,zeros(size(pert,1)),':k')
plot(t-W:t+W,pert(:,idx1),':k','LineWidth',2)
plot(t-W:t+W,pert(:,idx2),'k','LineWidth',2)
set(gca,'FontSize',fs)
box on
xlabel('time (day)','FontSize',fs)
ylabel('flow perturbation (mm/day)','FontSize',fs)
axis([t-W,t+W,-0.8,0.8])
set(hfig, 'Position', [0 0 450 400])
%file_name = [ '0' num2str(filename(1:end-4)) '_pert'] ;
%saveas(hfig,file_name,'fig');

% Associated conditional CDFs:
FU  = empiricalcdf(RMSEtU(:,t),YF{t});
RMSEtCik1 = RMSEtC{i,idx1} ;
FC1 = empiricalcdf(RMSEtCik1(:,t),YF{t});
RMSEtCik2 = RMSEtC{i,idx2} ;
FC2 = empiricalcdf(RMSEtCik2(:,t),YF{t});
hfig=figure;
hold on
plot(YF{t},FU,'r','LineWidth',2)
plot(YF{t},FC1,':k','LineWidth',2)
plot(YF{t},FC2,'k','LineWidth',2)
set(gca,'FontSize',fs,'XLim',[0.5,1.1])
box on
xlabel('RMSE (mm/day)','FontSize',fs)
ylabel('CDF','FontSize',fs)
set(hfig, 'Position', [0 0 450 400])
%file_name = [ '0' num2str(filename(1:end-4)) '_cdf'] ;
%saveas(hfig,file_name,'fig');

