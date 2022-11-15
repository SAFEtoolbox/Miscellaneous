% This script provide an example of the Direction of Change (DOC)
% assessment within the context of Time-Varying GSA, as described in:
%
% Wang, A., Pianosi, F., Wagener, T. (2020) A Diagnostic Approach to
% Analyze the Direction of Change in Model Outputs Based on Global
% Variations in the Model Inputs, Water Resources Research, 56(8)
% https://doi.org/10.1029/2020WR027153 
%
% This script was prepared by Anqi Wang, Hohai University, Nanjing, China, 
% in May 2020. Email to: wanganqi0718@163.com
%
% IMPORTANT NOTE!
% This script requires the HDMR decomposition software presented in:
%
% Ziehn, T., & Tomlin, A. S. (2009). GUI-HDMR A software tool 
% for global sensitivity analysis of complex models. 
% Environmental Modelling & Software, 24(7), 775-785. 
% https://doi.org/10.1016/j.envsoft.2008.12.002)
%
% The software can be requested from the authors through the website:
% http://www.gui-hdmr.de
% 
% MODEL AND STUDY AREA
% This workflows uses a dataset of Monte Carlo simulations
% of the rainfall-runoff Hymod model applied the Leaf catchment 
% in Mississipi, USA (see header of file LeafCatch.txt for more details).
% More information about Hymod can be found in the 'example' folder 
% of the SAFE Toolbox.
% Note that, with respect to that documentation of Hymod, two (of the five)
% model parameters are here defined in a slightly different 
% (but conceptually equivalent) way. Specifically, the parameter
% definitions, and ranges, used here are:
% 
% The inputs subject to SA are the 5 model parameters, and the scalar 
% output for SA is the total discharge with a moving window.
% SM = [0,400 ] mm
% beta = [0,2] (-)
% alfa = [0,1] (-)
% Rs = [10,150] days (this is the multiplactive inverse of the adimensional
% parameter % used in SAFE) 
% Rf = [1,10] days (same as above)

%% Step 1: set paths

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';
% 
% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath([ my_dir '/sampling'])
addpath([ my_dir '/util'])
addpath([ my_dir '/visualization'])
addpath([ my_dir '/example/hymod'])
addpath([ my_dir '/VBSA'])

%% Step 2:load input files

% Load MC simulations results:
load LHS_1000_s1.mat

% set number of input:
[N,M] = size(parsets);

% load rainfall and discharge:
load -ascii LeafCatch.txt;
Rain = LeafCatch(1:365*3,1);
Evap = LeafCatch(1:365*3,2);
Qobs = LeafCatch(1:365*3,3);

% remove warmup period:
Nwarm = 90;
Qobs(1:Nwarm) = [];
Rain(1:Nwarm) = [];
outMCS(:,1:Nwarm) = [];

% unit transfer m3/s to mm/day
Qsim = outMCS/1949/1000*24*60*60;

% sample size and time steps
[N,T] = size(Qsim);

%% Step 3:apply moving window as objective function

W = 15;  % W is the half window size, full window size is 2*W+1
Q1 = nan(N,T);
for i = 1:N   
    for t = 1:T
        Q1(i,t) = sum(Qsim(i,(max(1,t-W):min(T,t+W))));
    end
end

% Plot:
figure; i = 88; % choose (randomly) which output sample to plot
plot(Qsim(i,:),'k'); hold on; plot(Qobs,'.b'); plot(Q1(i,:)/(2*W+1),'r')
xlabel('Time (days)'); ylabel('Flow (mm/day)'); legend('Sim','Obs','Sim - moving average')

%% Step 4: rescale input and output between 0 and 1

% Define input ranges:
DistrPar  = { [ 0 400 ]; [ 0 2 ]; [ 0 1 ]; [ 10 150 ] ; [ 1 10 ] } ; 

% Rescale:
for i = 1:M
    input(:,i) = unifcdf(parsets(:,i),DistrPar{i}(1),DistrPar{i}(2));
end

% Rescale output:
Q = Q1/max(Q1(:));

%% Step 5: Approximation of the HDMR components fi(xi)

% Load the (pre-computed) coefficients of the orthonormal basis functions
% (functions phi in Eq. (8) of Anqi et al (2020)):
global ORTHOPOLYCOEFLEG
CoefOrthPolyLeg15;
% Here we load 15 coefficients, which means the polynomial order 'k' 
% in Eq. (8) will have to be 15 at most 
% (this is a reasonable assumption as we usually use much lower order)

% Set maximum polynomial order:
max_1st = 10; % max value of 'k' (see Eq. 8)
max_2nd = 3 ; % max value of the HDMR decomposition for calculating
% 2nd order Sobol' sensitivity index (not shown in paper) - this is needed
% because the HDMR code we are going to use for estimating first-order
% Sobol' sensitivity indices (see Step 7 below) also simoultaneously 
% estimates second-order sensitivity indices, so even if we are not going
% to use them we still need to estimate them

% Set number of iterations for variance reduction:
it_1st = 10;

% Estimate first-order coefficients of HDMR decomposition:

% connect to folder with GUI-HDMR software:
addpath('gui_hdmr_software_r2009') 
% (https://doi.org/10.1016/j.envsoft.2008.12.002)

for t = 1:T % loop over time steps
    
    % calculate f0 (Eq. 7), alpha (Eq. 8) and beta (not shown in the paper)   
    [alpha{t,1},f0(t,1)]=sub_alpha_1st(input,Q(:,t),max_1st);
    [beta_numbers] = sub_beta_numbers;
    [beta{t,1}]=sub_beta_2nd(input,Q(:,t),beta_numbers,max_2nd);
    
    % iteration process for variance-reduction:
    for s=1:it_1st 
        % calculate optimal order for 1st-order polynomials:    
        [opt_1st{t,1}] = sub_opt_order_1st_e1(f0(t,1),alpha{t,1},input,Q(:,t),max_1st);       
        % calculate analytical reference function h(x) 
        % with current optimal order and alpha:
        approx1st{t,1} = sub_comp_opt_1st(alpha{t,1},f0(t,1),opt_1st{t,1},input,Q(:,t));
        % calculate new alpha:
        alpha_new{t,1} = sub_get_new_alpha_ratio(input,Q(:,t),approx1st{t,1},alpha{t,1},max_1st);
        alpha{t,1}     = alpha_new{t,1};
        
    end
    
    % calculate optimal polynomial order:
    [opt_2nd{t,1},ii{t},jj{t},n_ij{t}] = sub_opt_order_2nd_e1_alt(f0(t,1),alpha{t,1},beta{t,1},input,Q(:,t),...
        beta_numbers,opt_1st{t,1},max_2nd);
    
end

%% Step 6: Compute derivatives of fi(xi) and DI indices

% Compute fi(xi) (Eq. 8):
for t = 1:T % loop over time steps
    
    for i = 1:M % loop over input factors
        x = input(:,i);
        fixi_sub = zeros(N,opt_1st{t,1}(i));
        for j = 1:opt_1st{t,1}(i)
            fixi_sub(:,j) = alpha{t,1}(i,j)*ortho_nom(j,x);
        end
        fixi{t,1}(:,i) = sum(fixi_sub,2);
    end
    
end

% Compute derivatives of fi(xi):
for t = 1:T
    
    for i = 1:M
    
        x = linspace(0,1,101);
        fixi_deri_sub = zeros(101,opt_1st{t,1}(i));
        for j = 1:opt_1st{t,1}(i)
            fixi_deri_sub(:,j) = alpha{t,1}(i,j)*ortho_nom_deri(j,x);
        end
        fixi_deri{t,1}(:,i) = sum(fixi_deri_sub,2);
    end


end

% Rearrange results (this will make it easier to plot results later on):
for i = 1:M
    deri_M{i,1} = zeros(101,T);
    for t = 1:T
       deri_M{i,1}(:,t) = fixi_deri{t}(:,i) ;
    end
end

% Transform derivative values into angle values:
for i = 1:M
    deri_M_rad{i,1} = atan(deri_M{i}); % DI (eq. 6) expressed in radiant
    deri_M_deg{i,1} = rad2deg(deri_M_rad{i,1}); % DI (eq. 6) expressed in degree
end

DI = deri_M_deg; % these are the Direction Indices

%% Step 7: Approximate first order Sobol' Indices Using the HDMR:

for t=1:T % loop over time-steps
    [D(t,1),Di{t,1},Dij,Si{t,1},Sij] = sub_sensitivity_indices(alpha{t,1},beta{t,1},opt_1st{t,1},...
        opt_2nd{t,1},input,Q(:,t));
end
    
%% Step 8: visualise results

% define parameter labels:
Ylegend = {'SM (mm)','Beta (-)','Alfa (-)','Rs (day)','Rf (day)'};

% Plot 1: Scatter plot of fi(xi) at one time step (chosen by user)

% choose time step:
pt = 400;
% plot:
for i = 1:M
    
    fig = figure;   
    set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);    
    plot(input(:,i),Q(:,pt),'x','Color',[0.5,0.5,0.5]);
    hold on
    plot(input(:,i),fixi{pt,1}(:,i)+f0(pt,1),'k.');
    
    % Customise:
    xtick1 = 0:0.2:1;
    for j = 1:length(xtick1)
        xlabel1{j} = num2str(xtick1(j)*(DistrPar{i}(2)-DistrPar{i}(1))+DistrPar{i}(1));
    end
    set(gca,'XTick',xtick1);
    set(gca,'XTicklabel',xlabel1);
    set(gca,'ylim',[0 1]);
    ytick1 = 0:0.2:1;
    for j = 1:length(ytick1)
        ylabel1{j} = num2str(round(ytick1(j)*max(Q1(:)),2));
    end
    set(gca,'YTick',ytick1);
    set(gca,'YTicklabel',ylabel1);
    posi = [0.3 0.3 0.5 0.5];
    set(gca,'Position',posi);   
    set(gca,'Fontsize',14);
    xlabel(Ylegend{i});
    ylabel('Discharge volume(mm/day)');    
    ax1 = gca;
    ax2 =  axes('Position',posi,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k');
    set(ax2,'ylim',[0 1]);
    set(ax2,'YTick',0:0.2:1);
    set(ax2,'xlim',[0 1]);
    set(ax2,'XTick',0:0.2:1);
    set(gca,'YTicklabel',ylabel1);
    xlabel('Rescaled input');
    ylabel('Rescaled output');
    set(ax2,'Fontsize',14);
    
end

% Plot 2 (similar to Fig. 2):
% Time series of DI indices:
plot_di(Si,DI,M,T,Rain,Qobs,DistrPar,Ylegend,W);




%% Supplementary: fix anomalous time step (see Sec. 2.3)

% set fix time step

t = 909;

% load sample input and output
sample_input=input;
sample_output=Q(:,t);

% set maximum polynomial order
max_1st = 10;
max_2nd = 3;

% calculate f0, alpha and beta
[alpha1,f01]=sub_alpha_1st(sample_input,sample_output,max_1st);
[beta_numbers] = sub_beta_numbers;
[beta1]=sub_beta_2nd(sample_input,sample_output,beta_numbers,max_2nd);

% calculate optimal polynomial order
[opt_1st1] = sub_opt_order_1st_e1(f01,alpha1,sample_input,sample_output,max_1st);
% display optimal order
[opt_1st_nonzero1] = sub_get_fac_fi_nonzero(opt_1st1,size(sample_input,2));
[approx1st1] = sub_comp_opt_1st(alpha1,f01,opt_1st1,sample_input,sample_output);
[opt_2nd1,i1,j1,n_ij1] = sub_opt_order_2nd_e1_alt(f01,alpha1,beta1,sample_input,sample_output,...
                     beta_numbers,opt_1st1,max_2nd);
% display optimal order
[opt_2nd_nonzero1] = sub_get_fac_fij_nonzero(opt_2nd1,i1,j1,n_ij1);
 

% calculate first and second order sensitivity indices
[D1,Di1,Dij1,Si1,Sij1] = sub_sensitivity_indices(alpha1,beta1,opt_1st1,...
                    opt_2nd1,sample_input,sample_output);

% replace matrix

Si{t} = Si1;
opt_1st{t} = opt_1st1;
alpha{t} = alpha1;
