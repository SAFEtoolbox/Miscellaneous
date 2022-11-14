
MM = [ 2, 3, 8, 20, nan ] ;
fun = {'Liu','Homma','Sobol','Morris'};
Nfun = length(fun)  ;

% % CODE TO LOAD DATA FROM AB.pawn.xlsx FILE GENERATE BUY PUY ET AL CODE
% TAKES A WHILE TO RUN. SKIP AND DIRECTLY LOAD THE .mat FILES IF AVAILABLE 
%
% k = 16384 ;
% ini = nan(size(MM)) ;
% ini(1) = 2 ; % row where first data is
% for i=2:Nfun+1
%     ini(i) = ini(i-1)+k*MM(i-1) ;
% end
% for i = 1:Nfun ;
% 
% M = MM(i) ;
% N = xlsread('AB.pawn.xlsx',[ 'C' num2str(ini(i)) ':C' num2str(ini(i+1)-1) ]);
% n = xlsread('AB.pawn.xlsx',[ 'D' num2str(ini(i)) ':D' num2str(ini(i+1)-1) ]);
% S = xlsread('AB.pawn.xlsx',[ 'H' num2str(ini(i)) ':H' num2str(ini(i+1)-1) ]);
% S = reshape(S,k,M) ;
% N = reshape(N,k,M) ;
% n = reshape(n,k,M) ;
% 
% eval(['N_' fun{i} ' =N(:,1);'])
% eval(['n_' fun{i} ' =n(:,1);'])
% eval(['S_' fun{i} ' =S;'])
%
% end

load all % workspace saved from previous Matlab session 
% where the above code was run to generate matlab structures

% Reproduce second row of Figure 2 in the paper (boxplots of PAWN indices
% when varying tuning parameters n,N)
figure(1)
sh = 0.01 ;
dx = (1-0.1*2-sh*(Nfun-1))/(sum(MM(1:end-1)));
Dx = dx.*MM(1:end-1)'       ;
dy = 1-0.1*2 ;
Dy = dy*ones(Nfun,1)        ;
pos_x = nan(Nfun,1); pos_x(1) = 0.1 ; for i=2:Nfun; pos_x(i)=pos_x(i-1)+Dx(i-1)+sh; end
pos_y = 1-0.1-Dy ;
h_pos = [pos_x , pos_y, Dx, Dy ] ;
for i = 1:Nfun ; 
    eval(['N=N_' fun{i} ';'])
    eval(['n=n_' fun{i} ';'])
    eval(['S=S_' fun{i} ';'])
    h = subplot(1,Nfun,i);
    boxplot( S ); set(gca,'YLim',[-0.3,1.3],'XLim',[0,MM(i)+1])
    h.Position = h_pos(i,:); 
    set(gca,'YTick',[0,0.5,1]); title(fun{i})
    if i == 1 ; ylabel('PAWN'); else set(gca,'YTickLabel',{}); end
    xlabel('inputs')    
end

% Generate Figure B in our Review
%
% Define threshols above which index values of first input factor X1
% are taken as 'outliers' (one threshold for each benchmark function) 
S_ul = [ 0.41, 0.53, 0.57, 0.21 ] ;
figure(2) % Figure B
tmp = unique(n) ;
for i = 1:Nfun ;  
    eval(['N=N_' fun{i} ';'])
    eval(['n=n_' fun{i} ';'])
    eval(['S=S_' fun{i} ';'])
    % Determine values of tuning parameters (n,N) that generated outliers:
    idx = S(:,1)>S_ul(i); 
    % Plot:
    subplot(1,Nfun,i); title(fun{i})
    for j=1:length(tmp); plot([tmp(j),tmp(j)],[min(N),max(N)],'k'); hold on; end
    plot(n(idx),N(idx),'xr');
    axis([min(n)-1,max(n)+1,0,max(N)+min(N)])
    if i == 1 ; ylabel('N'); else set(gca,'YTickLabel',{}); end
    xlabel('n')  
end


% Generate Figure C in our Review
%
Nc_min = [ 50, 80 ] ;
Nt = length(Nc_min) ;
sh = 0.01 ;
sv = 0.05 ; 
dx = (1-0.1*2-sh*(Nfun-1))/(sum(MM(1:end-1)));
Dx = dx.*MM(1:end-1)'       ;
dy = (1-0.1*2-sv*(Nt-1))/Nt ;
Dy = dy*ones(Nfun,1)        ;
pos_x = nan(Nfun,1); pos_x(1) = 0.1 ; for i=2:Nfun; pos_x(i)=pos_x(i-1)+Dx(i-1)+sh; end
pos_y = nan(1,Nt)  ; pos_y(1) = 1-0.1-Dy(1); for i=2:Nt;   pos_y(i)=pos_y(i-1)-Dy(i)-sv  ; end; 
tmp = repmat(pos_y,Nfun,1) ;
h_pos = [repmat(pos_x,Nt,1) , tmp(:), repmat(Dx,Nt,1), repmat(Dy,Nt,1) ] ;
k = 1 ;
figure(3)
for j = 1:Nt    
    for i = 1:Nfun ;
        eval(['N=N_' fun{i} ';'])
        eval(['n=n_' fun{i} ';'])
        eval(['S=S_' fun{i} ';'])       
        figure(3)
        h = subplot(length(Nc_min),Nfun,k);
        boxplot( S( N./n>Nc_min(j) ,:) ); set(gca,'YLim',[-0.3,1.3],'XLim',[0,MM(i)+1])
        h.Position = h_pos(k,:);
        set(gca,'YTick',[0,0.5,1]);
        k=k+1;        
        if i == 1 ; ylabel('PAWN'); else set(gca,'YTickLabel',{}); end
        if j == Nt; xlabel('inputs'); else set(gca,'XTickLabel',{}); end
        if i == 3 ; title(['Use only ' num2str(sum(N./n>Nc_min(j))) ' experiments where N_c = N/n > ' num2str(Nc_min(j)) ]) ; end        
    end
end







