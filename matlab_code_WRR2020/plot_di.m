
function [] = plot_di(Si,deri_M_deg,M,T,Rain,Qobs,DistrPar,Ylegend,W)


load colorbar1

SiM = nan(T,M);
for i = 1:M
    for t = 1:T
        SiM(t,i) = Si{t,1}(i);
    end
end
    

for i = 1:M
    
    fig1 = figure; %FP
    set(fig1,'defaultAxesColorOrder',[0 0 0; 0 0 0]); %FP
       
    %  plot Si
    SiM_p = SiM(:,i);
    SiM_p(:,2) = 1;
    SiM_p(:,3) = 0;
           
    h1 = subplot(2,1,1) ;
    set(h1,'position',[0.1300 0.805 0.80 0.028]);
    c = imagesc(transpose(SiM_p));
    colormap(h1,flip(gray));
    
    axis([0.5,T+0.5,0+0.5,1+0.5]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    ylabel('SI');
    set(gca,'Fontsize',14);
    
    hold on
    
    line1x = [W+0.5,W+0.5];
    line1y = [0,100];
    plot(line1x,line1y,'--','Color',[0.4 0.4 0.4],'LineWidth',2);
    hold on
    
    line2x = [T-W+0.5,T-W+0.5];
    line2y = [0,100];
    plot(line2x,line2y,'--','Color',[0.4 0.4 0.4],'LineWidth',2);
    % plot DI
    
    deri_M_deg_p = deri_M_deg{i};
    deri_M_deg_p = [deri_M_deg_p,90*ones(101,1),-90*ones(101,1)];
    
    h2 = subplot(2,1,2);
    set(h2,'position',[0.1300 0.35 0.80 0.4]);
    
    c = imagesc(flip(deri_M_deg_p,1));
    
    colormap(h2,clrs)
    
    axis([0.5,T+0.5,0+0.5,100+0.5]);
    
    ylabelnum = 0:25:100;
    
    yticknum = 100-ylabelnum;
    Ylabel = {};
    
    for j = 1:length(ylabelnum)
        
        Ylabel{j} = num2str(ylabelnum(j)/100*(DistrPar{1}(2)-DistrPar{1}(1))+DistrPar{1}(1));
        
    end
    
    set(gca,'YTick',flip(yticknum+0.5),'YTickLabel',flip(Ylabel));
    set(gca,'Fontsize',14);
    
    
    xlabel('Time step (day)');
    ylabel(Ylegend{i});
    
    hold on
    
    
    line1x = [W+0.5,W+0.5];
    line1y = [0,100];
    plot(line1x,line1y,'--','Color',[0.4 0.4 0.4],'LineWidth',2);
    
    hold on
    
    line2x = [T-W+0.5,T-W+0.5];
    line2y = [0,100];
    plot(line2x,line2y,'--','Color',[0.4 0.4 0.4],'LineWidth',2);
    
    % plot rainfall and discharge
    
    hold on
    
    Rainp = Rain/1.5/1.2/2;
    
    
    
    bar(Rainp,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0.5,0.5,0.5]);
    set(gca,'ydir','reverse');
    %     set(gca,'XTick',[]);
    set(gca,'Fontsize',12);
    %     set(gca,'ylim',[0 100]);
    %     set(gca,'YTick',0:50:250);
    %     ylabel('Rianfall(mm/day)');
    
    yyaxis right
    
    Qobsp = Qobs/50/1.2/2;
    
    plot(Qobsp,'LineWidth',1.2);
    
    set(gca,'ylim',[0 1]);
    set(gca,'YTick',[]);
    %     set(gca,'YTickLabel',[]);
    %     ylabel('Discharge(mm/day)');
    set(gca,'Fontsize',12);
    
end