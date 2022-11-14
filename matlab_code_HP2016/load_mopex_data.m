function [prec,ept,temp,flow,time,dy,ly] = load_mopex_data(filename,ini_date,fin_date)
%
% [prec,ept,temp,flow,time,dy,ly] = load_mopex_data(filename,ini_date,fin_date)
%
% EXAMPLE:
% filename = '03451500.txt'
% ini_date =[ 1948 10 1  ];
% fin_date =[ 1951 09 30 ];

% Load data from txt file:
data   = load(filename,'-ascii');

% Extract data corresponding to the above 'date':
t_start= find(ismember(data(:,1:3),ini_date,'rows')>0);
if isempty(t_start); error('ini_date too early, no data available before [%d %d %d]',data(1,1),data(1,2),data(1,3)); end
t_end  = find(ismember(data(:,1:3),fin_date,'rows')>0);
if isempty(t_end); error('fin_date too late, no data available after [%d %d %d]',data(end,1),data(end,2),data(end,3)); end
time   = t_start:t_end          ;
prec   = data(time,4)           ; % mm/day
ept    = data(time,5)           ; % mm/day
flow   = data(time,6)           ; % mm/day
temp   = mean(data(time,7:8),2) ; % C

% % Plot data:
%Nsteps = length(time) ;
dy = repmat([ 31 30 31 31 28 31 30 31 30 31 31 30 ],1,3); dy = [1 cumsum(dy)];
ly = {'O','N','D','J','F','M','A','M','J','J','A','S'};
% figure
% subaxis(4,1,1,'SpacingHoriz',0.01,'SpacingVert',0.01)
% plot(prec,'k'); ylabel('prec (mm/day)')
% set(gca,'XTick',dy,'XTickLabel',{},'XLim',[1,Nsteps])
% subaxis(4,1,2,'SpacingHoriz',0.01,'SpacingVert',0.01)
% plot(temp,'k'); ylabel('temp (C)')
% set(gca,'XTick',dy,'XTickLabel',{},'XLim',[1,Nsteps])
% hold on; plot(zeros(Nsteps,1),':k')
% subaxis(4,1,3,'SpacingHoriz',0.01,'SpacingVert',0.01)
% plot(ept,'k'); ylabel('ept (mm/day)')
% set(gca,'XTick',dy,'XTickLabel',{},'XLim',[1,Nsteps])
% subaxis(4,1,4,'SpacingHoriz',0.01,'SpacingVert',0.01)
% plot(flow,'k'); ylabel('flow (mm/day)')
% set(gca,'XTick',dy,'XTickLabel',ly,'XLim',[1,Nsteps])

