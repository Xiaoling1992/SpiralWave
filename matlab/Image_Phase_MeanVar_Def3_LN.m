clear all;

ColMap=[
       0 0 1
       0 1 0
       1 0 0      
       ];

% ColMap=[
%        0 1 0
%        1 0 0      
%        ];


means=[6, 10, 17, 28, 46, 77, 129, 216, 360, 600, 1000, 1669, 2783, 4643, 7744, 12918, 21547, 35941, 59951, 100000];
variances=[2, 4, 6, 11, 20, 34, 61, 108, 190, 336, 595, 1051, 1857, 3282, 5800, 10250, 18116, 32016, 56583, 100000];

data= load(['/home/xiaoling/SpiralWave/data/Spi41_LN_Phase.dat']);
phase_Pre=data(1:length(means),:);
phase_Pre(phase_Pre==0)   =0;   %death
phase_Pre( phase_Pre==3 )=1;    %break->radial
phase_Pre( (phase_Pre==2) )=2;  %BS-> spiral
phase_Pre(  (phase_Pre==8)  )=2; %SS-> spiral
phase_Pre( (phase_Pre==9) )=1;  %secondary -> radial  
phase_Pre( phase_Pre==5 )=1;  %radial -> radial


figure
imagesc(phase_Pre);
hold on

colormap(ColMap);
%cb =colorbar('Ticks',[1, 2 ,3],'TickLabels',{'D','R' ,'O'}, 'fontsize',20);
%set(cb,'position',[.9 0.5 0.05 0.4])
%c.Label.String = 'phase';  
%   caxis([0,9])


%set(gca,'xtick',[1:3: length(phis)],'xticklabel',phis([1:3: length(phis)]),'ytick',1:3:length(EpsOffs),'yticklabel',EpsOffs(1:3:length(EpsOffs) ));
y1=(1:1:length(means) )';
x1=means';
p = polyfit(log(x1),y1,1);
x1=[10, 1000, 100000];
y1 = polyval(p,log(x1));

y2=(1:1:length(variances) )';
x2=variances';
p = polyfit(log(x2),y2,1);
x2=[10, 1000, 100000];
y2 = polyval(p,log(x2));
set(gca,'xtick',y1,'xticklabel', {'10^1','10^3','10^5'},'FontSize',15, 'ytick',y2,'yticklabel', {'10^1','10^3','10^5'},'FontSize',15);
xlabel('standard deviation','fontsize',25);
ylabel('mean','fontsize',25);

 
%     set(gca,'xtick',[])
%     set(gca, 'ytick',[])
    
 set(gca,'ydir','normal')
 
% plot(find(phis==0.46), find(EpsOffs==530), 'kX', 'LineWidth', 3, 'markersize' ,15,'MarkerEdgeColor','k' )
%  plot(find(phis==0.49), find(EpsOffs==19), 'k*', 'LineWidth',3,'markersize' ,15,'MarkerEdgeColor','k' )
 
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);
print(gcf,['../PapFig/Image_Phase_MeanVariance_LN_Def3.png'],'-dpng','-r400');






