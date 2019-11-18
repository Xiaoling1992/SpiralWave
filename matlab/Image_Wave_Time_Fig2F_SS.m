clear all;
height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;

ColMap=colormap(parula(101));

%as=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8];


t=0:2:400;
tau=240;
phi=0.48;
time=[130, 138, 146, 152 ];
ColMap=zeros(141, 3);
ColMap(:, 2)=linspace(0,1, 141);
ColMap(:, 3)=linspace(0, 1, 141);


dd='/home/xiaoling/SpiralWave/data/';                 
matrix=load([dd,'Spi10Rho75Phi640t14c10_wave_1.dat']);  

for i=1:1:length(time)
    
    figure
    image_time= matrix( find(t==time(i), 1)*height-height+1: 1: find(t==time(i), 1)*height, :);
     
    imagesc(image_time);
    colormap(ColMap);  
%     colorbar('Ticks',[0, 1],...
%          'TickLabels',{'non-firing','firing'})
    caxis([-0.4, 1])
%    c.Label.String = 'phase';  
 
%    xlabel('\phi','fontsize',15);
%    ylabel('\tau','fontsize',15);
   % set(gca,'ydir','normal');
    set(gca,'xtick',[])
    set(gca, 'ytick',[])
    text(60, 10, ['t=', num2str(time(i))], 'FontSize',40, 'Color', 'w')
    
%      if(i==1)
%        annotation('textarrow', [0.2, 0.3],[0.5, 0.5], 'Color', 'r', 'LineWidth',3); 
%     elseif(i==2)
%        annotation('textarrow', [0.4, 0.3]+0.05,[0.4, 0.4], 'Color', 'r', 'LineWidth',3); 
%     elseif(i==3)
%         annotation('textarrow', [0.2, 0.3], [0.5, 0.5]+0.05, 'Color', 'r', 'LineWidth',3); 
%     else
%         annotation('textarrow', [0.4, 0.3]+0.05,[0.3 0.3]+0.05, 'Color', 'r', 'LineWidth',3); 
%     end
    %text(80, 8, ['t=', num2str(time(i))], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi10Rho75Phi640t14c10_T',num2str( time(i) ),'.png'],'-dpng','-r500');
    
    
end

