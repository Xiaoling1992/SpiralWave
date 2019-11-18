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
time=[4, 8, 12, 18];
ColMap=zeros(141, 3);
ColMap(:, 2)=linspace(0,1, 141);
ColMap(:, 3)=linspace(0, 1, 141);


dd='/home/xiaoling/SpiralWave/data/';                 
matrix=load([dd,'Spi423_bRho0u2200s57c10_wave.dat']);  

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
    text(75, 10, ['t=', num2str(time(i))], 'FontSize',30, 'Color', 'w')
    %text(80, 8, ['t=', num2str(time(i))], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi423_bRho0u2200s57c10_T',num2str( time(i) ),'.png'],'-dpng','-r400');
    
    
end

