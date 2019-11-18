clear all;
height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;


%as=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8];


t=0:2:400;
tau=16;
tau_off=5;
phi=0.48;
time=[0];
ColMap=[
       1 1 1;
       0 0 0;
       0 1 1;
       1 0 0];


dd='/home/xiaoling/SpiralWave/data/';                 
matrix=load([dd,'Spi46_bd20t16tf5c10_wave.dat']);  
map=load([dd,'Spi46_bd20t16tf5c10_map.dat']);  
map=map([1:height] +height,:);
for i=1:1:length(time)
    
    figure
    image_time= ones(height, width);
    
    image_time(map==tau_off)= 2; 
    image_time(1:40, 51)=3;
    image_time(1:40, 52)=4;
     
    imagesc(image_time);     
    colormap(ColMap);
%     colorbar('Ticks',[0, 1],...
%          'TickLabels',{'non-firing','firing'})
    caxis([1, 4])
%    c.Label.String = 'phase';  
 
%    xlabel('\phi','fontsize',15);
%    ylabel('\tau','fontsize',15);
   % set(gca,'ydir','normal');
    set(gca,'xtick',[])
    set(gca, 'ytick',[])
    
    
    title(['t=', num2str(time(i))], 'FontSize',25, 'FontWeight', 'bold','interpreter', 'latex')
    
    text(53, 30, '$$\leftarrow v=0.2$$',  'Color','k','FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    text(28, 30, '$$u=1 \rightarrow$$',  'Color','k','FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    text(61, 50, 'off cells',  'Color','k','FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    text(43, 80, 'on cells',  'Color','k','FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %text(80, 8, ['t=', num2str(time(i))], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi46_bd20t16tf5c10_Ini',num2str( time(i) ),'.png'],'-dpng','-r400');
   
    
    
end

savefig(gcf,'../PapFig/Spi46_bd20t16tf5c10_Ini.fig')

