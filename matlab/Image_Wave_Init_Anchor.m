clear all;
height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;

ColMap=colormap(parula(101));

%as=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8];


t=0:2:400;
tau=16;
tau_off=5;
phi=0.48;
time=[0];
ColMap=zeros(142, 3);
ColMap(:, 2)=linspace(0,1, 142);
ColMap(:, 3)=linspace(0, 1, 142);
ColMap(142,:)=[1,0,0];  %Initiation w=0.2


dd='/home/xiaoling/SpiralWave/data/';                 
matrix=load([dd,'Spi46_bd20t16tf5c10_wave.dat']);  
map=load([dd,'Spi46_bd20t16tf5c10_map.dat']);  
map=map([1:height] +height,:);
for i=1:1:length(time)
    
    figure
    image_time= matrix( find(t==time(i), 1)*height-height+1: 1: find(t==time(i), 1)*height, :);
    
    image_time(map==tau_off)= -0.4; 
    image_time(1:40, 52)=1.03;
     
    imagesc(image_time);     
    colormap(ColMap);
%     colorbar('Ticks',[0, 1],...
%          'TickLabels',{'non-firing','firing'})
    caxis([-0.41, 1.03])
%    c.Label.String = 'phase';  
 
%    xlabel('\phi','fontsize',15);
%    ylabel('\tau','fontsize',15);
   % set(gca,'ydir','normal');
    set(gca,'xtick',[])
    set(gca, 'ytick',[])
    
    
    title(['t=', num2str(time(i))], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    
    text(53, 30, '$$\leftarrow v=0.2$$',  'Color',[0.9,0.9,0.9],'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    text(28, 30, '$$u=1 \rightarrow$$',  'Color',[0.9,0.9,0.9],'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    text(60, 60, 'off cells',  'Color',[0.9,0.9,0.9],'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %text(80, 8, ['t=', num2str(time(i))], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi46_bd20t16tf5c10_Ini',num2str( time(i) ),'.png'],'-dpng','-r400');
   
    
    
end

savefig(gcf,'../PapFig/Spi46_bd20t16tf5c10_Ini.fig')

