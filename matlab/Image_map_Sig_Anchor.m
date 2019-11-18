clear all;
height=100;
width=100;
ColMap=[
        0, 0, 0;
        0,1,1
        ];
data= load(['/home/xiaoling/SpiralWave/data/Spi7_n6_2Rho0t240p580c10_map.dat']);
data=data(height+1: height+height,:);
figure
imagesc(data);
hold on

plot([40, 40 ,40, 30, 20, 10, 20, 20, 30, 40],[70, 80, 90, 90, 90, 80, 85, 80, 75, 70], 'linewidth',3,'Color','m')
colormap(ColMap);
%c =colorbar('Ticks',[5, 240],'TickLabels',{'off','on'});
%c.Label.String = 'phase';  
%   caxis([0,9])
%caxis([-0.4, 1])
%c.Label.String = 'phase';  
 
%    xlabel('\phi','fontsize',15);
%    ylabel('\tau','fontsize',15);
  
    set(gca,'xtick',[])
    set(gca, 'ytick',[])
    
    % set(gca,'ydir','normal')
    %set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);
print(gcf,['../PapFig/Image_Spi7_n6_2Rho0t240p580c10_F4.png'],'-dpng','-r400');
