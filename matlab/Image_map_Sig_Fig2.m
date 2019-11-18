clear all;
height=100;
width=100;
ColMap=[
        0, 0, 0;
        1,1,1
        ];
data= load(['/home/xiaoling/SpiralWave/data/Spi7_bRho0t530p46c10_map.dat']);
data=data(height+1: height+height,:);
figure
imagesc(data);
hold on
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
    %title('$$\rho=0$$','fontsize',35, 'FontWeight','bold','interpreter','latex' )
    
    % set(gca,'ydir','normal')
    %set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);
print(gcf,['../PapFig/Map_Spi7_bRho0t530p46c10.png'],'-dpng','-r400');
