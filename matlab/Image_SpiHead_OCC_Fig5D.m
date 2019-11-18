clear all;

height=100;
width=100;
ColMap=[1,1,1; 0 0 0]
OCC=load('/home/xiaoling/SpiralWave/data/Spi7_n5Rho0t240p420c10_OCCin.dat');
Y=floor(OCC/width);
X=rem(OCC, width);

Y=Y+1;
X=X+1;

SpiHead=[50    44    45    45    44    48    48  51  53  50;
         18    19    27    39    49    47    41  31  23   18];
     
figure;
biofilm=zeros(height, width);

for i=1:length(X)
    biofilm(Y(i), X(i))=1;
end

imagesc(biofilm);
hold on;
colormap(ColMap);

plot(SpiHead(1,:), SpiHead(2, :), '-O', 'color', 'c', 'LineWidth', 4, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
title('$$\phi=0.42$$', 'fontsize',25,'interpreter','latex','FontWeight', 'bold')
 set(gca,'xtick',[])
 set(gca, 'ytick',[])
   %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   %print(gcf,['../PapFig/Spi7_n5Rho0t240p420c10_OCC.png'],'-dpng','-r500');