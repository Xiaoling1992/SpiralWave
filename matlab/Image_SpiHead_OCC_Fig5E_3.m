%Mark different off cell clusters with different colors. 
clear all;

height=100;
width=100;
ColStep=1001;
ColMap=[1 1 1;  %Environment
        1 0 0;
        1 1 0;
       
        0  1 0;
       255/255,165/255,0/255;
        0  0  1;
        
        1  0  1;
        ]

% ColMap(2:end,1)=linspace(0.8, 0, ColStep-1 );
% ColMap(2:end,2)=linspace(0.8, 0, ColStep-1);
% ColMap(2:end,3)=linspace(0.8, 0, ColStep-1 );
OCC=load('/home/xiaoling/SpiralWave/data/Spi7_n5Rho0t240p540c10_OCCin_2.dat');
Y=floor(OCC(:,1)/width);
X=rem(OCC(:,1), width);
Y=Y+1;
X=X+1;

ColCod=OCC(:,2);     %the No. of cluster in the whole biofilm;
CluMin=min(ColCod);
CluMax=max(ColCod);
ColCod=ColCod- min(ColCod); 
NewClu=cumsum ((diff(ColCod) >0 ));  %the No. of cluster in the picked up clusters
ColCod(2:end)=NewClu;
ColCod=rem(ColCod, 6)+1;

% k=(1000-1)/(CluMax-CluMin);
% b=1- k*CluMin;
% ColCod= round ( (ColCod*k+b), 0); 

%Define the number of cluster: mark the cells in the same cluster with the
%same number, which increases. 



SpiHead=[38  43    47   40     39   31     24    20   25     31    39  
;    77  79    78   72      67  68     71    76   79     81    80  
];
     
figure;
biofilm=zeros(height, width);

for i=1:length(X)
    biofilm(Y(i), X(i))= ColCod(i);
end

imagesc(biofilm);
hold on;
colormap(ColMap);

plot(SpiHead(1,:), SpiHead(2, :), '-O', 'color', 'c', 'LineWidth', 4, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
title('$$\phi=0.54$$', 'fontsize',25,'interpreter','latex','FontWeight', 'bold')

 set(gca,'xtick',[])
 set(gca, 'ytick',[])
   %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi7_n5Rho0t240p540c10_OCC_2.png'],'-dpng','-r500');