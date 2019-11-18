clear all;
h=100;
w=100;

ColMap=[0.9,0.9,0.9];
 
map=load('/home/xiaoling/SpiralWave/data/Spi7_n6_2Rho0t240p500c10_map.dat');
map=map([1:h]+h,:);
map1=map(75:78, 81:84);
data=ones(100, 100);
% [MatX, MatY]=meshgrid(1: size(data,2), 1:size(data, 1));
% MatR=sqrt( (MatX- (w+1)/2).^2 + (MatY- (h+1)/2).^2 );
% data( abs(MatR-20)<1.5 )=2;

 figure()
    imagesc(data);
    hold on
    colormap(ColMap);
    set(gca, 'XTickLabels',[], 'YTickLabels',[])
    
    dim = [.32 .32 .4 .4];
    annotation('ellipse',dim,'Color','c','LineWidth',7)
    
    plot(50.5, 50.5, 'O', 'MarkerSize', 10, 'MarkerFaceColor','c', 'MarkerEdgeColor','c')
    text(42, 42, 'trigger', 'FontSize', 20)
    
    text(34, 15, 'wave propagation', 'FontSize', 20)

    
%     annotation('textarrow', [0.7, 0.8]-0.03,[0.7, 0.8]-0.03, 'Color', 'g', 'LineWidth',3);
%     annotation('textarrow', [0.3, 0.2]+0.07,[0.3, 0.2]+0.07, 'Color', 'g', 'LineWidth',3);
%      annotation('textarrow', [0.3, 0.2]+0.05,[0.7, 0.8]-0.05, 'Color', 'g', 'LineWidth',3);
%       annotation('textarrow', [0.7, 0.8]-0.05,[0.3, 0.2]+0.05, 'Color', 'g', 'LineWidth',3);
    %annotation('textarrow', [0.3, 0.3, -0.1, -0.1], 'Color', 'g', 'LineWidth',7,'TextRotation',45);
    
    
    plot([80.5, 84.5, 84.5, 80.5, 80.5], [80.5, 80.5, 84.5, 84.5, 80.5]-6, '-k','LineWidth', 3)
    
    
     %title(name(i),'position',[1 1]);
     
    print(gcf,['../PapFig/Cartoon_WF_BioFilm_Fig0.png'],'-dpng','-r400');
    
%     figure()
%     ColMap=[0, 0, 0;
%             1, 1 , 1];
%         
%     for x=1:size(map1, 2)
%         for y=1:size(map1, 1)
%             fill([x, x+1, x+1, x], [y, y, y+1, y+1], ColMap( (map1(y,x)==240)+1,:),'Edgecolor',[0.5, 0.5, 0.5],'LineWidth',3);  
%             hold on
%         end
%     end
%         
%     colormap(ColMap);
%    set(gca, 'XTickLabels',[], 'YTickLabels',[])
%     print(gcf,['../PapFig/Cartoon_BioFilm_BW_Fig0.png'],'-dpng','-r400');
%     

