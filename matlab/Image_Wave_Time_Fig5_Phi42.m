height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;

ColMap=colormap(parula(101));

%as=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8];


t=0:2:400;
tau=240;

repeat=11; %Which spiral. 
times= 254: -6: 188;


ColMap=zeros(141, 3);
ColMap(:, 2)=linspace(0,1, 141);
ColMap(:, 3)=linspace(0, 1, 141);

matrix=matrix_Pre([1:length(t)*height]+ (repeat-1)*length(t)*height ,:);

nvth=nvth_Pre([1:height]+(repeat-1)*height, :);
dlmwrite(['../data/Spi7_Rho0t240p420c10R', num2str(repeat), '.dat'],nvth,'delimiter','\t','precision',3)


% dd='/home/xiaoling/SpiralWave/data/';   
% matrix_Pre=load([dd,'Spi7_R40Rho0t240p420c10_wave.dat']);  
% nvth_Pre=load([dd,'Spi7_R40Rho0t240p420c10_map.dat']); 


for i=1:1:length(times)
    
    time=times(i);
    figure
    image_time= matrix( find(t==time, 1)*height-height+[1: 1: height], :);
     
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
    
    text(80, 8, ['t=', num2str(time)], 'FontSize',20, 'FontWeight', 'bold','interpreter', 'latex')
    %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../Fig5/Spi7_R',num2str(repeat),'Rho0t240p420c10_T',num2str( time ),'.png'],'-dpng','-r400');
    
    
end

