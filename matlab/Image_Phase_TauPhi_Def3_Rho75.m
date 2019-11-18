clear all;

ColMap=[
       0 0 1
       0 1 0
       1 0 0      
       ];


taus=[5, 6, 7, 8, 9, 10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data=zeros( length(taus)*length(phis), 20);
data(1:5*length(phis), :)= load(['/home/xiaoling/SpiralWave/data/Spi10_Phase_TauPhi_Rho75_3.dat']);
data(5*length(phis)+1:end, :)= load(['/home/xiaoling/SpiralWave/data/Spi10_Phase_TauPhi_Rho75_1.dat']);
phase_Pre=zeros(length(taus)*length(phis), 3); % Death; Radial; Break up spiral; Secondary;
phase_Pre2=zeros(length(taus)*length(phis),1);

for i=1:length(taus)*length(phis)
    phase_Pre(i,1)= sum( data(i,:)==0 );
    phase_Pre(i,2)= sum( data(i,:)==1 );  %normal
    phase_Pre(i,3)= sum( data(i,:)==8 ); %spiral
 
    
   
    if (phase_Pre(i,3)>0 )
         
         phase_Pre2(i)=3;
    else        
       phase_Pre2(i)=1* ( phase_Pre(i,1)> phase_Pre(i,2) )+2* ( phase_Pre(i,1)<= phase_Pre(i,2) );  
    end

end

phase=reshape(phase_Pre2, [length(phis), length(taus)]);  %reshape will read the data and then write data in the format of "fortran"
phase=phase';

figure
imagesc(phis, 1:length(taus), phase);
hold on

colormap(ColMap);
%cb =colorbar('Ticks',[1, 2 ,3],'TickLabels',{'D','R' ,'O'}, 'fontsize',20);
%set(cb,'position',[.9 0.5 0.05 0.4])
%c.Label.String = 'phase';  
%   caxis([0,9])



%xlabel('\phi','FontSize',30);
%h=ylabel('\tau_{on}','FontSize',30);
%get(h)
%   h=get(gca,'xlabel')
%   set(h, 'FontSize', 100) 
%   set(h,'FontWeight','bold') %bold font

y=taus'; % x is uniform distributed in log space.
y1=(find(y==10):1:length(y) )'; 
p = polyfit(log(y(y>=10) ),y1,1);
y2=[10, 100, 1000];
y_ticks = polyval(p,log(y2));

%set(gca,'xtick',[find(phis==0.13),find(phis==0.4), find(phis==0.7), find(phis==1)],'xticklabel',{'0.13', '0.4', '0.7', '1'},'FontSize',15,'ytick',y_ticks,'yticklabel', {'10^1','10^2','10^3'}, 'FontSize',15);
set(gca,'ytick',y_ticks,'yticklabel', {'10','10^2','10^3'}, 'FontSize',15); 
set(gca,'ydir','normal')
 
xlabel('$$\phi$$','FontSize',30,'interpreter','latex');
ylabel('$$\tau_{on}$$','FontSize',30,'interpreter','latex');
 
%     set(gca,'xtick',[])
%     set(gca, 'ytick',[])
text( 0.37-0.02, find(taus==452), 'E','FontSize',20, 'FontWeight','bold','color','k' )
text( 0.64-0.02, find(taus==14), 'F', 'FontSize',20,'FontWeight','bold','color','k' )
 
%print(gcf,['../PapFig/Image_Phase_TauPhi_Rho75_Spi10_3.png'],'-dpng','-r500');

