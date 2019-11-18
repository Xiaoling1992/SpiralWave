clear all;

ColMap=[
       0 0 1
       0 1 0
       1 0 0      
       ];

% ColMap=[
%        0 1 0
%        1 0 0      
%        ];


EpsOffs=[0.003, 0.004, 0.005, 0.006, 0.008, 0.011, 0.014, 0.018, 0.023, 0.030, 0.039, 0.050, 0.065, ...
	0.083, 0.108, 0.139, 0.180, 0.232, 0.300, 0.387, 0.500, 0.646, 0.834, 1.077, 1.391, 1.797, 2.321, 2.998, 3.871, 5.000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data= load(['/home/xiaoling/SpiralWave/data/Spi7_b_Phase_Rho0_tau240.dat']);
phase_Pre=zeros(length(EpsOffs)*length(phis), 3); % Death; Radial; Break up spiral; Secondary;
phase_Pre2=zeros(length(EpsOffs)*length(phis),1);

for i=1:length(EpsOffs)*length(phis)
    phase_Pre(i,1)= sum( data(i,:)==0 );
    phase_Pre(i,2)= sum( data(i,:)==3 )+sum( data(i,:)==5 );  %Break up and normal
    phase_Pre(i,3)= sum( data(i,:)==2 )+sum( data(i,:)==8 )+sum( data(i,:)==9 ); %other: BS, SS, and secondary
 
    
   
    if (phase_Pre(i,3)>0 )
         
         phase_Pre2(i)=3;
    else        
        [ tem,phase_Pre2(i)]=max( phase_Pre(i, :) );
    end

end

phase=reshape(phase_Pre2, [length(phis), length(EpsOffs)]);
phase=phase';

figure
imagesc(phis, 1:length(EpsOffs),phase);
hold on

colormap(ColMap);
%cb =colorbar('Ticks',[1, 2 ,3],'TickLabels',{'D','R' ,'O'}, 'fontsize',20);
%set(cb,'position',[.9 0.5 0.05 0.4])
%c.Label.String = 'phase';  
%   caxis([0,9])

xlabel('$$\phi$$','fontsize',25,'FontWeight', 'bold','interpreter','latex');
ylabel('$$\epsilon_{off}$$','fontsize',25,'FontWeight', 'bold','interpreter','latex');
%set(gca,'xtick',[1:3: length(phis)],'xticklabel',phis([1:3: length(phis)]),'ytick',1:3:length(EpsOffs),'yticklabel',EpsOffs(1:3:length(EpsOffs) ));
y=(1:1:length(EpsOffs) )';
x=EpsOffs';
p = polyfit(log(x),y,1);
x1=[0.01, 0.1, 1];
y1 = polyval(p,log(x1));
set(gca,'ytick',y1,'yticklabel', {'10^{-2}','10^{-1}','10^{0}'});


 
%     set(gca,'xtick',[])
%     set(gca, 'ytick',[])
    
 set(gca,'ydir','normal')
 
% plot(find(phis==0.46), find(EpsOffs==530), 'kX', 'LineWidth', 3, 'markersize' ,15,'MarkerEdgeColor','k' )
%  plot(find(phis==0.49), find(EpsOffs==19), 'k*', 'LineWidth',3,'markersize' ,15,'MarkerEdgeColor','k' )
 
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);
print(gcf,['../PapFig/Image_Phase_EpsoffPhi_Tau240_Other.png'],'-dpng','-r400');






