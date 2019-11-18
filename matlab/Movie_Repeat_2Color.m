clear all;
height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;

ColMap=zeros(282,3);
ColMap(1:141,3)=linspace(0.1,1,141); %off cells
ColMap(142:282,1)=linspace(0.1,1,141); %on cells
ColMap(142:282,2)=linspace(0.1,1,141); %on cells

a=0;
c=10;
rho=0;

tau=240;
tau_off=5;
phis=0.32: 0.02 : 0.58;
t=[0:2:400];


dd='/home/xiaoling/SpiralWave/data/';
matrix_Pre=load([dd,'Spi7_R40Rho0t240p420c10_wave.dat']);  
nvth_Pre=load([dd,'Spi7_R40Rho0t240p420c10_map.dat']); 
       
          
for i=0:1:length( nvth_Pre(:,1) )/height-1;           
           % fprintf('a=%f bH=%f bL=%f c=%f\n',a,bH,bL, c);
 
  matrix=matrix_Pre( [1: length(t)*height] +i*length(t)*height, : );
  nvth=nvth_Pre([1:height]+i*height, :);
  dlmwrite(['../data/Spi7_Rho0t240p420c10R', num2str(i+1), '.dat'],nvth,'delimiter','\t','precision',3)
%matrix=load([dd,'Fit12Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'.dat']);            
%matrix=load([dd,'Fit12Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'K',num2str(K*1000),'dt',num2str(dt*100),'.dat']);
%nvth1=load([dd,'Fit10Phi50Rho100a',num2str(a*100),'b',num2str(b*10000),'c',num2str(c*10),'cL',num2str(cL*10),'Mp.dat']);
%nvth=nvth1(2*height+1:3*height,:);
%matrix=load('/home/xiaoling/Fitzhugh/data/Fit8a1l10b8.dat');

 Tau01=(nvth~=tau_off);
%order=1:1:length(t);


%gifname=['../figure/','PC1Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'K',num2str(K*1000),'dt',num2str(dt*100),'.gif'];
figure(1000)

vidobj = VideoWriter(['../PapFig/Spi7_Rho0t240p420c10R', num2str(i+1), '.avi'],'Motion JPEG AVI');
vidobj.FrameRate = 10;
vidobj.Quality = 75;

open(vidobj);
%a=max(max(data))
%b=min(min(data))
for it=1:1:length(t)
   wavepeak=matrix((it-1)*height+window_s:(it-1)*height+window_e,1:width);
   wavepeak=wavepeak+1.41*Tau01;
   imagesc(1:width,1:window_h,wavepeak);
   colormap(ColMap)


   t_str=['t=' num2str(t(it))]; %Generate the time string
   title(['Time: ' t_str, ' \tau=240', ' \phi=0.42']);
    fr = getframe(gcf);
    writeVideo(vidobj,fr);
end
close(vidobj)

end

%      
% for it=1:1:length(t)
% 
% wavepeak=matrix((it-1)*height+window_s:(it-1)*height+window_e,1:width);
% wavepeak=wavepeak+1.41*Tau01;
% 
% imagesc(1:width,1:window_h,wavepeak);
% caxis=[-0.4, 2.82]
% clim=[-0.4, 2.82]
% 
% 
% %clim([-0.4,1]);
% colorbar;
% 
% t_str=['t=' num2str(t(it))]; %Generate the time string
% title(['Time: ' t_str,' a=', num2str(a), ' tau_{on}=',num2str(tau)]);
% 
%       drawnow
%       frame = getframe(1000);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%        if it == 1;
%           imwrite(imind,cm,gifname,'gif', 'Loopcount',inf);
%        else
%           imwrite(imind,cm,gifname,'gif','WriteMode','append');
%       end
% end



    


      
    











