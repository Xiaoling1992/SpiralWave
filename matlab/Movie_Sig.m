clear all;
height=100;
width=100;

window_s=1;
window_e=height;
window_h=window_e-window_s+1;

ColMap=colormap(parula(101));

%as=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8];


t=0:2:300;
tau=240;
phi=0.42;
ColMap=zeros(141,3);
ColMap(1:141,2)=linspace(0,1,141);
ColMap(1:141,3)=linspace(0,1,141);

dd='/home/xiaoling/SpiralWave/data/';                  
           % fprintf('a=%f bH=%f bL=%f c=%f\n',a,bH,bL, c);
 matrix=load([dd,'Spi10Rho75Phi640t14c10_wave_1.dat']);  
 %matrix=load([dd,'Fit12Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'.dat']);            
%matrix=load([dd,'Fit12Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'K',num2str(K*1000),'dt',num2str(dt*100),'.dat']);
%nvth1=load([dd,'Fit10Phi50Rho100a',num2str(a*100),'b',num2str(b*10000),'c',num2str(c*10),'cL',num2str(cL*10),'Mp.dat']);
%nvth=nvth1(2*height+1:3*height,:);
%matrix=load('/home/xiaoling/Fitzhugh/data/Fit8a1l10b8.dat');
%nvth=load('/home/xiaoling/Fitzhugh/data/Fit8a1l10b8Mp.dat');


%order=1:1:length(t);

vidobj = VideoWriter('../PapFig/Spi10Rho75Phi640t14c10_1.avi','Motion JPEG AVI');
vidobj.FrameRate = 10;
vidobj.Quality = 75;
open(vidobj);
%a=max(max(data))
%b=min(min(data))
for it=1:1:length(t)
   wavepeak=matrix((it-1)*height+window_s:(it-1)*height+window_e,1:width);
   imagesc(1:width,1:window_h,wavepeak);
   colormap(ColMap);
   caxis([-0.4,1]);
   colorbar;

   t_str=['t=' num2str(t(it))]; %Generate the time string
    title(['Time: ' t_str, ' \tau=14', ' \phi=0.64']);
    fr = getframe(gcf);
    writeVideo(vidobj,fr);
end
close(vidobj)


% %gifname=['../figure/','PC1Phi',num2str(phi*100),'Rho',num2str(rho*100),'a',num2str(a*100),'b',num2str(bH*10000),'bL',num2str(bL*10000),'c',num2str(c*10),'K',num2str(K*1000),'dt',num2str(dt*100),'.gif'];
%  gifname=['../PapFig/Spi423_bRho0u2200s40c10.mp4'];
% 
% figure(1000)
%      
% for it=1:1:length(t)
% 
% wavepeak=matrix((it-1)*height+window_s:(it-1)*height+window_e,1:width);
% imagesc(1:width,1:window_h,wavepeak);
% caxis([-0.4,1]);
% colorbar;
% 
% t_str=['t=' num2str(t(it))]; %Generate the time string
% title(['Time: ' t_str, ' \tau=22', ' \phi=0.40']);
% 
%       drawnow
%       frame = getframe(1000);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%        if it == 1;
%           imwrite(imind,cm,gifname,'', 'Loopcount',inf);
%        else
%           imwrite(imind,cm,gifname,'','WriteMode','append');
%       end
% end

       


      
    











