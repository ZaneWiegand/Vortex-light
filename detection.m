clc;
clear;
Lambda=6.328e-7;
wzero=3e-3;
f=0.2;
DT=0.2;
a=(4e-3)/pi;
b=4.5e-3;
L1=5;
L2=1;

D=2e-2;
M=1000;
dx=D/M;
du=1/D;
X=-D/2:dx:D/2-dx;
Y=X;
[x,y]=meshgrid(X,Y);
U=-1/(2*dx):1/D:(1/(2*dx)-1/D);
V=U;
[u,v]=meshgrid(U,V);
r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
XT=U*Lambda*DT;
YT=XT;
[xt,yt]=meshgrid(XT,YT);

FAI1=2*pi*a/(Lambda*f).*(y.*theta-x.*log(sqrt(x.^2+y.^2+eps)/b)+x);
FAI2=-2*pi*a*b/(Lambda*f).*exp(-xt./(a+eps)+eps).*cos(-yt./(a+eps)+eps);
ULA=1*sqrt(2/pi)./wzero.*exp(-r.^2/(wzero.^2)).*sqrt(1/ factorial(abs(L1))).*(r.*sqrt(2)./wzero).^abs(L1).*exp(1i*L1*theta)+1*sqrt(2/pi)./wzero.*exp(-r.^2/(wzero.^2)).*sqrt(1/factorial(abs(L2))).*(r.*sqrt(2)./wzero).^abs(L2).*exp(1i*L2*theta);
ULA1=ULA.*exp(1i*FAI1); 
ULA2=ifftshift(ULA1);
ULAT0=fft2(ULA2);
ULAT=fftshift(ULAT0);

ULB=ULAT.*exp(1i*FAI2);
ULB1=ifftshift(ULB);
ULBT0=fft2(ULB);
ULBT=fftshift(ULBT0);

figure;
subplot(1,3,1)
surf(x,y,abs(ULA))
%colormap('gray')
  shading interp;%插值平滑
    axis equal; %横纵坐标刻度标尺一致
   axis([-(D/2),(D/2),-(D/2),(D/2)]); %取值范围
   view(90,90 )
    box on;
    grid off;
    
    subplot(1,3,2)
surf(x,y,abs(ULAT))
%colormap('gray')
  shading interp;%插值平滑
    axis equal; %横纵坐标刻度标尺一致
   axis([-(D/2),(D/2),-(D/2),(D/2)]); %取值范围
   view(90,90 )
    box on;
    grid off;
    
    subplot(1,3,3)
surf(x,y,abs(ULBT))
%colormap('gray')
  shading interp;%插值平滑
    axis equal; %横纵坐标刻度标尺一致
   axis([-0.00015,0.00015,-0.00015 ,0.00015]) %取值范围
   view(90,90 )
    box on;
    grid off;
    
    hold on;