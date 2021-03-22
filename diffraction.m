clear;
N=100;
lambda=6.328e-7;
k=2*pi/lambda;
w0=0.5e-3;
[x0,y0]=meshgrid(linspace(-3*w0,3*w0,N));
[theta,r] = cart2pol(x0,y0);
E=zeros(N);
[x,y]=meshgrid(linspace(-3*w0,3*w0,N));
s=1;
for L=-4:1:4
     subplot(3,3,s)
     E0=exp(1i*L*theta);
     %矩孔衍射
     Z=0.1;
     t=RectGrating(E0,0.1,0.1);
     %双缝衍射
     %Z=0.25;
     %t=DoubleslitGrating(E0,0.1);
     for a=1:N
            for b=1:N
                E(a,b)=exp(1i*k*Z)/1i/lambda/Z*sum(sum(E0.*t.*exp(1i*k/2/Z*((x(a,b)-x0).^2+(y(a,b)-y0).^2))));
            end
     end
     I=E.*conj(E);
     pcolor(x,y,I);
     axis square;
     title(['L = ',num2str(L)]);
     s=s+1;
end