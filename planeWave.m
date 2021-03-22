clear;
N=100;
lambda=6.328e-7;
k=2*pi/lambda;    %波数
x=linspace(-10,10,N);
y=linspace(-10,10,N);
[X,Y]=meshgrid(x,y);
[theta,r]=cart2pol(X,Y);
s=1;
for L=-2:0.5:2
    subplot(3,3,s)
    E1=exp(-1i*k*X);
    E2=exp(1i*L*theta);
    E=E1+E2;
    I=E.*conj(E);
    pcolor(X,Y,I);
    axis square;
    title(['L = ',num2str(L)]);
    s=s+1;
end