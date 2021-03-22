function [result]=RectGrating(E0,xp,yp)
[Nx,Ny]=size(E0);
result=zeros(Nx,Ny);
X_begin=ceil((1-xp)*Nx/2);
X_end=ceil((1+xp)*Nx/2);
Y_begin=ceil((1-yp)*Ny/2);
Y_end=ceil((1+yp)*Ny/2);
for a=X_begin:X_end
    for b=Y_begin:Y_end
        result(a,b)=1;
    end
end
end