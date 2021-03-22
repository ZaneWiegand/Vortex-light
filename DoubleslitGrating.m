function [result]=DoubleslitGrating(E0,p)
[Nx,Ny]=size(E0);
result=zeros(Nx,Ny);
W=round(p*Nx/2);
result(Nx/2+W,:)=1;
result(Nx/2-W,:)=1;
end