function [P,vy,vx,Dx1,Dx2,Dy1,Dy2,H1,H2]=JF_blocos(Ni,Nj,bi,bj,p,h1,h2,par)
%Ni numero de elementos em x
%Nj numero de elementos em y
%bi tamanho do bloco em x
%bj tamanho do bloco em y

%p porosidade
phi=0.3*ones(Ni,Nj); %ser� o valor na rua
xi=1; %x inicial
yi=1; %x final
xf=1+bi;
yf=1+bj;
% (xi<=Ni && yi<=Nj)

%velocidade
Vx=0;%.1*par.v; % velocidade na dire��o x
Vy=0;%.1*par.v; % velocidade na dire��o y

%ser� o valor na rua
vx=par.v*ones(Ni,Nj); % velocidade na dire��o x
vy=3*par.v*ones(Ni,Nj); % velocidade na dire��o y 


%Difus�o
DDx1=0.3*par.D;  DDx2=DDx1; % coeficiente de difus�o na dire��o x
DDy1=0.3*par.D; DDy2=DDy1; % coeficiente de difus�o na dire��o y

%ser� o valor na rua
Dx1=1*par.D*ones(Ni,Nj);  Dx2=Dx1; % coeficiente de difus�o na dire��o x
Dy1=1*par.D*ones(Ni,Nj); Dy2=Dy1; % coeficiente de difus�o na dire��o y
                    

%O Veneno ser� dado por h1 e h2
H1=zeros(Ni,Nj);
H2=zeros(Ni,Nj);

%cria os blocos
while (yi<=Nj)
    
while (xi<=Ni)

for j=yi:yf
    for i=xi:xf
    
        phi(i,j)=p;
        
        vx(i,j)=Vx;
        vy(i,j)=Vy;
        
        Dx1(i,j)=DDx1;
        Dx2(i,j)=DDx2;
        Dy1(i,j)=DDy1;
        Dy2(i,j)=DDy2;
        H1(i,j)=h1;
        H2(i,j)=h2;
        
    end

end
xi=xf+bi;
xf=xf+2*bi;
    if xf>Ni
       xf=Ni;
    end
end
xi=1;
xf=1+bi;

yi=yf+bj;
yf=yf+2*bj;
    
    if yf>Nj
       yf=Nj;
    end
       

end

P=phi;