function [vy,vx]=wind(Ni,Nj,vx,vy,t,N1,par)

% par.v=5*10^-2;
V=0;%0.1*par.v;

Vx=par.v;
Vy=3*par.v;
t1=round(N1/1);
% fileID = fopen('vento.txt');%[Estacao, Data, Hora, Dir, Vel] 
% file=fopen('vento2.txt')
% [A] = fscanf(file, '%f, %f, %f, %f, %f', [5 Inf])
% fclose(file)
% A=A';
% A

if (t<=t1)%/1.4)
    for j=1:Nj
        for i=1:Ni
            if(vx(i,j)~=V)
                vx(i,j)=Vx;
%                 vx(i,j)=abs(sin(t));
            end
        
            if(vy(i,j)~=V)
                vy(i,j)=Vy;
%                 vy(i,j)=abs(sin(t));
            end
        end
    end
else
    for j=1:Nj
        for i=1:Ni
            vx(i,j)=V;
            vy(i,j)=V;
        
        end
    end
end