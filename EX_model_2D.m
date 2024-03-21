%%%%%%%%%%%%%%%Modelo Geral%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%Constantes
par.D=1.25*10^-2; %difusão dispersão aleatoria
par.v=5*10^-2; %advecção vento
par.g=0.2; %gama taxa de maturação: fase aq>>fase alada
par.r=30; %taxa de oviposição
par.k1=25; %capacidade suporte dos mosquitos
par.k2=100; %capacidade suporte da fase aquatica
par.m1=4.0*10^-2; %taxa de mortalidade dos mosquitos
par.m2=1.0*10^-2; %taxa de mortalidade da fase aquatica


Ni=40; %numero de pontos
Nj=Ni;
NN=Ni*Nj;
xi=0; xf=5; yi=0; yf=5; %malha espacial
dx=(xf-xi)/(Ni-1);
X=xi:dx:xf;
dy=(yf-yi)/(Nj-1);
Y=yi:dy:yf;
ti=0; tf=6; %tempo
% dt=(tf-ti)/(10*N);
dt=0.01;
N1=tf/dt;
W=dt/(dx*dy);

%===========POROSIDADE E BLOCOS==========

Adv=zeros(Ni,Nj);

bi=5; %tamanho do bloco na direção x
bj=8;  %tamanho do bloco na direção y

[P,vx1,vy1,Dx1,Dx2,Dy1,Dy2,H1,H2]=P_blocos(Ni,Nj,bi,bj,1,0,0,par);
[vx,vy]=wind(Ni,Nj,vx1,vy1,1,N1,par);

%===========CONDIÇÕES INICIAIS==========

M0=5;  %fase alada
A0=15; %fase aquática

M=zeros(Ni,Nj);
A=zeros(Ni,Nj);
MT=zeros(Ni,Nj);
AT=zeros(Ni,Nj);

for j=1:7
    for i=15:20
        M(i,j)=M0;
        A(i,j)=A0;
    end
end

%===========CONDIÇÕES DE CONTORNO=======

%M_x(xi,y,t)=0=M_x(xf,y,t)
%MC(i-1,j)=MC(i,j);
%M_y(yi,y,t)=0=M_y(yf,y,t)

%%%Plotando as condições iniciais

for j=1:Nj
        for i=1:Ni
            PA(i,j)=P(i,j)*A(i,j);            
            PM(i,j)=P(i,j)*M(i,j);            
        end
end

    [X1 ,Y1]=meshgrid(40*X,40*Y);
    R1=[0:40:200];
    figure(1);
     colormap(jet);
%     [X1 ,Y1]=meshgrid(X,Y);
    surf(X1,Y1,M);
    view(59,28)
    xlabel('x(m)','fontsize',14);
set(gca,'XTick',R1);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',R1);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
c = colorbar('fontsize',14);
c.Label.String = 'Total Population';
%,'Position',ylp,'VerticalAlignment','middle','HorizontalAlignment','right');
%set(ylh,'Rotation',0,'Position',ylp+[-0.3 20 0],'VerticalAlignment','middle','HorizontalAlignment','right');
%      title('Initial Condition (\psi x M)','fontsize',14);
%    colorbar('YTickLabel',{'Low density',[],[],[],[],[],[],[],[],'High density'})
%     set(gca,'YTicklabelMode','manual')
    
% ylh=get(gca,'ylabel');
% gyl=get(ylh);
% ylp=get(ylh,'Position');
% set(ylh,'Rotation',0,'Position',ylp,'VerticalAlignment','middle','HorizontalAlignment','right');
    
    figure(2);
     colormap(jet);
%     [X1 ,Y1]=meshgrid(X,Y);
    surf(X1,Y1,A);
    view(59,28)
    set(gca,'Zlim',[0 60]);
    xlabel('x(m)','fontsize',14);
set(gca,'XTick',[0:40:200]);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',[0:40:200]);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
c = colorbar('fontsize',14);
c.Label.String = 'Total Population';
%,'Position',ylp,'VerticalAlignment','middle','HorizontalAlignment','right');
%     title('Initial Condition (\psi x A)','fontsize',14);
    %colorbar('YTickLabel',{'Low density','High density',A0})
%    colorbar('YTickLabel',{'Low density',[],[],[],[],[],[],[],'High density'})
%     colorbar('YTickLabel',{0,'Low density',A0,'High density'})
   
%Para escolher a condição de fronteira
  FR=0; %neumann dM/dx=0
% FR=1; %dirichlet
% FN=5;
    
       
for n=2:N1
    [vx,vy]=wind(Ni,Nj,vx1,vy1,n,N1,par);
    for j=1:Nj
        for i=1:Ni
           
          if (FR==0)
              
            %condicão de fronteira
            if (i==1)
               M(1,j)=M(2,j);
               %A(1,j)=A(2,j);
            end
            if (j==1)
               M(i,1)=M(i,2);
               %A(i,1)=A(i,2);
            end     
            if (i==Ni)
               M(Ni,j)=M(Ni-1,j);
               %A(N,j)=A(N-1,j);
            end
            if (j==Nj)
               M(i,Nj)=M(i,Nj-1);
               %A(i,N)=A(i,N-1);
            end 
         else
             if (FR==1)
                 %condicão de fronteira
                if (i==1)%mudei aki
%                  if (j>=15 && j<=25) %(i>=15 && i<=25)
%                    M(1,j)=M(2,j)+FN*dx;
                   %A(1,j)=A(2,j);
%                     else
                    M(1,j)=M(2,j);
%                  end
                end
                if (j==1)
                    M(i,1)=0;%M(i,2);
                    %A(i,1)=A(i,2);
                end     
                if (i==Ni)
                    M(Ni,j)=M(Ni-1,j);
                    %A(N,j)=A(N-1,j);
                end
                if (j==Nj)
                    M(i,Nj)=M(i,Nj-1);
                    %A(i,N)=A(i,N-1);
                end 
             end
           end
        
            if ((i~=1)&&(j~=1)&&(i~=Ni)&&(j~=Nj))
                    
                    [Z1,Z2]=MT_AT(M(i,j),A(i,j),M(i-1,j),M(i+1,j),M(i,j-1),M(i,j+1),dx,dy,dt,W,vx(i,j),vy(i,j),Dx1(i,j),Dx2(i,j),Dy1(i,j),Dy2(i,j),i,j,P(i,j),H1(i,j),H2(i,j),par);
                    MT(i,j)=Z1(i,j);
                    AT(i,j)=Z2(i,j);
                    
            end
         end
    end

    M = MT;
    A = AT;

    if (mod(n,300)==0)
%     figure(n);
%     hold off;
% %     colormap(gray);
%         [X1 ,Y1]=meshgrid(X,Y);
%         surf(X1,Y1,M);
%         view(59,28)
% %         set(gca,'Zlim',[0 10]);
%         xlabel('x','fontsize',14);
%         ylabel('y','fontsize',14);
%         zlabel('Population (M)','fontsize',14);
%         colorbar%('YTickLabel',{0,'Low density',[],[],[],[],[],[],[],'High density',20})
% %     set(gca,'YTicklabelMode','manual')
        
%         %title('População no tempo','fontsize',14);
%         
%      figure(n+1);
% %      colormap(gray);
%         [X1 ,Y1]=meshgrid(X,Y);
%         surf(X1,Y1,A);
%         view(59,28)
% %         set(gca,'Zlim',[0 60]);
%         xlabel('x','fontsize',14);
%         ylabel('y','fontsize',14);
%         zlabel('Population (A)','fontsize',14);
%         colorbar%('YTickLabel',{0,'Low density',[],[],[],[],[],'High density',900})
% %     set(gca,'YTicklabelMode','manual')
%         %title('População no tempo t=...','fontsize',14);
%         
        
    for j=1:Nj
        for i=1:Ni
            PA(i,j)=P(i,j)*A(i,j);
            PM(i,j)=P(i,j)*M(i,j);            
        end
    end
    figure(n+2);
%     hold off;
     colormap(jet);
%         [X1 ,Y1]=meshgrid(X,Y);
        surf(X1,Y1,PA);
%         xlabel('x(m)','fontsize',14);
%         ylabel('y(m)','fontsize',14);
%         zlabel('Accumulation Term (\psi)','fontsize',14);
%         title('Total Population (\psi x A)','fontsize',14);
xlabel('x(m)','fontsize',14);
set(gca,'XTick',[0:40:200]);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',[0:40:200]);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
c = colorbar('fontsize',14);
c.Label.String = 'Total Population';
%        colorbar('YTickLabel',{'Low density',[],[],[],[],[],'High density'})
%     set(gca,'YTicklabelMode','manual')
        %saveas(gcf,'fig1.pdf')
%         surf(Y1,X1,M);

%     
%     for j=1:Nj
%         for i=1:Ni
%             PM(i,j)=P(i,j)*M(i,j);            
%         end
%     end
%         

 figure(n+3);
%     hold off;
     colormap(jet);
        %[X1 ,Y1]=meshgrid(X,Y);
        surf(X1,Y1,PM);
%         xlabel('x(m)','fontsize',14);
%         ylabel('y(m)','fontsize',14);
%         zlabel('Accumulation Term (\psi)','fontsize',14);
%         title('Total Population (\psi x M)','fontsize',14);
xlabel('x(m)','fontsize',14);
set(gca,'XTick',[0:40:200]);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',[0:40:200]);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
c = colorbar('fontsize',14);
c.Label.String = 'Total Population';
%colorbar('YTickLabel',{'Low density',[],[],[],[],[],[],[],'High density'})
%     set(gca,'YTicklabelMode','manual')
    end     
end

% ylh=get(gca,'ylabel');
% gyl=get(ylh);
% ylp=get(ylh,'Position');
% set(ylh,'Rotation',0,'Position',ylp,'VerticalAlignment','middle','HorizontalAlignment','right');
             
% figure(n+4)
% %colormap(jet);
% hold on;
% RR=[P(2,2) P(39,39)];
% pie(RR);
% legend('Blocks','Streets');


    figure(n+4)
    hold on;
     colormap(jet);
% axis off;
%         [X1 ,Y1]=meshgrid(X,Y);
        surf(X1,Y1,P);
%        set(gca,'Xticklabel',[]);
%         xlabel('x(m)','fontsize',14);
%         set(gca,'Yticklabel',[]);
%         ylabel('y(m)','fontsize',14);
xlabel('x(m)','fontsize',14);
set(gca,'XTick',[0:40:200]);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',[0:40:200]);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
zlabel('Population area support(\psi)','fontsize',14);
%set(ylh,'Rotation',0,'Position',ylp+[-0.3 0 0],'VerticalAlignment','middle','HorizontalAlignment','right');
%         title('Domain (blocks and streets)','fontsize',14);
%         colorbar
%         legend(P(2,2),'Blocks',P(39,39)'Streets');
        %saveas(gcf,'fig2.pdf')

       
        figure(n+5);
     colormap(jet);
%         [X1 ,Y1]=meshgrid(X,Y);
        quiver(X1,Y1,vy,vx);        
        xlabel('x(m)','fontsize',14);
set(gca,'XTick',[0:40:200]);
ylabel('y(m)','fontsize',14);
set(gca,'YTick',[0:40:200]);
ylh=get(gca,'ylabel');
gyl=get(ylh);
ylp=get(ylh,'Position');
set(ylh,'Rotation',0);
        title('Wind direction','fontsize',14);
