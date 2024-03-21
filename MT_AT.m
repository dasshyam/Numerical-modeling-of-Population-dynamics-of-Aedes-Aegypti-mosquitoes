function [MT,AT]=MT_AT(Mi,Ai,Mxminus,Mxplus,Myminus,Myplus,dx,dy,dt,W,vx,vy,Dx1,Dx2,Dy1,Dy2,i,j,p,h1,h2,par)

                    
                    
                    M(i,j)=Mi;
                    A(i,j)=Ai;
                    M(i-1,j)=Mxminus;
                    M(i+1,j)=Mxplus;
                    M(i,j-1)=Myminus;
                    M(i,j+1)=Myplus;
                    
%                     Dx1=0.5*par.D;  Dx2=Dx1; % coeficiente de difusão na direção x
%                     Dy1=0.5*par.D; Dy2=Dy1; % coeficiente de difusão na direção y
                    
                    
                    if (vx>0)
                        M1=M(i,j);
                        M2=M(i-1,j);
                        else
                        M1=M(i+1,j);
                        M2=M(i,j);
                    end
            
                      
                    if (vy>0)
                        M3=M(i,j);
                        M4=M(i,j-1);
                        else
                        M3=M(i,j+1);
                        M4=M(i,j);
                    end
                    
                    %termos fontes
                    psi1=par.g*Ai*(1-Mi/par.k1)-(par.m1+h1)*Mi; %Fase Alada
                    psi2=par.r*(1-Ai/par.k2)*Mi-(par.m2+par.g+h2)*Ai; %Fase Aquática               
                    
                    Adv(i,j)=dy*(vx*M1-vx*M2)+dx*(vy*M3-vy*M4);                

                    MT(i,j) = M(i,j)+(1/p)*(dt*psi1-W*Adv(i,j)+W*(Dx1*(M(i+1,j)-M(i,j))/dx + Dx2*(-M(i,j)+M(i-1,j))/dx + Dy1*(M(i,j+1)-M(i,j))/dy+Dy2*(-M(i,j)+M(i,j-1))/dy));
                    AT(i,j) = A(i,j)+(dt/p)*psi2; %add 1/p