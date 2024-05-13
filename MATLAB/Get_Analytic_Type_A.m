function [z_analytic,phi_analytic]=Get_Analytic_Type_A(Zmin,Zmax,angle)
%% Get the data of z and Potential



%% Input



zl=Zmax-Zmin;
Tph=2.2; % eV
Te=12; % eV
lambda=1.378; % m
lr=zl*lambda; % physical length, m


%% Analytic

p=typeA_in2(Tph,Te,angle);
phi=fsolve(@(t)funs(t,p),[10,-1]);
p=Updatep(phi,p);% update the parameter set p since some parameters are obtained through fsolve

phi_barrier=SolveZ(p,p.n_e_inf,p.phi_0,p.phi_m)+Zmin*lambda;

solinitB=bvpinit(linspace(Zmin*lambda,phi_barrier,200),[1.0,-0.01]);
solB=bvp4c(@(t,y)odefcnB(t,y,p),@(ya,yb)bcfcnB(ya,yb,p),solinitB);
xB=linspace(Zmin*lambda,phi_barrier,200);
phi_tempB=deval(solB,xB);
phiB = phi_tempB(1,:);

solinitC=bvpinit(linspace(phi_barrier,Zmax*lambda,200),[1.0,-0.01]);
solC=bvp4c(@(t,y)odefcnC(t,y,p),@(ya,yb)bcfcnC(ya,yb,p),solinitC);
xC=linspace(phi_barrier,Zmax*lambda,200);
phi_tempC=deval(solC,xC);
phiC = phi_tempC(1,:);

phiAnalytic=[phiB,phiC];
xBC=[xB,xC];

%% Return

z_analytic=xBC/lambda;
phi_analytic=phiAnalytic/Tph;

% plot(phi_analytic,z_analytic)




%% bcfcnB

function res=bcfcnB(ya,yb,p)

res=[ya(1)-(p.phi_0); yb(1)-p.phi_m];

end


%% bcfcnC

function res=bcfcnC(ya,yb,p)

res=[ya(1)-p.phi_m; yb(1)-0];

end



%% funs

function F=funs(phi,p)
phi_0=phi(:,1);
phi_m=phi(:,2);
a=1-2*erf(-p.V_d/p.Vth_e)+erf(real(sqrt(-phi_m/p.Te))-p.V_d/p.Vth_e);

b=p.n_p_0*0.5*exp(-phi_0/p.Tp).*(1-erf(real(sqrt(-phi_m/p.Tp))))-p.n_i_inf;
n_e_inf=-b./a*2;

E_e_f_inf = n_e_inf/p.n_p_0_90d*p.Te/p.Tp*(( 1 - erf(real(sqrt((-phi_m)/p.Te))-p.V_d/p.Vth_e))...
          - exp((phi_m)/p.Te).*(1-erf(-p.V_d/p.Vth_e))...
          +1/sqrt(pi)*p.Vth_e/p.V_d*exp((phi_m)/p.Te-p.V_d^2/p.Vth_e^2)...
          .*(exp(2*p.V_d/p.Vth_e*real(sqrt((-phi_m)/p.Te)))-1)...
                            );                        

E_e_r_inf=2*n_e_inf/p.n_p_0_90d*p.Te/p.Tp*(( erf(real(sqrt((-phi_m)/p.Te))-p.V_d/p.Vth_e)+erf(p.Vd/p.Vth_e))...
          -1/sqrt(pi)*p.Vth_e/p.V_d*exp(phi_m/p.Te-p.V_d^2/p.Vth_e^2).*(exp(2*p.V_d/p.Vth_e*real(sqrt((-phi_m)/p.Te)))-1)...
                            );

E_p_f_inf = p.n_p_0/p.n_p_0_90d*( exp(-phi_0/p.Tp)*(1-erf(real(sqrt((-phi_m)/p.Tp))))...
            -exp((phi_m-phi_0)/p.Tp)*(1- 2/sqrt(pi)*real(sqrt((-phi_m)/p.Tp))));
        
E_p_c_inf = 2*p.n_p_0/p.n_p_0_90d*( exp((-phi_0)/p.Tp)*(erf(real(sqrt((-phi_m)/p.Tp))))...
            -2/sqrt(pi)*exp((phi_m-phi_0)/p.Tp)*real(sqrt((-phi_m)/p.Tp)));
      
E_i_inf=2*p.n_i_inf*p.Te/(p.Tp*p.n_p_0_90d)*p.Mi^2*(1-real(sqrt(1-2*phi_m/(p.Te*p.Mi^2))));

F(:,1)=E_e_f_inf+E_e_r_inf+E_p_f_inf+E_i_inf;

F(:,2)= p.n_p_0*exp((phi_m-phi_0)/p.Tp)-n_e_inf.*(...
      sqrt(p.Te/p.Tp).*exp(-(real(sqrt(-phi_m/p.Te))-p.V_d/p.Vth_e).^2)+...
      p.V_d/p.Vth_p*sqrt(pi)*erfc(-p.V_d/p.Vth_e+real(sqrt(-phi_m/p.Te))) )...
      +p.n_i_inf*sqrt(2*pi*p.Te/p.Tp*p.me/p.mi)*p.Mi;
end

%% odefcnB

function dydt=odefcnB(t,y,p)

% dydt=zeros(2,1);

dydt(1)=y(2);
dydt(2)=ODEfunB(y(1),p);

end





%% odefcnC

function dydt=odefcnC(t,y,p)

% dydt=zeros(2,1);

dydt(1)=y(2);
dydt(2)=ODEfunC(y(1),p);

end

%% ODEfun

function d2=ODEfun(phi,p)
n_e_f = p.n_e_inf*0.5*exp((phi-p.phi_inf)/(p.Te)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_e_r = p.n_e_inf*exp((phi-p.phi_inf)/(p.Te)).*...
    ( erf(p.V_d/p.V_t) + erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_p_f = p.n_p_0*0.5*exp((phi-p.phi_0)/(p.Tp)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) );

n_i = p.n_i_inf*(1-2*p.e_*phi/(p.V_d^2*p.mi)).^(-0.5);% seems error here, V_d or v_i_inf ?

n_p_c=p.n_p_0*exp((phi-p.phi_0)/p.Tp).*erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) ;

d2 = -p.e_/p.epsilong_0*(n_i-n_e_f-n_e_r-n_p_f)*1e6;

end

%% ODEfunB


function d2=ODEfunB(phi,p)

n_e_f = p.n_e_inf*0.5*exp((phi-p.phi_inf)/(p.Te)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_e_r = p.n_e_inf*exp((phi-p.phi_inf)/(p.Te)).*...
    (erf(real(sqrt((phi-p.phi_m)/(p.Te)))-p.V_d/p.V_t)-erf(-p.V_d/p.V_t));

n_p_f = p.n_p_0*0.5*exp((phi-p.phi_0)/(p.Tp)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) );

n_i = p.n_i_inf*(1-2*p.e_*phi/(p.V_d^2*p.mi)).^(-0.5);% seems error here, V_d or v_i_inf ?

n_p_c=p.n_p_0*exp((phi-p.phi_0)/p.Tp).*erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) ;

d2 = -p.e_/p.epsilong_0*(n_i-n_e_f-n_p_f-n_p_c)*1e6;

end



%% ODEfunC

function d2=ODEfunC(phi,p)

n_e_f = p.n_e_inf*0.5*exp((phi-p.phi_inf)/(p.Te)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_e_r = p.n_e_inf*exp((phi-p.phi_inf)/(p.Te)).*...
    (erf(real(sqrt((phi-p.phi_m)/(p.Te)))-p.V_d/p.V_t)-erf(-p.V_d/p.V_t));

n_p_f = p.n_p_0*0.5*exp((phi-p.phi_0)/(p.Tp)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) );

n_i = p.n_i_inf*(1-2*p.e_*phi/(p.V_d^2*p.mi)).^(-0.5);

d2 = -p.e_/p.epsilong_0*(n_i-n_e_f-n_e_r-n_p_f)*1e6;

end


%% SolveZ

function r=SolveZ(p,n_e_inf,phi_0,phi_m)

a=@(y)1./sqrt(...
            n_e_inf./p.n_p_0_90d.*p.Te/p.Tp.*(exp(y/p.Te).*( 1 - erf(real(sqrt((y-phi_m)/p.Te))-p.V_d/p.Vth_e))...
          - exp((phi_m-p.phi_inf)/p.Te).*(1-erf(-p.V_d/p.Vth_e))+...
          1/sqrt(pi).*p.Vth_e/p.V_d.*exp((phi_m)/p.Te-p.V_d^2/p.Vth_e^2)...
          .*(exp(2.*p.V_d/p.Vth_e.*real(sqrt((y-phi_m)/p.Te)))-1))+...%E_e_f_inf
                            2.*n_e_inf/p.n_p_0_90d.*p.Te/p.Tp*(exp(y/p.Te).*( erf(real(sqrt((y-phi_m)/p.Te))-p.V_d/p.Vth_e)+erf(p.Vd/p.Vth_e))...
          -1/sqrt(pi).*p.Vth_e/p.V_d.*exp((phi_m)/p.Te-p.V_d^2/p.Vth_e^2).*(exp(2*p.V_d/p.Vth_e.*real(sqrt((y-phi_m)/p.Te)))-1))+...%E_e_r_inf
                            2.*p.n_i_inf.*p.Te/(p.Tp.*p.n_p_0_90d).*p.Mi^2*(real(sqrt(1-2.*y/(p.Te.*p.Mi^2)))-real(sqrt(1-2.*phi_m/(p.Te*p.Mi^2))))+...%E_i_inf
                            p.n_p_0/p.n_p_0_90d.*( exp((y-phi_0)/p.Tp).*(1-erf(real(sqrt((y-phi_m)/p.Tp))))...
            -exp((phi_m-phi_0)/p.Tp).*(1- 2/sqrt(pi)*real(sqrt((y-phi_m)/p.Tp))))+...%E_p_f_inf
            2.*p.n_p_0/p.n_p_0_90d.*( exp((y-phi_0)/p.Tp).*(erf(real(sqrt((y-phi_m)/p.Tp))))...
            -2/sqrt(pi).*exp((phi_m-phi_0)/p.Tp).*real(sqrt((y-phi_m)/p.Tp)))...%E_p_c_inf
        );

r=integral(a,phi_m,phi_0);

end


%% STmatrix

function M=STmatrix(z1,z2,z3)
% this function fetch one element of each column of z to form a row, and
% then return a matrix M. If z is a 2 by 3 matrix, then M is a 2^3 by 3
% matrix.
if nargin==3
n1=length(z1);
n2=length(z2);
n3=length(z3);
Nrow=n1*n2*n3;
q=1;
M=zeros(Nrow,3);
for i=1:n1
    for j=1:n2
        for k=1:n3
            M(q,:)=[z1(i),z2(j),z3(k)];
            q=q+1;
        end
    end
end
end
if nargin==2
    n1=length(z1);
    n2=length(z2);
    Nrow=n1*n2;
    q=1;
    M=zeros(Nrow,2);
    for j=1:n1
        for k=1:n2
            M(q,:)=[z1(j),z2(k)];
            q=q+1;
        end
    end
end

end
%% typeA_in2




function p=typeA_in2(Tp,Te,SEA)

%% some constants
me=9.109e-31; % kg
eV2J = 1.6e-19; % eV to J
e_=1.6e-19; % C
k=1.38e-23;% 
epsilong_0=8.85e-12;
%% input variables
% Te = 12; % eV
% Tp = 2.2; % eV
% Te=12*1.6e-19; % J
% Tp=2*1.6e-19; % J
Vd_sw = 468e3; % m/s, solar wind average drifting velocity
% SEA = 90; % sun elevation angle, degree
n_p_0_90d = 64; % at 90d incidance angle, 1/cm3, used as density reference
n_p_0 = n_p_0_90d*sind(SEA); % 1/cm3
n_i_inf=8.7; % 1/cm3

%%
Vd=Vd_sw*sind(SEA);
V_d=Vd;
Vt_p=sqrt(Tp*eV2J/me); % thermal velocity of photoelectron
V_ref = Vt_p;

Vth_e=sqrt(2*Te*eV2J/me); % thermal velocity of solar wind electrons
                          % same as 'Vt_e' as in Jeong p.16, Eq. 2.17
V_t=Vth_e;
Vth_p=sqrt(2*Tp*eV2J/me);

T_hat_e=Te/Tp; % use Tp as normalization reference
beta=Te/Tp; % Jeong dissertation p.18, Eq. 2.28

% normalize
phi_hat_inf=0; % use infinity position as potental reference
phi_inf=0;
n_hat_p_0=1*sind(SEA); % use photoelectron density as normalization reference

V_hat_d=Vd/V_ref;
V_hat_t=Vth_e/V_ref;

n_hat_i_inf = n_i_inf / n_p_0_90d;

mi=1836*me; % proton
Cs=sqrt(Te*eV2J/mi);
Mi=Vd/Cs;

V_hat_inf = Vd / V_ref;
m_hat_i=1836;
x_max=137.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These three parameters should be updated after they are determined with
% other uer-defined functions.

phi_hat_0=0;
phi_hat_m=0;
n_e_inf=0;
phi_m=0;
phi_0=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=struct('Te',Te,'Tp',Tp,'beta',beta,'phi_hat_inf',phi_hat_inf,'Vd',Vd,...
'me',me,'Vth_e',Vth_e,'n_hat_p_0',n_hat_p_0,'n_p_0',n_p_0,'n_p_0_90d',n_p_0_90d,'n_i_inf',n_i_inf,...
'n_hat_i_inf',n_hat_i_inf,'mi',mi,'Cs',Cs,'Mi',Mi,'phi_hat_0',phi_hat_0,...
'phi_hat_m',phi_hat_m,'n_e_inf',n_e_inf,'V_hat_inf',V_hat_inf,'Vt_p',Vt_p,...
'V_hat_d',V_hat_d,'m_hat_i',m_hat_i,'V_hat_t',V_hat_t,'e_',e_,'x_max',x_max,...
'T_hat_e',T_hat_e,'phi_inf',phi_inf,'k',k,'phi_0',phi_0,'phi_m',phi_m,'V_d',V_d,...
'V_t',V_t,'epsilong_0',epsilong_0,'Vth_p',Vth_p);

end


%% Updatep

function p=Updatep(phi,p)

phi_0=phi(1);
p.phi_0=phi_0;
phi_m=phi(2);
p.phi_m=phi_m;
a=1-2*erf(-p.V_d/p.Vth_e)+erf(real(sqrt(-phi_m/p.Te))-p.V_d/p.Vth_e);
b=p.n_p_0*0.5*exp(-phi_0/p.Tp)*(1-erf(real(sqrt(-phi_m/p.Tp))))-p.n_i_inf;
n_e_inf=-b/a*2;
p.n_e_inf=n_e_inf;
end

end