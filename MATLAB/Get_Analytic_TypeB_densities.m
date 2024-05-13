function [z_analytic,swi_analytic,swe_analytic,phe_analytic,tot_analytic]=Get_Analytic_TypeB_densities(Zmin,Zmax,angle)
%% Get the data of z and density



%% Input



Tph=2.2; % eV
Te=12; % eV
lambda=1.378; % m
nl=200; % number of the sample points along z
n0=64; % normalize


%% Analytic

p=typeb_in(Tph,Te,angle);

% find the surface potential, as well as n_e_inf
[p.n_e_inf,p.phi_hat_0]=BonCon(p);
% p.n_e_inf
p.phi_hat_0
% p.phi_hat_0 = -p.phi_hat_0;

p.phi_hat_m=p.phi_hat_0; % Type C photoelectron sheath
p.phi_0=p.phi_hat_0*p.Tp;
p.phi_m=p.phi_inf;

%% now solve the ODE to get the sheath profile

solinit=bvpinit(linspace(Zmin*lambda,Zmax*lambda,nl),[1.0,-0.01]);
sol=bvp4c(@(t,y)odefcn(t,y,p),@(ya,yb)bcfcn(ya,yb,p),solinit);
%% ODE solved
x=linspace(Zmin*lambda,Zmax*lambda,nl);
phi_temp=deval(sol,x);
phi = phi_temp(1,:);

%% Return


z_analytic=x/lambda;

%%




%% electrons
ne=p.n_e_inf/2.*exp((phi-p.phi_inf)/p.Te).*(1-erf(real(sqrt((phi-p.phi_m)/p.Te))-p.V_d/p.Vth_e));
swe_analytic=-ne/n0;

%% ions
ni=p.n_i_inf./sqrt(1-2.*phi/(p.Te.*p.Mi^2));
swi_analytic=ni/n0;


%% phe

nphe=p.n_p_0/2.*exp((phi-p.phi_0)/p.Tp).*(1-erf(real(sqrt((phi-p.phi_m)/p.Tp))));
npheca=p.n_p_0.*exp((phi-p.phi_0)/p.Tp).*erf(real(sqrt((phi-p.phi_m)/p.Tp)));
np=nphe+npheca;
phe_analytic=-np/n0;

%% density

nn=ni-ne-np;
tot_analytic=nn/n0;


% testing
% plot(swe_analytic,z_analytic)



%% bcfcn
function res=bcfcn(ya,yb,p)

res=[ya(1)-(p.phi_0); yb(1)-0];

end


%% bisection


function  m = bisection(f, low, high, tol)
% low
% high
% Evaluate both ends of the interval
y1 = feval(f, low);
y2 = feval(f, high);
i = 0; 
% Display error and finish if signs are not different
if y1 * y2 > 0
error ('Initial interval not correct!')
end 

% Work with the limits modifying them until you find
% a function close enough to zero.
while (abs(high - low) >= tol)
    i = i + 1;
    % Find a new value to be tested as a root
    m = (high + low)/2;
    y3 = feval(f, m);
    
    if y3 == 0
        
        fprintf('Root at x = %f \n\n', m);
        return
    end
%     fprintf('%2i \t %f \t %f \t %f \n', i-1, low, high, m);   

    % Update the limits
    if y1 * y3 > 0
        low = m;
        y1 = y3;
    else
        high = m;
    end
end 
end

%% BonCon

function [n_e_inf,phi_hat_0]=BonCon(p)

phi_hat_0=MyRootFinder(@(x)Confun(x,p),0,100,1e-8,1e3);
a=1-erf(-p.Vd/p.Vth_e);
b=p.n_p_0/2*exp(-phi_hat_0)-p.n_i_inf;

n_e_inf = -b /a * 2;

end

%% Confun

function result=Confun(phi_hat_0,p)

a = 1 - erf(-p.Vd/p.Vth_e);

b = p.n_p_0/2 * exp(-phi_hat_0) - p.n_i_inf;

n_e_inf= - b./a*2;

c = ...
    (   sqrt(p.beta) * exp(-(-p.Vd/p.Vth_e).^2)...
        + p.Vd/p.Vth_p*sqrt(pi)*erfc(-p.Vd/p.Vth_e)...
    );
d = p.n_i_inf * sqrt(2*pi*p.beta*p.me/p.mi)*p.Mi;

result = p.n_p_0*exp(-phi_hat_0) - c.*n_e_inf + d;

end

%% MyRootFinder

function r=MyRootFinder(f,low,high,tol,num)
%This rootfinder is especially definded for type c.
% f is the function handle, [low,high] is the initial interval
% tol is the error tol
% num is the number that [low,high] will be discretized into.
% first, we evaluate f at linspace (low,high,num) and roughly locate a
% valid initial interval for bisection. If more than one intervals are
% valid, then the leftmost is selected. Since it is hard to determine a
% valid initial interval for bisection, this function is going to do this
% hard work.
%--Xinpeng Wei, 07/16/2019

xnodes=linspace(low,high,num);
F=feval(f,xnodes);

if F(1)>0 %From the left to the right, find the first negative value.
   indexR=find(F<0,1,'first');
   indexL=indexR-1;
end
if F(1)==0
    r=low;
    return
end
if F(1)<0 %From the left to the right, find the first positive value.
   indexR=find(F>0,1,'first');
   indexL=indexR-1;
end
if sum(F>0)==num || sum(F<0)==num
    error('No proper initial interval for bisection! Consider increasing num')
end
disp('The root has been bracked into a valid interval')
low=xnodes(indexL);
high=xnodes(indexR);
r=bisection(f,low,high,tol);
end

%% odefcn

function dydt=odefcn(t,y,p)


dydt(1)=y(2);
dydt(2)=ODEfun(y(1),p);

end



%% ODEfun


    function d2=ODEfun(phi,p)
n_e_f = p.n_e_inf*0.5*exp((phi-p.phi_inf)/(p.Te)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_e_r = p.n_e_inf*exp((phi-p.phi_inf)/(p.Te)).*...
    ( erf(p.V_d/p.V_t) + erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_p_f = p.n_p_0*0.5*exp((phi-p.phi_0)/(p.Tp)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) );

n_i=p.n_i_inf*(1-2*p.e_*phi/(p.V_d^2*p.mi)).^(-0.5);


n_p_c=p.n_p_0*exp((phi-p.phi_0)/p.Tp).*erf(real(sqrt((phi-p.phi_m)/(p.Tp))));

d2=-p.e_/p.epsilong_0*(n_i-n_e_f-n_p_c-n_p_f)*1e6;


end

%% Pdensity

function [n_i,n_e_f,n_e_r,n_p_f]=Pdensity(phi,p)

    n_e_f = p.n_e_inf*0.5*exp((phi-p.phi_inf)/(p.Te)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_e_r = p.n_e_inf*exp((phi-p.phi_inf)/(p.Te)).*...
    ( erf(p.V_d/p.V_t) + erf(real(sqrt((phi-p.phi_m)/(p.Te))) - p.V_d/p.V_t) );

n_p_f = p.n_p_0*0.5*exp((phi-p.phi_0)/(p.Tp)).*...
    ( 1 - erf(real(sqrt((phi-p.phi_m)/(p.Tp)))) );

n_i=p.n_i_inf*(1-2*p.e_*phi/(p.V_d^2*p.mi)).^(-0.5);


end



%% typeb_in

function p=typeb_in(Tp,Te,SEA)

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

end