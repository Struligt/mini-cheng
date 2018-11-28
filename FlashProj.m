function [Vfracf, xif, Kif] = FlashProj(P,T,zi)

%Inputs: P in psia, T in Celcius
%Outputs:

% this performs a flash calculation for mixture
% Peng-Robinson EOS

%remove clear when calling function
%clear all
%close all

hold off
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);
set(0,'defaultaxeslinewidth',1);
set(0,'defaultlinelinewidth',1);


%***************** Input Variables *********************************

%***************** Debugging Variables *********************************
%P = 100;
%T = -59.44;
% SPECIFY FEED MOLE FRACTION
%zi = [0.62 0.38];% 1]; %Tc level %this algorithm s ability to converge sensitive to this feed composition 

zi=zi/sum(zi);

R= 8.314;
T= 273.15 + T; 
P=6894.757 * P; %converted to Pa

%***************** Mixture Defintion *********************************
%for ethane(1), methane(2), propane (3)
nC= 3; %number of components
%critical parameters to be filled in
%inputs for generalized EOS below:
Tbi=[-88.6 -161.4 -42.2];
Tbi=Tbi+273.15;
Tci=[305.32 190.564 369.83];
Pci=[4.872e6 4.599e6 4.248e6];
wi= [0.0995 0.0115 0.1523]; %acentric factor w=-1.0-log(Pvap @Tr=0.7 /Pc)

Ki =[0.2 2.5 0.1];
Vfrac=0.5; %start first guess here of vapor fraction

flashloops=100;
tolerance = 1e-7; %tolerance for N-R root solving of vapor fraction Vfrac
tol = 1e-7; %tolerance for equilibrium test = sum(log(fivap/filiq))
convtol = 1e-20;
eqmtest=1000; %dummy variable to store equilibrium test = ""

%***************** EOS & EOS Mixing Parameters*********************************

binp=zeros(nC); %binary interaction parameter fudge factor for better fitting of PR in mixture

%PR constants ac(Tc), b(Tc), kappa(w) 
aci=0.45724*(R^2).*(Tci.^2)./Pci;
bi=0.07780*R.*Tci./Pci;
kappai=0.37464+1.54226.*wi-0.26992.*(wi.^2);

alphai=(1+kappai.*(1-sqrt(T./Tci))).^2;
ai=aci.*alphai;%a(T) =ac(Tc)*alpha(T,kappa(w)),
Ai = ai*P/((R*T)^2);
Bi = bi*P/(R*T);

%************* FLASH CALCULATION ****************************************
n=1; %counter of loop
while (n < flashloops) && ((eqmtest*eqmtest) > (tol*tol))

%A=== given Ki and feed mole fractions, determine mixture vapor fracn =====

g0=sum(zi.*(Ki-1));
g1=1-sum(zi./Ki);
if g0*g1<=0
    %then 2phases so solve for Vfrac by N-R
    i=1; %counter of iterations of N-R        
    obj = sum(zi.*(Ki-1)./(1+Vfrac.*(Ki-1))); %objective fxn of N-R root solving
    while i < 50 && (obj*obj)>(tolerance*tolerance) %perform N-R root solving                
        dobjdx = -sum((zi.*(Ki-1).^2)./((1+Vfrac.*(Ki-1)).^2));
        Vfrac = Vfrac - obj / dobjdx;
        obj = sum(zi.*(Ki-1)./(1+Vfrac.*(Ki-1)));
        i=i+1;
    end %while
    np=2;
    %'2 phases'
else
    if g0<=0 %feed is all liquid
        Vfrac=0;
        np=1;
    else
        if g1>=0 %feed is all vapor
            Vfrac=1;
            np=1;
        else %neither VLE, nor L, nor V:
            %print error message
            'Error in # phase calculation'
            break
            %quit program
        end
    end
end %if

xi=zi./(1+Vfrac.*(Ki-1)); 
yi=Ki.*xi;
%normalization of mole fraction
xi=xi./sum(xi); %liquid mole fractions
yi=yi./sum(yi); %vapor mole fractions

%B======== calculation of fugacities===================================== 

%Calculation of mixing params for P-R EOS
bvapmix = sum(yi.*bi); bliqmix = sum(xi.*bi);

avapmix= sum(sum(((yi.*sqrt(ai))'*(yi.*sqrt(ai))).*(1-binp)));
aliqmix= sum(sum(((xi.*sqrt(ai))'*(xi.*sqrt(ai))).*(1-binp)));

Avapmix = avapmix*P/((R*T)^2); 
Aliqmix = aliqmix*P/((R*T)^2);
Bvapmix = bvapmix*P/(R*T); 
Bliqmix = bliqmix*P/(R*T);

%calculate V for vapor
    A=Avapmix;  B=Bvapmix;
    Z = roots([1 (B-1) (A-3*B^2-2*B) (-1*(A*B-B^2-B^3))]);
    % reduce the 3 roots of Z by first removing any imaginary and negative roots and then taking minimum(liquid) and maximum(vapor)
    Z = sortv(Z); Zvap=max(Z);
    Vvap = max(Z)*(R*T/P);
%calculate V for liquid
    A=Aliqmix;  B=Bliqmix;
    Z = roots([1 (B-1) (A-3*B^2-2*B) (-1*(A*B-B^2-B^3))]);
    % reduce the 3 roots of Z by first removing any imaginary and negative roots and then taking minimum(liquid) and maximum(vapor)
    Z = sortv(Z); Zliq=min(Z);
    Vliq = min(Z)*(R*T/P);
 
%now calculate mixture fugacities for each component
fivap=P.*yi.*exp( (bi./bvapmix).*(Zvap-1)-log(Zvap-bvapmix*P/(R*T))-(avapmix/(2*sqrt(2)*bvapmix*R*T))...
             *(2*((sqrt(ai)'*sqrt(ai).*(1-binp))*yi')'/avapmix - bi/bvapmix)...
             *log((Zvap+(1+sqrt(2))*bvapmix*P/(R*T))/(Zvap+(1-sqrt(2))*bvapmix*P/(R*T))));
        
        
filiq=P.*xi.*exp( (bi./bliqmix).*(Zliq-1)-log(Zliq-bliqmix*P/(R*T))-(aliqmix/(2*sqrt(2)*bliqmix*R*T))...
            *(2*((sqrt(ai)'*sqrt(ai).*(1-binp))*xi')'/aliqmix - bi/bliqmix)...
            *log((Zliq+(1+sqrt(2))*bliqmix*P/(R*T))/(Zliq+(1-sqrt(2))*bliqmix*P/(R*T))));
        
%============== testing for equilibria ================================== 
        
%testing phase equilibria
eqmtest0 = eqmtest;
eqmtest = sum(log(fivap./filiq));    
convtest = abs(eqmtest-eqmtest0);

%store values for later plot of convergence and eqm test for debugging
    if n == 1
        nplot=n;convplot=convtest;eqmplot=eqmtest*eqmtest;
    else
        nplot=[nplot n]; convplot=[convplot convtest]; eqmplot=[eqmplot eqmtest*eqmtest];
    end
        
if (eqmtest*eqmtest)<=(tol*tol)
    %'Flash converged to VLE mixture'
    break %quit while loop and script
else    
    if n==1 
        Ki = Ki.*(filiq./fivap);%.^1.5; 
    else
        if convtest < convtol
            if eqmtest >0
                %'Liquid feed'                
                %break
            else
                %'Vapor feed'
                 %break
            end%if    
        else
            Ki=Ki.*(filiq./fivap);%.^1.5;            
        end %if    
    end %if    
end %if
n=n+1;
end %while

%output to function
npf=np;Vfracf=Vfrac;nf=n;Kif=Ki; xif=xi;

%output results to screen
% sprintf('# of phases in feed is %d, feed vapor fraction %d; with following compositions(trials taken to reach this: %d)',np,Vfrac,n)
% zi = zi+0; %sprintf(' For feed with mole ratio of %d', zi)
% xif = xif+0
% yif = yi+0
% Kif = Kif+0
% Vvap = Vvap+0
% Vliq = Vliq+0
% fugratio = filiq./fivap

% %plot convergence
% figure(1)
% clf
% hold off
% [AX,hline1,hline2] = plotyy(nplot,convplot,nplot,eqmplot,'semilogy','plot');
% hold on
% %set axes colors
% set(gca,'Color','w','ycolor','k'); %background color white and y axes black
% set(get(gca,'Ylabel'),'Color','k'); %ylabel to black
% set(gcf,'Color','w');%background color white
% %set axes
% set(hline1,'LineStyle','-.','color','k')
% axes(AX(1))
% ylabel('Covergence of Eqm Metric')
% 
% set(hline2,'LineStyle','--','color','r', 'linewidth',1.0)
% axes(AX(2))
% ylabel('Equilibrium Metric')
% legend('convergence of Eqm Metric (sum ln(fivap/filiq))','Eqm Metric = sum of ln(fivap/filiq)')

