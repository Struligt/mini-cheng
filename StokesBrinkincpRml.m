function StokesBrinkincpRml

% 0411 added pressure ca line 164
% 041: made R fixed, a varying so that we can plot for various k

%modelling of sphere translating with U in Stokes medium in spherical cavity in a larger
%Brinkman medium

close all
clear all
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);

%global parameters
eta=1000;%fluid viscosity

nlines =50; %nlines for contour plot of streamfxn
region = 2.0; %multiplication of R to graph for plot of streamfxn
numplot = 30;
U = 1; %only affects streamfxn calc

%default params (later overridden)
R = 1; %sphere radius                   !! shud not be duplicated w fxn below!!
abyRdef = 0.8
kRdef = 1/(1-abyRdef) %default R / brinkman screening length
fscale = ones(1,numplot);

%define plotting region for drag force
abyR=logspace(-2, 0, numplot); %sphere radius/cavity radius
kR=logspace(-2,2,numplot); %R / brinkman screening length
kRmini=logspace(-1,1,3)

%=========================================================================%
%plot of Stokes scaled drag ag abyR
for j=1:length(kRmini)
    kRuse = kRmini(j);
    for n=1:length(abyR)
        abyRn=abyR(n);
        paramv = paramgen(abyRn,kRuse,eta);
        A=paramv(1,1); %B=paramv(2); C=paramv(3); E=paramv(4); beta=paramv(5); lam = paramv(6); % !!check: lam = 0??
        a = R*abyRn;
        fscale(n)=(2/3)*((eta*a)^-1).*A;
    end %for n
figure(2);
loglog(abyR,fscale,'-r')
hold on
%also plot Brinkman drag (see problem #2, p105) & Stokes Drag lines
end %for j
fStokes = ones(1,length(abyR));
fBrink = ((kRuse.*abyR).^2)/9+kRuse.*abyR+1;
loglog(abyR,fBrink,'--b',abyR,fStokes,':b');
axis([min(abyR) max(abyR)   min(fscale)*0.5  max(fscale)*2])
legend('drag','Stokes', 'Brinkman (k max)')
%grid on

Title('Sphere Drag in Brinkman Stokes Cavity, various k','FontSize',16,'Interpreter','latex');
xlabel('$a/R$','FontSize',16,'Interpreter','latex');
ylabel('$f_d / (6 \pi \eta U a)$','FontSize',16,'Interpreter','latex');
set(gca,'FontSize',16,'FontName','Times')
print('-depsc','f vs abyR Brinkcavity.eps')


%=========================================================================%
%plot of Stokes scaled drag ag k*R
for n=1:length(kR)
    dumkR=kR(n)
    paramv = paramgen(abyRdef,kR(n),eta);
    A=paramv(1); %B=paramv(2); C=paramv(3); E=paramv(4); beta=paramv(5); lam = paramv(6); % !!check: lam = 0??
    fscale(n)=(2/(3*eta*a))*A;
end %for
figure(3);
loglog(kR,fscale,'ok')
hold on
Title('Sphere Drag in Brinkman Stokes Cavity (a/R=0.1)','FontSize',16,'Interpreter','latex');
xlabel('$ R / l_B \equiv k R $','FontSize',16,'Interpreter','latex');
ylabel('$f_d / (6 \pi \eta U a)$','FontSize',16,'Interpreter','latex');
set(gca,'FontSize',16,'FontName','Times')
print('-depsc','f vs kR Brinkcavity.eps')

%=========================================================================%
%calculate streamfxns
R=1; %again
abyRdef = 0.8;
a=R*abyRdef;
kRdef = R/(R-a) %default R / brinkman screening length
k=kRdef/R;
% kRdef = 10^2 %default R / brinkman screening length
% k=kRdef/R;

paramv = paramgen(abyRdef,kRdef,eta); %generate and assign params
A=paramv(1); B=paramv(2); C=paramv(3); E=paramv(4); beta=paramv(5); lam = paramv(6); % !!check: lam = 0??

%draw streamfunctions
% kind thanks to HW3 files from where some of this code taken -
% data generation
x=linspace(-region*R,region*R,401); % $$ inserted 
y=x;
[xx,yy]=meshgrid(x,y);
[tt,rr]=cart2pol(xx,yy);
%generate streamfxn values for all points on plane
for m=1:size(x,2)
    for n=1:size(x,2)
        if(rr(m,n)<a)
            zz(m,n)=0;
        else                
            if (rr(m,n)<R) %ie if less than cavity size, compute streamfunction for inner Stokes region
                r=rr(m,n);
                theta = tt(m,n);
                psi = U*(sin(theta)^2)*(0.5/eta*(A*r+0.1*B*r^4)+0.5*C*r^2-E*r^-1);
                zz(m,n) = psi;            
            else %if in Brinkman region, compute streamfunction for outer Brinkman region
                r=rr(m,n);
                theta = tt(m,n);
                fr = beta*h1f(k,r)*(r^-1) - beta*h1rf(k,r) + 3* lam*r^-3;
                gr=  beta*h1f(k,r)*(r^-1) + beta*h1rf(k,r) - lam*r^-3;
                psi = 0.5*U*(r^2)*(fr+gr)*(sin(theta)^2);
                zz(m,n) = psi;
%                zz(m,n)=0;
            end %if
        end %if     
    end %for
end %for
%plot streamfxn contourlines
figure(4);
contour(xx,yy,zz,nlines,'k');
hold on
circle([0,0],a,50,'k-')
circle([0,0],R,100,'k:')
Title('Sphere Moving in Brinkman Stokes Cavity (a/R=0.8; k(R-a)=1)','FontSize',14,'Interpreter','latex');
xlabel('$x/a$','FontSize',16,'Interpreter','latex');
ylabel('$y/a$','FontSize',16,'Interpreter','latex');
set(gca,'FontSize',16,'FontName','Times')
axis square
axis tight
axis(region*a*[-R,R,-R,R]) % $$ inserted 
print('-depsc','strmlinesBrinkcavity1.eps')



%=========================================================================
%this function generates the needed parameters for calculation of fluid
%flow, drage force etc

function params = paramgen(abyRf,kRf,eta)
R = 1; %sphere radius
a=R*abyRf; %cavity radius
k=kRf/R; %reciprocal Brinkman screening length ie k = 1/lb
U = 1; %sphere velocity

%Hankel functions evaluated at r=R
h = 1i*exp(-k*R)*((k*R)^-1+(k*R)^-2);%first Hankel function
hr = -1i*k*exp(-k*R)*((k*R)^-1 + 2*(k*R)^-2 + 2*(k*R)^-3);%derivative of first Hankel function
hrr = 1i*(k^2)*exp(-k*R)*((k*R)^-1 + 3*(k*R)^-2 + 6*(k*R)^-3 + 6*(k*R)^-4);%2nd derivative of first Hankel function

%solve for params using BCs 
BCa =   [0.5/(eta*a)    -(a^2)/(10*eta)     0   -3*(a^-3)    0                            0]; %no slip along sphere along position vector
BCb =   [0.5/(eta*a)    2*(a^2)/(10*eta)    1   a^-3         0                            0]; %no slip along sphere along flow
BCc =   [0.5/(eta*R)    -(R^2)/(10*eta)     0   -3*(R^-3)    -(h*(R^-1)-hr)       -3*(R^-3)]; %equal velocities at interface, here in direction of position vector
BCd =   [0.5/(eta*R)    0.2*(R^2)/eta       1   R^-3         -(h*(R^-1)+hr)            R^-3]; %equal velocities at interface, here in direction of flow
% BCe =   [-2/(eta*(R^2)) R/(10*eta)        0   18*(R^-4)    (4*h*(R^-2)-4*hr*(R^-1)+hrr)      18*(R^-4)]; %equal traction in direction of position vector   
% BCf =   [0              3*R/(10*eta)      0   -6*(R^-4)    -hrr                    -6*(R^-4)];
BCe =   [-3/(eta*(R^2)) -9*R/(10*eta)       0   18*(R^-4)    (4*h*(R^-2)-4*hr*(R^-1)+hrr)      (1*(k^2)*(R^-2) + 18*(R^-4))]; %equal traction in direction of position vector   
BCf =   [0              3*R/(10*eta)        0   -6*(R^-4)    -hrr                               -6*(R^-4)];

BC =    [BCa; BCb; BCc; BCd; BCe; BCf];
BCv =   [0 1 0 0 0 0]';
%solve for params
%params = BC \ BCv;
params = inv(BC)*BCv;
return

function h1f = h1f(k,r) %hankel
    h1f = 1i*exp(-k*r)*((k*r)^-1+(k*r)^-2);%first Hankel function 
return

function h1rf = h1rf(k,r) %hankel
    h1rf = -1i*k*exp(-k*r)*((k*r)^-1 + 2*(k*r)^-2 + 2*(k*r)^-3);%derivative of first Hankel function
return
% function h1rrf = h1rrf(k,r) %hankel
%     h1rrf = i*k^2*exp(-k*R)*((k*R)^-1 + 3*(k*R)^-2 + 6*(k*R)^-3 + 6*(k*R)^-4)%2nd derivative of first Hankel function
% return

function H=circle(center,radius,NOP,style)
%------------------------------------------------------
% H=CIRCLE(CENTER,RADIUS,NOP,STYLE)
% This routine draws a circle with center defined as
% a vector CENTER, radius as a scaler RADIS. NOP is 
% the number of points on the circle. As to STYLE,
% use it the same way as you use the rountine PLOT.
% Since the handle of the object is returned, you
% use routine SET to get the best result.
%
% Usage Examples,
%
% circle([1,3],3,1000,':'); 
% circle([2,4],2,1000,'--');
%
% Zhenhai Wang <zhenhai@ieee.org>
% Version 1.00
% December, 2002
%-------------------------------------------------------

if (nargin <3),
error('Please see help for INPUT DATA.');
elseif (nargin==3)
style='b-';
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
H=plot(X,Y,style,'linewidth',2);
axis square;
