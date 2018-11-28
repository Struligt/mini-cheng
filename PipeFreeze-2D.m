
function PFFEM2D
%transient home pipe freeze problem in 2D
%Note : 
    %1. importance of specifying DBCs for stationary problem at least. If no DBC imposed, then won't converge to solution. 
    %(probably as then trivial solution applies)
    %2. stationary solution exists, but not all impositions of BC s work.
    %For example, Robin BC (convective heat transfer) BC imposition does
    %not give convergence if imposed on y==0, x>0 but does if imposed on
    %x==xmin
    %3. make sure that signs are correct for the line integral around bdry
    %: ie if integral of dx but normal du/dy points down,  then -integral dx
    %of du/dy ; if x direction is negative x, then also -integral etc.
% :: PFFEM2D.m

%Transient Solution

clear all; close all; hold off;
set(0,'defaultaxesfontsize',10); set(0,'defaulttextfontsize',10); set(0,'defaultaxeslinewidth',1); set(0,'defaultlinelinewidth',1);

%Constants
W = [0.2777778 0.444444 0.27777778]; %[5/18 4/9 5/18] 
GP = 0.5*[-0.7746 0 0.7746]+0.5; %Gaussian points over [0,1] 
method = 1; %default - can't change this

%== process params ======================================================
R = 0.75/2.54/2;
ksoil = 0.25;
kconcrete = 0.9;

hstart = 0.5;
hend = hstart+0.01;

cthick = 0.1; %thickness of concrete
pdepth = 1.5; %pipe depth

hhole = 5; %10[ W/m2.K] for non-circulating air
hgnd = ksoil/pdepth;
hin = 0.05; %assume high here as pipe insulated 

L = 5;%about 5m
%xhole = 0.5;
alphaw = 1e-7;
kappaw = 0.6;

%BCs
Tin = 17; %17C
Tmn = 4; %4C 
%IC
Tstart = 5;

%outside temperature
Tav = 10; Tvar = 20; %*later to be amended to be sinusoidal* 
Tinf = Tav; %initial (tbe changed later)
Tmax = Tav+Tvar; Tmin = Tav-Tvar;

tc = L^2/alphaw;
freq = 2*pi()/(365*24*60*60/tc);

%adimensionalize above
Ar = (L/R); %aspect ratio
xhstart = hstart/L;
xhend =  hend/L;

Tfus = (0-Tmin)/(Tmax-Tmin);

htc = 1; %heat diffusivity of soil TBReviewed*******
Bi = 1;
%========================================================================
%=========define time space==============================================
%time parameters
tstart = 0;
tfinal = 0.25;
disp(['tstart= ', num2str(tstart),' ,tend= ',num2str(tfinal*tc/(3600*24*7)),' wks']) 
DT = 0.01;%0.1 %OBS : set this too high and your output may be unintelligible
EPS = .0001;%.001 %LTE tolerance
tsteptol = 0.8;
%========================================================================
DTO=DT;
t=tstart;

%===parametric data======================================================
convtol = 1e-5; %convergence tolerance
maxNR = 30; %maximum number of N-R iterations
v1 = 2;%nonuniform grid spacing if~=1
%========================================================================

%variables tb used later
L1 = 0; M1 =0;%global node counter
err = 1; %Euclid norm of convergence
residerr = 1;
NBC = 0;
x1 = 0;
I1 = 0;
elem = 0; % element : used in plotting part

%explicit marching fwd
%==Define Mesh=====================================================
%define solution space
xmax=+4;
xmin=-xmax;   
ymin = -4; ymax = 2;

NEX1 = 8; %# of elements for top of the "L" in x direction
NEX2 = NEX1*xmax/abs(xmin);
NEX = NEX1+NEX2 %bottom of the "L" (full length)
NEY1 = 16; %the "L" height
NEY = NEY1;
NEY2 = ceil(NEY1*(-ymin/(ymax-ymin))) %rounded up to ensure integer # of elements

NE1 = NEX1*NEY1;
NE2 = NEX2*NEY2;
NE = NE1 + NE2;

NX1 = NEX1+1; NX2 = NEX2; NX = NX1+NX2; maxNX = NX; %# nodes in x direction
NY1 = NEY1+1;NY2 = NEY2+1; maxNY = NY1; %NY = NY1+NY2;% OBS NY2
N1 = NX1*NY1;
N2 = NX2*NY2;
N=N1 + N2; %NX*NY; %# nodes = # of global basis functions

xcoords = (linspace(xmin, xmax, maxNX).^1)
% xcoords = [-(abs(linspace(xmin,0,maxNX).^v1)) linspace(0,xmax,maxNX).^v1]%(linspace(xmin, xmax, maxNX).^1); 
% xcoords = sort(unique(xcoords)) %get rid of duplicate 0
ycoords = (linspace(ymin, ymax, maxNY).^1);
%==================================================================

numlocal = 4; %# of local basis functions ==4 for linear case
N2a = [0 1 NY1 NY1+1]; %used in NOP function: gives stepping of local nodes in sawtooth convention of naming nodes as you go around element 
N2b = [0 1 NY2 NY2+1];
W1 = zeros(1,numlocal);

%a. Compute Mesh
xnodes = zeros(1, N);ynodes = zeros(1, N);
for i =1:N1 %fill in each nodes coordiate
    xnodes(i) = xcoords(1+floor((i-1)/NY1));
    ynodes(i) = ycoords(i-floor((i-1)/NY1)*NY1);
end %for :
ilast = 1+floor((N1-1)/NY1);
for i =1:N2 %fill in each nodes coordiate
    xnodes(N1+i) = xcoords((ilast+1)+floor((i-1)/NY2));
    ynodes(N1+i) = ycoords(i-floor((i-1)/NY2)*NY2);
end %for :
xnodes = xnodes; %debug check
ynodes = ynodes; %OK! 

%node locater array : NOP(i,j) is a NE X 4 (4 local linear basis functions
%per element) giving global nodes for each global element i & local node j
NOP = zeros(NE,4);

elch = NE1+NEY2;
for i =1:elch %fill in first section & first adjoining elements
    I1 = i+floor((i-1)/NEY1);    
    for j=1:numlocal        
        NOP(i,j)= I1+N2a(j);                 
    end% loop over local basis
end% loop over NE
lastI1 = I1;
for i =(elch+1):NE %fill in 2nd section & first adjoining elements
    I1 = i+floor((elch-1)/NEY1)+(NEY1+1-NEY2)+floor(((i-elch)-1)/NEY2);    
    for j=1:numlocal        
        NOP(i,j)= I1+N2b(j);                 
    end% loop over local basis
end% loop over NE
NOP = NOP; %debug check

%b. Plot mesh & nodes
% figure(1);hold on;
% %a. Plot mesh
% for i = 1: NX, plot(xcoords(i)*ones(1,length(ycoords)), ycoords,':r'), end;
% for j = 1: NY1, plot(xcoords, ycoords(j)*ones(1,length(xcoords)),':r'), end;
% %- plot coords of each node
% for i =1: N
% %    text(xnodes(i),ynodes(i),['(',num2str(xnodes(i),'%1g'), ', ',num2str(ynodes(i),'%1g'),')'],'FontSize',10,'Interpreter','latex')
% end;
% %- label each node
% for i =1: NE
%     text(xnodes(NOP(i,1))+0.5,ynodes(NOP(i,1))+0.5,[num2str(i)],'FontSize',12,'Interpreter','latex') %label element #
%     for j =1:numlocal 
%         text(xnodes(NOP(i,j)),ynodes(NOP(i,j)),[num2str(NOP(i,j))],'FontSize',10,'Interpreter','latex') %label nodes
%     end;    
% end; %for i     
% pause

%==Step 2: Define BCs & Initial Guess===================================================
NBCx1min = 0; %at x=xmin, no heat transfer as no gradient
NBCx1max = 0*1;%at x=0, heat transfer
NBCx2max = 0; 
NBCy1max = 0*1;%top of L
NBCy2max = 0*1;%bottom top sectopm of L
NBCymin = 0;
us = ones(N,1)*0.5; %initial nodal values == unknown basisfxn wts tb solved for===========

% << here define initial guess <<<< OBS: defining it HERE makes a big differnce in
% proper functioning of calculation of ux and uy
for i = 1: NX1 %fix top BC
    us(i*NY1) = sin(xnodes(i*NY1)*pi()/(-xmin));
end %for i 
%=========================================================================

disp(['linear basis fxns used, NE=',num2str(NE),' :'])

Bi = 1;%----

% SOLVE GFEM
NRiter=0; %N-R iteration counter
while err > convtol
NRiter = NRiter+1; %N-R iteration counter

%Step 3: Calculate general wted residual vector and Jacobian of Stiffness wrt weights
Jij=zeros(N); %jacobian (==stiffness matrix if linear) tb recalculated for kth iteration 
Fi=zeros(N,1); %wted residual vector (==load matrix if linear) tb recalculated for kth iteration
for i = 1: NE %loop over all elements & calculate wted residual vector Fi and Jacobian Jij for NR converging
    for j=1:4, %store the global nodes relavent to the ith element for later use
        W1(j)=NOP(i,j); 
    end %for j: 
    dx = xnodes(W1(3))-xnodes(W1(1));   %<<<<< changed Fortran scratch code 
    dy = ynodes(W1(2))-ynodes(W1(1));    %<<<<<     ""
    for j = 1:length(GP) %perform Gaussian integral over element i                 
        for k =1:length(GP)
        [phi, phix, phiy] =  blinbasfxn(GP(j),GP(k),dx, dy); % <<

        %calculate function value at each of the GPs 
        u = 0;ux = 0; uy=0; %initialize basis fxn weights and basis fxn weight derivatives
        for L=1:numlocal %calculate FEM's approx function at relevant GP                           
            u = u+phi(L)*us(W1(L)); %FEM's approx function value at GP = sum(ui*phi)
            ux = ux+phix(L)*us(W1(L));
            uy = uy+phiy(L)*us(W1(L));            
        end %loop over local basis functions 
        
        
        for L = 1:numlocal %calculate using 2X1 and 2X2 Jij and Fi                          
            L1 = W1(L); %global node
            %in the following 2 loops the particular GFEM form of the actual GE has to be specified          
            %==============================================================
            Fi(L1) = Fi(L1)+dx*dy*W(j)*W(k)*(Bi*(ux*phix(L)+uy*phiy(L))); %missing Neuman BC terms as to be added in below            !!0*dcon*phi(L)+ 
            for m = 1:numlocal 
                M1 = W1(m); %global node                                              
                Jij(L1, M1) = Jij(L1, M1)+dx*dy*W(j)*W(k)*(Bi*(phix(m)*phix(L)+phiy(m)*phiy(L))); % !!! +0*phi(m)*phi(L)/DT+
            end % for m    
            %=============================================================
        end % L: filling of 2X1 and 2X2 Jij and Fi                             
        end % k: Gaussian integration in x == zeta direction
   end % j : Gaussian integration in y == eta direction   
%disp(['i=',num2str(i),': dx=',num2str(dx),', dy=',num2str(dy),',u=',num2str(u),', ux=',num2str(ux)]) %debug
end %i : loop over all the elements
%Fi = Fi, Jij = Jij

%Step 4: impose Dirichlet and Neumann BCs on J & F
%=======================================================================
% fix Dirichlet BCs into our Galerkin equation
%  
% for i = 1: NX1 %fix top BC
%     us(i*NY1) = sin(xnodes(i*NY1)*pi()/(-xmin));
%     Fi(NY1*i) = 0; 
%     Jij(NY1*i,:)=zeros(1,N);
%     Jij(NY1*i,NY1*i)=1;
% end %for i 

%find other nodes that lie on BC
% DBC = 0;
% BCnodeis = [find(ynodes == ymax)]; %find(ynodes == 0 & xnodes >=0) 
% BCnodeis = sort(unique(BCnodeis));
% for i = 1: length(BCnodeis) %apply BC at these nodes
%     Fi(BCnodeis(i)) = 0; 
%     Jij(BCnodeis(i),:) = zeros(1,N);
%     Jij(BCnodeis(i),BCnodeis(i)) = 1;
%     us(BCnodeis(i))= DBC;
% end

DBC = 0.5;
BCnodeis = [find(ynodes == ymin)]; %find(ynodes == 0 & xnodes >=0) 
BCnodeis = sort(unique(BCnodeis));
for i = 1: length(BCnodeis) %apply BC at these nodes
    Fi(BCnodeis(i)) = 0; 
    Jij(BCnodeis(i),:) = zeros(1,N);
    Jij(BCnodeis(i),BCnodeis(i)) = 1;
    us(BCnodeis(i))= DBC;
end

% DBC = 0;
% BCnodeis = [find(ynodes == ymax) find(xnodes == 0 & ynodes >=0) find(ynodes == 0 & xnodes >=0)]; %find(ynodes == 0 & xnodes >=0) 
% BCnodeis = sort(unique(BCnodeis));
% for i = 1: length(BCnodeis) %apply BC at these nodes
%     Fi(BCnodeis(i)) = 0; 
%     Jij(BCnodeis(i),:) = zeros(1,N);
%     Jij(BCnodeis(i),BCnodeis(i)) = 1;
%     us(BCnodeis(i))= DBC;
% end

% imposition of below Robin BC s doesn't converge : not sure why
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% % add convective heat transfer Robin BCs to our wted residual vector Fi
Tinf = 1;
BCnodeis = [find(xnodes == xmin)];% & xnodes >=0)]; %set of boundary nodes making up SINGLE line 
BCnodeis = sort(unique(BCnodeis)); %apply BC at these nodes
for i = 1: (length(BCnodeis)-1) 
    BCIs = [BCnodeis(i) BCnodeis(i+1)];%node #s along top of L BC
    dy = ynodes(BCIs(2))-ynodes(BCIs(1));
    for j = 1:length(GP) %calculate line integral of NBCs that we include in the Fi residual
        [phi] =  lbasfxn(GP(j));            
        u=0; %initialize the approx fxn
        for L=1:2
            u = u+phi(L)*us(BCIs(L)); %FEM's approx function value at GP = sum(ui*phi)
        end %for         
%         x1 = xnodes(BCIs(1))+dy*GP(j); %x coordinate of GP(j)
%         Tinf = sin(x1*pi()/xmax);
        RBC = -Bi*(u-Tinf); %convective heat transfer        
        for L =1:2
           Fi(BCIs(L))=Fi(BCIs(L))-dy*W(j)*phi(L)*(RBC);
        end                
    end % j Gaussian integral
end %for i
% 
% % add convective heat transfer Robin BCs to our wted residual vector Fi
% Tinf = 0;
% BCnodeis = [find(xnodes == 0 & ynodes >=0)]; %***set of boundary nodes making up SINGLE line 
% BCnodeis = sort(unique(BCnodeis))
% %apply BC at these nodes
% LI = 0; %line integral to be calculated : initialize
% for i = 1: (length(BCnodeis)-1) 
%     BCIs = [BCnodeis(i) BCnodeis(i+1)];%node #s along top of L BC
%     dy = ynodes(BCIs(2))-ynodes(BCIs(1));
%     for j = 1:length(GP) %calculate line integral of NBCs that we include in the Fi residual
%         [phi] =  lbasfxn(GP(j));            
%         u=0; %initialize the approx fxn
%         for L=1:2
%             u = u+phi(L)*us(BCIs(L)); %FEM's approx function value at GP = sum(ui*phi)
%         end %for         
%         y1 = ynodes(BCIs(1))+dy*GP(j); %x coordinate of GP(j)
%         %Tinf = 0;%sin(x1*pi()/xmax);
%         RBC = -Bi*(u-Tinf); %convective heat transfer        
%         for L =1:2
%            Fi(BCIs(L))=Fi(BCIs(L))+dy*W(j)*phi(L)*(RBC);
%         end                
%     end % j Gaussian integral
% end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%=========================================================================


%Step 5: invert to get next iteration's unknown basis nodal values == basis fxn weights 
A = inv(Jij)*(-Fi);
us = us + A; %output weights
%Step 6: calculate error == Euclid norm of diff of N-R iterations of C
err = sqrt(A'*A);
if NRiter> maxNR, break, end
disp(['N-R iterations : ', num2str(NRiter),' ,convergence eror norm of : ',num2str(err)])%' , wted resid vector norm of ',num2str(err)])
end % while

disp(['Final N-R iterations : ', num2str(NRiter),' ,convergence eror norm of : ',num2str(err)])%' , wted resid vector norm of ',num2str(err)])
if NRiter>maxNR, disp(['OBS: No N-R convergence in ', num2str(maxNR),' iterations']), end
% 

%Contour map
%simple plot of basis wt fxns only
figure(1);
meshnum = NEY1;
ti = xmin:((xmax-xmin)/meshnum):xmax; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(xnodes,ynodes,us,XI,YI);
contour(XI,YI,ZI,50);
hold on
plot3(xnodes,ynodes,us,'.k','MarkerSize',4)


% %Step 6: Plot Solution us for different y's
% % simple plot of basis wt fxns only
% % usplot = zeros(NY,NX); %us to plot
% % for j =1:NY
% %     for i =1: NX
% %         usplot(j,i) = us(j+(i-1)*NY);
% %     end %i        
% % end %j 
% % figure(2);
% % plot(xcoords,usplot,'-')        
% 
% %better plot of GFEM solution for any specified x and y
% %=========================================================================
% x=linspace(xmin,xmax,100); %specify here range over which to plot
% y = linspace(ymin, ymax, 6);
% %=========================================================================
% ua = zeros(length(y),length(x)); %GFEM approximate solution tb filled in 
% for m =1:length(y)
%     iy = max(find(ycoords < y(m))); %bottom node for each y in specified vector to plot 
%     if isempty(iy), iy = 1; end;
%     dy = ycoords(iy+1)-ycoords(iy);
%     etay = (y(m) - ycoords(iy))/dy;
%     for n=1:length(x)    
%         ix = max(find(xcoords < x(n))); %LHS node for each x in specified vector to plot         
%         if isempty(ix), ix = 1; end;
%         dx = xcoords(ix+1)-xcoords(ix);
%         etax = (x(n) - xcoords(ix))/dx;
%         [phi, phix, phiy] =  blinbasfxn(etax,etay,dx, dy);    
%         ua(m,n) = 0; %approximate value of function tb calculated for at (y,x) (OK confusing - indexing to be changed)
%         elem = (ix-1)*NEY+iy; %element number  
%         for L=1:numlocal %calculate FEM's approx function at relevant GP              
%             ua(m,n) = ua(m,n)+phi(L)*us(NOP(elem,L)); 
%         end %loop over local basis functions    
%     end % loop over all points specified to be plotted
% end %for m 
% figure(1); 
% plot(x,ua);
% legend(num2str(y'),'Location','Best');
% hold on; 
% text(x(5),ua(1,5),['y = ', num2str(y(1),'%g')],'FontSize',12,'Interpreter','latex');
% text(x(5),ua(1,length(x-10)),['NR iterations=',num2str(NRiter),', NR-error= ' num2str(err,'%.1g')],'FontSize',12,'Interpreter','latex')
% title(['GFEM solution of 2D NL PDE, used order ',num2str(method),', NEX=',num2str(NEX),', NEY=',num2str(NEY)],'FontSize',10,'Interpreter','latex');
% xlabel('x','FontSize',12,'Interpreter','latex'); 
% ylabel('Numerical Solution Ua(x,y)','FontSize',12,'Interpreter','latex');
% grid on; 
% print('-depsc',['CHEE662HW2a_',num2str(method),'o_',num2str(NEX),'x',num2str(NEY),'mesh.eps'])

end %function GFE2b
%==========================================================================
%Local Basis Functions:
%A. 1D linear piecewise Lagrangian polynomials
function [phi] = lbasfxn(GP) %returns value of LINEAR local basis fxn phi L = 1,2 
    phi(1) = 1- GP;
    phi(2) = GP;    
end %function phi

function [phi, phix, phiy] =  blinbasfxn(GP1,GP2,dx, dy) %returns value of BILINEAR local basis fxn phi L 
    phi(1) = (1- GP1)*(1-GP2);
    phi(2) = (1-GP1)*GP2;
    phi(3) = GP1*(1-GP2);
    phi(4) = GP1*GP2;
    
    phix(1) = -(1-GP2)/dx;
    phix(2) = -GP2/dx;
    phix(3) = (1-GP2)/dx;
    phix(4) = GP2/dx;

    phiy(1) = -(1-GP1)/dy;
    phiy(2) = (1-GP1)/dy;
    phiy(3) = -GP1/dy;
    phiy(4) = GP1/dy;
end %function bilinear basis function
