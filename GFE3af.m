function GFE3af

% sample 2D GFE solution to linear 2D PDE 
%solves linear 2D PDE Del2(u)+1=0 with D (no slip) BCs

% GFE3ae but now added calculation of rel error of nodes

clear all; close all; hold off;

%Constants
W = [0.2777778 0.444444 0.27777778]; %[5/18 4/9 5/18] 
GP = 0.5*[-0.7746 0 0.7746]+0.5; %Gaussian points over [0,1] 
method = 1; %default - can't change this

%===parametric data======================================================
roundtol = 1e12; %1/tolerance for which Matlab rounds to nearest signficant integer
convtol = 1e-6; %convergence tolerance
maxNR = 30; %maximum number of N-R iterations
Ri = 1; %dimensionless Ringskog number = (viscosity * max velocity) / (- pressure gradient * lc^2)
meshnum = 20; 
relerrsQ = [];
relerrs = [];%relative errors of nodal points
relerr2s = [];%approximated relative errors
%========================================================================


%==Define Mesh=====================================================
%define solution space
xmin=0; xmax=1;%xmax=2;  
ymin = 0; ymax = xmax; %xmax and ymax have to be equal in this mesh
Rmin = 0.5; Rmax = 1; %adimensionalized radius
%Rmin = 1; Rmax = 2; %this is really only maximum radius along theta =0 or theta =pi/2
Tmin = 0; Tmax = pi/2;

NEgs = [4 8 16 32 64]; %OBS : (128+1)^2 array does not work: not enough memory
for iNEg = 1:length(NEgs) 

%variables tb used later
L1 = 0; M1 =0;%global node counter
err = 1; %Euclid norm of convergence
x1 = 0;
I1 = 0;
detJt = 0; %determinant of transformation Jacobian
Qin = 0; Q = 0; %total flow rate
    
NEg= NEgs(iNEg); %global variable # of elements
NER = NEg; %# of elements along radial direction
NET = NEg; %# of elements along theta(angular) direction btw 0 and pi/2
NE = NER*NET;
NR = NER+1; %# nodes in radial direction
NT=NET+1;
N=NT*NR; %# nodes = # of global basis functions

Tcoords = linspace(Tmin, Tmax, NT); 
%fill in R coords
Rcoords = zeros(NT, NR);
for i=1:NT
    if Tcoords(i) <= pi/4
        Rcoords(i,:)= linspace(Rmin, (ymax-ymin)/cos(Tcoords(i)), NR);
    else
        Rcoords(i,:)= linspace(Rmin, (xmax-xmin)/sin(Tcoords(i)), NR);
        %Rcoords(i,:)= linspace(Rmin, (xmax-xmin)/cos(pi/2-Tcoords(i)), NR);
    end % if
end 

numlocal = 4; %# of local basis functions ==4 for linear case
N2 = [0 1 NR NR+1]; %used in NOP function: gives stepping of local nodes in sawtooth convention of naming nodes as you go around element 
W1 = zeros(1,numlocal);
phix = zeros(1, numlocal); phiy = zeros(1, numlocal);
%==================================================================

%a. Compute Mesh
Tnodes = zeros(1, N); Rnodes = zeros(1, N);
xnodes = zeros(1, N); ynodes = zeros(1, N);
for i =1:N %fill in each node's coordinates     
    Tnodes(i) = Tcoords(1+floor((i-1)/NR));
    Rnodes(i) = Rcoords(1+floor((i-1)/NR),i-floor((i-1)/NR)*NR);        
    xnodes(i) = xmax - Rnodes(i)*sin(Tnodes(i)); 
    ynodes(i) = ymax - Rnodes(i)*cos(Tnodes(i)); 
end %for :

%!!!!!necessary for later BC checking : this is a bug with computer memory and
%Matlab !!!!!!! 
xnodes = round(xnodes*roundtol)/roundtol;
ynodes = round(ynodes*roundtol)/roundtol;
Rnodes = round(Rnodes*roundtol)/roundtol;

%OK
%node locater array : NOP(i,j) is a NE X 4 (4 local linear basis functions
%per element) giving global nodes for each global element i & local node j
NOP = zeros(NE,4);
for i =1:NE    
    I1 = i+floor((i-1)/NER);
    for j=1:numlocal        
        NOP(i,j)= I1+N2(j);                 
    end% loop over local basis
end% loop over NE

%b. Plot mesh & nodes
if iNEg == 1
    figure(1);hold on;
    plot(xnodes, ynodes,'or');
    title(['Mesh used with NER=',num2str(NER),', NET=',num2str(NET)],'FontSize',12,'Interpreter','latex');
    xlabel('x','FontSize',12,'Interpreter','latex'); 
    ylabel('y','FontSize',12,'Interpreter','latex');
    grid on;
    for i =1: NE
        for j =1:numlocal 
            text(xnodes(NOP(i,j)),ynodes(NOP(i,j)),[num2str(NOP(i,j))],'FontSize',12,'Interpreter','latex');
        end;
    end; %for i
    print('-depsc',['CHEE662HW3af_meshused.eps'])
end % if 
    %OK so far 10.26


% %==Step 2: Define BCs&Initial Guess======================================
DBC = 0; %Dirichlet BC no slip
initguess = DBC;
us = ones(N,1); %initial nodal values == unknown basisfxn wts tb solved for
us = us*initguess; % << here define initial guess
%==========================================================================

disp(['linear basis fxns used, NE=',num2str(NE),' :'])

% SOLVE GFEM
% axpt = zeros(1,numlocal); aypt=zeros(1, numlocal); %empty vector to store asymmetric nodal coords 

NRiter=0; %N-R iteration counter
while err > convtol
NRiter = NRiter+1; %N-R iteration counter
Qin = 0; %reinitialize total flow rate

%Step 3: Calculate general wted residual vector and Jacobian of Stiffness wrt weights
Jij=zeros(N); %jacobian (==stiffness matrix if linear) tb recalculated for kth iteration 
Fi=zeros(N,1); %wted residual vector (==load matrix if linear) tb recalculated for kth iteration
for i = 1: NE %loop over all elements
    for j=1:4, %store the global nodes relavent to the ith element for later use
        W1(j)=NOP(i,j); 
%         axpt(j) = xnodes(W1(j));aypt(j) = ynodes(W1(j));
    end %for j: 
    % (dx, dy calculated here before)
    for j = 1:length(GP) %perform Gaussian integral over element i                 
        for k =1:length(GP)
        [phi, phic, phie] =  blinbasfxn(GP(j),GP(k),1, 1); % << now calcing derivatives wrt xi and eta
        %perform isoparametric mapping
            %calculate elements of Transformation Jacobian
            xc = 0; xe = 0; yc = 0; ye = 0;
            for h = 1:4
                xc = xc + xnodes(W1(h))*phic(h);
                xe = xe + xnodes(W1(h))*phie(h);
                yc = yc + ynodes(W1(h))*phic(h);
                ye = ye + ynodes(W1(h))*phie(h);
            end% for h
            detJt = xc*ye - xe*yc; %determinant
            for h = 1:numlocal
                phix(h) = (ye*phic(h)-yc*phie(h))/detJt;
                phiy(h) = (-xe*phic(h)+xc*phie(h))/detJt; %XXXXXXXXXXXXchanged(-xe*phic(h)-xc*phie(h))/detJt;
            end% for h
        
        %calculate function value at each of the GPs 
        u = 0; ux = 0; uy=0; %initialize basis fxn weights and basis fxn weight derivatives
        for L=1:numlocal %calculate FEM's approx function at relevant GP                           
            u = u+phi(L)*us(W1(L)); %FEM's approx function value at GP = sum(ui*phi)
            ux = ux+phix(L)*us(W1(L));
            uy = uy+phiy(L)*us(W1(L));            
        end %loop over local basis functions 
        
        for L = 1:numlocal %calculate using 2X1 and 2X2 Jij and Fi                          
            L1 = W1(L); %global node
            %in the following 2 loops the particular GFEM form of the actual GE has to be specified          
            %==============================================================
            Fi(L1) = Fi(L1)+detJt*W(j)*W(k)*(-Ri*(ux*phix(L) + uy*phiy(L))+(1*phi(L))); %$$$$ !! note + sign : for N-R must take negative later            
            for m = 1:numlocal 
                M1 = W1(m); %global node                                              
                Jij(L1, M1) = Jij(L1, M1)+detJt*W(j)*W(k)*(-Ri*(phix(m)*phix(L)+phiy(m)*phiy(L))-(0*phi(L))); %$$$$  
            end % for m    

            %=============================================================
        end % L: filling of 2X1 and 2X2 Jij and Fi                             
        Qin = Qin + detJt*W(j)*W(k)*u;
        end % k: Gaussian integration in x == zeta direction
   end % j : Gaussian integration in y == eta direction   
end %i : loop over all the elements
%Fi = Fi, Jij = Jij

%Step 4: impose Dirichlet and Neumann BCs on J & F
%=======================================================================
% fix Dirichlet BCs into our Galerkin equation
% ==>> 
%find nodes that lie on BC
BCnodeis = [find(xnodes == xmax) find(ynodes == ymin) find(xnodes == xmin) find(ynodes == ymax) find(Rnodes == Rmin)]; 
BCnodeis = sort(unique(BCnodeis));
for i = 1: length(BCnodeis) %apply BC at these nodes
    Fi(BCnodeis(i)) = 0; 
    Jij(BCnodeis(i),:) = zeros(1,N);
    Jij(BCnodeis(i),BCnodeis(i)) = 1;
    us(BCnodeis(i))= DBC;
end
%Fi=Fi, Jij = Jij

% % here, no Neumann BCs to our wted residual vector Fi

% %=========================================================================

%Step 5: invert to get next iteration's unknown basis nodal values == basis fxn weights 
A = -Jij\Fi; %inv(Jij)*(-Fi);
us = us + A; %output weights new iterative solution

%Step 6: calculate error == Euclid norm of diff of N-R iterations of C
err = sqrt(A'*A);
if NRiter> maxNR, break, end  
end % while

disp(['N-R iterations : ', num2str(NRiter),' ,convergence eror norm of : ',num2str(err)])%' , wted resid vector norm of ',num2str(err)])
if NRiter>maxNR, disp(['OBS: No N-R convergence in ', num2str(maxNR),' iterations']), end

%Step 7: Calculate Q
Q = 0;
for i = 1: NE %loop over all elements
    for j=1:4, %store the global nodes relavent to the ith element for later use
        W1(j)=NOP(i,j); 
    end %for j: 
    for j = 1:length(GP) %perform Gaussian integral over element i                 
            for k =1:length(GP)
            [phi, phic, phie] =  blinbasfxn(GP(j),GP(k),1, 1); % << now calcing derivatives wrt xi and eta
            %perform isoparametric mapping
                %calculate elements of Transformation Jacobian
                xc = 0; xe = 0; yc = 0; ye = 0;
                for h = 1:4
                    xc = xc + xnodes(W1(h))*phic(h);
                    xe = xe + xnodes(W1(h))*phie(h);
                    yc = yc + ynodes(W1(h))*phic(h);
                    ye = ye + ynodes(W1(h))*phie(h);
                end% for h
                detJt = xc*ye - xe*yc; %determinant            
            %calculate function value at each of the GPs 
            u = 0; ux = 0; uy=0; %initialize basis fxn weights and basis fxn weight derivatives
            for L=1:numlocal %calculate FEM's approx function at relevant GP                           
                u = u+phi(L)*us(W1(L)); %FEM's approx function value at GP = sum(ui*phi)
            end %loop over local basis functions         
            Q = Q + detJt*W(j)*W(k)*u;                                        
        end % k: Gaussian integration in x == zeta direction
   end % j : Gaussian integration in y == eta direction   
end %i : loop over all the elements
%Q = Q
%Qin = Qin
%end Step 7


%Step 8: calculate relative error
if iNEg ~= 1 %if iNEg >1, store nodal diffs, then calculate the relative error using Gaussian integration    
    usdiffs=[];
    for inode =1:N        
        %first test to see if nodes match
        xnodematch = find(oxnodes == xnodes(inode)); %possible xnodes that match
        ynodematch = find(oynodes == ynodes(inode));    
        a=size(xnodematch); 
        b= size(ynodematch);        
        if (a(2)~=0) && (b(2)~=0) %if there are matches in both x and y
            for ixmatch = 1:length(xnodematch) %check to see which one gives same node
                nodematch = find(xnodematch(ixmatch) == ynodematch); %index for which both x & y matches have same node #
                match=size(nodematch);            
                if match(2) %if match not empty, calculate difference btw solutions       
                    %relerr = relerr + (us(inode) - origus(nodematch))^2;
                    %store us
                    usdiffs = [usdiffs (us(inode)-origus(nodematch))];
                end % if
            end % for
        end % if     
    end % for inode
    relerrQ = sqrt((oQ-Q)^2);
    relerrsQ = [relerrsQ relerrQ];

    %calculate relative error now using above stored nodal diffs
    relerr = 0;
    for i = 1: oNE %loop over all elements
        for j=1:4, %store the global nodes relavent to the ith element for later use
            W1(j)=oNOP(i,j); 
        end %for j: 
        for j = 1:length(GP) %perform Gaussian integral over element i                 
                for k =1:length(GP)
                [phi, phic, phie] =  blinbasfxn(GP(j),GP(k),1, 1); % << now calcing derivatives wrt xi and eta
                %perform isoparametric mapping
                    %calculate elements of Transformation Jacobian
                    xc = 0; xe = 0; yc = 0; ye = 0;
                    for h = 1:4
                        xc = xc + oxnodes(W1(h))*phic(h);
                        xe = xe + oxnodes(W1(h))*phie(h);
                        yc = yc + oynodes(W1(h))*phic(h);
                        ye = ye + oynodes(W1(h))*phie(h);
                    end% for h
                    detJt = xc*ye - xe*yc; %determinant            
                %calculate function value at each of the GPs 
                usdiff = 0;%initialize fxn to calculate at specific (xi,eta)
                for L=1:numlocal %calculate FEM's approx function at relevant GP                           
                    usdiff = usdiff+phi(L)*usdiffs(W1(L)); %FEM's approx function value at GP = sum(ui*phi)
                end %loop over local basis functions         
                relerr = relerr + detJt*W(j)*W(k)*(usdiff^2);                                        
            end % k: Gaussian integration in x == zeta direction
       end % j : Gaussian integration in y == eta direction   
    end %i : loop over all the elements
    relerr = sqrt(relerr);
    relerr2 = sqrt(usdiff'*usdiff*detJt);
    relerr2s = [relerr2s relerr2]; 
    disp(['NEg : ', num2str(NEg),' ,Q = ', num2str(Q),', relerrQ = ',...
        num2str(relerrQ),', relerr = ', num2str(relerr),', relerr2 = ', num2str(relerr2),', last detJt = ', num2str(detJt)])
    
    
end %if iNEg
%store nodal coordinates and fxn values for next run
oQ = Q;
oxnodes = xnodes;
oynodes = ynodes;
origus = us;
oNE = NE;
oNOP = NOP;
disp([' '])

end %for iNEg

%Step 6: Plot Solution
%simple plot of basis wt fxns only
figure(2);
ti = xmin:((xmax-xmin)/meshnum):xmax; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(xnodes,ynodes,us,XI,YI);
contour(XI,YI,ZI,25);
hold on
%mesh(XI,YI,ZI), hold on
%plot3(xnodes,ynodes,us,'ok','MarkerSize',4)
%view(0,0)
maxua = max(max(ZI))
%maxua = max(us)
[imaxx, imaxy] = find(ZI == maxua); 
plot3(XI(imaxx,imaxy),YI(imaxx,imaxy),ZI(imaxx,imaxy),'kx','MarkerSize',12)
plot3(diag(XI),diag(YI),diag(ZI),':k','LineWidth',1.5)
text(XI(imaxx,imaxy)*1.1,YI(imaxx,imaxy),ZI(imaxx,imaxy),['$w_{z,max}$ = ', num2str(ZI(imaxx,imaxy),'%.3g')],'FontSize',10,'Interpreter','latex');
text(xmax*0.7,ymax*0.62,ZI(imaxx,imaxy)*0.5,['$Q_{v,tot}$ =', num2str(Q,'%.3g')],'FontSize',12,'Interpreter','latex');
circle([xmax,ymax],Rmin,50,'k-')

title(['GFEM solution for axial velocity, NER=',num2str(NER),', NE$\theta$ =',num2str(NET)],'FontSize',10,'Interpreter','latex');
xlabel('x/L','FontSize',12,'Interpreter','latex'); 
ylabel('y/L','FontSize',12,'Interpreter','latex');

view(0,90)
%print('-depsc',['CHEE662HW3ac_soln_3Dcontour.eps'])
print('-depsc',['CHEE662HW3af_soln_3Dcontour_',num2str(NEg),'.eps'])
%pause
view(62,4)
zlabel('w* (dimensionless axial velocity)','FontSize',11,'Interpreter','latex');
%print('-depsc',['CHEE662HW3ac_soln_3Dplot.eps'])
print('-depsc',['CHEE662HW3af_soln_3Dplot_',num2str(NEg),'.eps'])
%pause

%now plot to ensure solution is symmetric
%========================================
%choose theta s to plot to check symmetry of solution about pi/4 axis
iTs = (1:NET/4);
%========================================
diffiT = zeros(length(iTs), NR); 
for i=1:length(iTs)
    iTmirror = (NT-iTs(i)+1);
    nodes = (iTs(i)-1)*NR + (1:NR);
    nodesmirror = (iTmirror-1)*NR + (1:NR);
    diffiT(i,:) = us(nodes)-us(nodesmirror);
end;
figure(3);
plot((1:NR), diffiT,'--');
legend(num2str(iTs'),'Location','best');
title('GFEM solution differences for similar nodes','FontSize',10,'Interpreter','latex');
xlabel('LRHS - ULHS node solution','FontSize',12,'Interpreter','latex'); 
ylabel('Radial Node','FontSize',12,'Interpreter','latex');
grid on;
print('-depsc',['CHEE662HW3af_solnsymm.eps'])


%plot relative errors
figure(4);
loglog(NEgs(2:length(NEgs)), relerrsQ, 'o','MarkerSize',10); hold on;
loglog(NEgs(2:length(NEgs)), relerr2s, 'x','MarkerSize',10)
legend('Rel Error for Q','Rel Error for w s');
grid on;
title('GFEM Root-square of Euclidian Error Norm of w and Q ','FontSize',10,'Interpreter','latex');
xlabel('number of R or $\theta$ elements','FontSize',12,'Interpreter','latex'); 
ylabel('$\|e(Q)\|_2 \equiv \sqrt{(Q_v - Q_{v,prev})^2}$  or $\|e(w_z)\|_2$ ','FontSize',10,'Interpreter','latex');
print('-depsc',['CHEE662HW3af_errconv.eps'])


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
THETA=linspace(pi,1.5*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
H=plot(X,Y,style,'linewidth',1);
axis square;
end %circle