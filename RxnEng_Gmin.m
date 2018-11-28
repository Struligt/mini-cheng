function HW1Gmin
%Solves for eqm mole fractions by minimization of free energy


%define mass balance matrix A, where A.delN =! 0 
% [(SO2, CH4, S2, H2O, CO2, H2S, H2); (H,C,O,S)]

%comments: 
% - used logarithmic scales (here used in system of equations to solve) 
% since orders of variables vary wildly (as the mole fracns do here)

% - initial ratio info enters in the Aelem*yfracns = b (elemental abundance
%vector)

% - constraint of sum of mole fractions ==1 implicit

clear all; close all

global Gf Aelem R T b N M%P initratio
%parameters
conv = 4.1868; %J/cal : kcal to kJ
R = 8.314;
P = 1; %in atm
T0 = 298.15; 
Ts = [400, 600, 800];
%Ts = 400;

N = 7;
M = 4;
React = N-M;
A = zeros(M,N);

%********* input **********************************************
Asp = [' SO2       ', ' CH4       ', ' S2       ', ' H2O       ', ' CO2       ', ' H2S       ', ' H2       '];
Asp2 = ['SO2', 'CH4', 'S2', 'H2O', 'CO2', 'H2S', 'H2'];
Aelem = [0 4 0 2 0 2 2; 0 1 0 0 1 0 0; 2 0 0 1 2 0 0; 1 0 2 0 0 1 0]; %ie (H,C,O,S)
initratio = 2; %initial ratio of moles SO2 to CH4
b = [4 1 initratio*2 initratio*1]./(initratio+1) %b is elemental abundance vector ['H', 'C', 'O', 'S']

%from Perry's ChemEng Handbook:
Gf0 = [-71.68, (-5.049*1e-3*1e7/1e3)/conv, 19.36, -54.6351, -94.26, -7.85, 0]; %kcal/mol
Hf0 = [-70.94, (-7.452*1e-3*1e7/1e3)/conv, 31.02, -57.7979, -94.052, -4.77, 0];
Gf0 = conv*Gf0*1e3; %in J
Hf0 = conv*Hf0*1e3; %in J

% initial guess of mole fractions SINGLE PHASE 
if initratio == 2
    %logyinit = [-4 -5 -1 -2 -1 -5 -10]; %works for feed ratio = stoichiometric ratio 
    logyinit = [-1 -46 -5 -1 -1 -1 -19];            %OK i cheated here and used the soln using Kequil
elseif initratio < 2
    logyinit = [-35 -1 -28 -1 -1 -1 -8]; %works for feed ratio <> stoichiometric ratio with SO2 limiting reactant
else 
    logyinit = [-1 -100 -1 -2 -1 -5 -10];
end %if    
yinit = 10.^logyinit;

loggaminit = [1 1 1 1]*0; %initial guess of Lagrange multipliers
%**************************************************************
logsolinit = [logyinit, loggaminit]; %soln vector to solve for (N+M x 1)


% A = rref(Aelem); %reduced row echelon form for mass balance matrix
% muij = [-A(:,M+1:N); eye(React,React)]; %muij (NxR) is the table of stoichiometric coeffs for each reaction (colmn)
% Grxn0 = Gf * muij;
% Hrxn0 = Hf * muij
% lnKstd = - Grxn0/(R*T0)
% Kstd = exp(lnKstd);
% %assuming Enthalphy of formation does not change with temperature...
% lnKs = zeros(1, React); %lnK for given T

resultstore = zeros(length(Ts),N+3);

for i =1:1 %length(Ts)    
    %calculate partial molar Gibbs Free Energies at specified temperature
    Gf = Gf0 + (Gf0-Hf0)*log(Ts(i)/T0) %partial molar Gibbs free energy at specified temperature
    T = Ts(i);    
    %     lnKs  =  lnKstd - Hrxn0/R*(1/Ts(i) - 1/T0);    
    %     Ks = exp(lnKs);

    %use Matlab to solve for the equil concn         
    options=optimset('Display','final');%optimset('Display','iter'); % Option to display output
    [logyLM,fval] = fsolve(@eqnsetGFmin,logsolinit,options);     %[y,fval] = fsolve(@eqnset,yinit,options)  % Call optimizer
    
    yLM = 10.^logyLM
    
    disp(['For T = ',num2str(Ts(i)),', equilibrium gas mole fracns are: '])
    disp(num2str(yLM,'%.2e\t'))
    disp(['check: sum of mole fracns = ',num2str(sum(yLM),'%10.3e\n')])%'%10.3e\n')])
    disp(['relative amount of H2S to S2:',num2str(yLM(6)/yLM(3),'%10.5e\n')])
    fprintf('\n')
    %resultstore(i,:) = [ Ts(i), y, y(6)/y(3), sum(y)];
    
end % for
% 
% disp(['Results: , for init ratio SO2/CO2 = ', num2str(initratio)])
% resultcols = [' T(K)       ', Asp, '       H2S/S2 ', '       checksum '];
% disp(resultcols)
% disp(num2str(resultstore,'%.2e\t'))
% fprintf('\n')
% 
% Create pretty graph
% bar(log10(resultstore(:,2:1+N))')
% title(['Species Mole Fraction for Varying Temperature, InitRatio = ',num2str(initratio)],'FontSize',10,'Interpreter','latex');
% xlabel('Species','Interpreter','latex')
% ylabel('log10(Mole Fracn)','Interpreter','latex')
% for i=1:N        
%         text(i,10,[Asp2(i)], 'FontSize',10,'Interpreter','latex')  
% end
% grid on
% print('-depsc',['CHEE641HW1_RxnEqmBar_init',num2str(initratio),'.eps'])


end %function

function F = eqnsetGFmin(logyLM)
 %Gibbs Free Energy minimization
global Gf Aelem T b N M R %initratio muij 
gam = 10.^logyLM((N+1):(N+M)); %Lagrange multipliers
y = 10.^logyLM(1:N);%mole fracns
F = [ 
      Gf(1) + R*T*log(10)*logyLM(1) + gam*Aelem(:,1); %Lagrangian fxn minima wrt species moles
      Gf(2) + R*T*log(10)*logyLM(2) + gam*Aelem(:,2);
      Gf(3) + R*T*log(10)*logyLM(3) + gam*Aelem(:,3);
      Gf(4) + R*T*log(10)*logyLM(4) + gam*Aelem(:,4);
      Gf(5) + R*T*log(10)*logyLM(5) + gam*Aelem(:,5);
      Gf(6) + R*T*log(10)*logyLM(6) + gam*Aelem(:,6);
      Gf(7) + R*T*log(10)*logyLM(7) + gam*Aelem(:,7);
      Aelem(1,:)*y' - b(1);
      Aelem(2,:)*y' - b(2);
      Aelem(3,:)*y' - b(3);
      Aelem(4,:)*y' - b(4);    
      
      ];      

end %function
 

% function cumprodexp = cumprodexp(y,x)
% temp=1;
% for i =1:length(x)
%    temp=temp*(y(i)^x(i)); 
% end    
% cumprodexp=temp;    
% end %function
