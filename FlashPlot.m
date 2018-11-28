function FlashPlot(T)
%uses FlashProj function
% eg FlashPlot((-75-32)/1.8) if T is in Farheneit 
%( function [Vfracf, xif, Kif] = FlashProj4(P,T,zi)  )
%define array of input feeds, at given temp and pressure:

hold off

%graphing colors and linestyles
cv = ['r' 'g' 'b' 'm' 'k'];
colorv = [cv cv cv]; %color call for plots which should contain enough colors for lines to be plotted
stylev = ['-' '--' ':' '-.'];
stylev = [stylev stylev stylev];

%==============================inputs====================================
Pv=[100 200 400 600 800]; %in psia
nC=3; %must agree with FlashProj
dz=0.10;
Compi={'Ethane', 'Methane', 'Propane'};
Ploop =length(Pv);

%load experimental data files
x1=importdata('x1.mat');
Kbyx1=importdata('Kbyx159C.mat');

%feed fractions define
zbase=[(1*dz):dz:(1-1.5*dz)];
zM=ones(length(zbase),nC);
zM(:,1)=zbase;
zM(:,2)=0.225;
zM(:,3)=abs(1-zM(:,2)-zM(:,1));

%normalize
for i=1:length(zM(:,1))
    zM(i,:)=zM(i,:)/sum(zM(i,:));
end
zMorig=zM

%initialize plots
for i=1:nC
    figure(i); clf
    %title('Flash Calculation for %d,' Compi)
    title(['Flash Calculation for ', Compi(i)],'FontSize',12); %insert title   
    xlabel('Liquid Molar Fraction of Ethane');
    ylabel('K-factor : y(i)/x(i) for component i');
    grid on;
    
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontsize',12);
    set(0,'defaultaxeslinewidth',1);
    set(0,'defaultlinelinewidth',1.5);
    set(gca,'Color','w','ycolor','k'); %background color white and y axes black
    set(get(gca,'Ylabel'),'Color','k'); %ylabel to black
    set(gcf,'Color','w');%background color white

end %for

for j=1:Ploop
    P= Pv(j);
%for each pressure, run flash calculation
%size(zM)
Vfracm=[];
xim=[];
Kim=[];
    if j==3 %adjust feed fraction to allow convergence - fudge
        zM(:,2)=0.225+0.20;
        zM(:,3)=abs(1-zM(:,2)-zM(:,1));
        for i=1:length(zM(:,1)) %normalize
           zM(i,:)=zM(i,:)/sum(zM(i,:));
        end
    end %if j==3
    if j==4 %revert to before
        zM=zMorig;
    end %if j==4
    
    for i=1:length(zM(:,1))
        [Vfracf,xif, Kif] = FlashProj(P,T,zM(i,:));
        Vfracm=[Vfracm; Vfracf];
        xim=[xim; xif];
        Kim=[Kim; Kif];
    end %for i
sprintf('Results for pressure (psia) of %d :', P)
Vfracm=Vfracm
xim=xim
Kim=Kim

%plot data
    for i=1:nC
        figure(i);
        hold on
        xlength=length(xim(:,1))-j+1; %contracted to match look of P-R paper : can readily be removed
        plot(xim(1:xlength,1),Kim(1:xlength,i),':','color',colorv(j),'linewidth',2);        
    end %for i
end %for j

for j=1:Ploop
    for i=1:nC %show legend
        figure(i);
        legend('100 psia', '200 psia', '400 psia', '600 psia', '800 psia');
    end %for i
end %for j

lenx1=length(x1(1,:));
for j=1:1 %plot imported data for P=100psia and T=-75F
    for i=1:nC 
        figure(i);
        plot(x1(j,:),Kbyx1(j,((i-1)*lenx1+1):(i*lenx1)),'o','color',colorv(j)),end%if
    end %for i
  end %for j