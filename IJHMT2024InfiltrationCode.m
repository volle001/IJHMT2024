%Vaughan Voller, University of Minnesota, volle001@umn.edu
%August 2024

%Code for infiltration into a partially saturated porous media--
%provided as supporting material for the following paper

%Anomalous Infiltration in Partially Saturated Porous Media
%Vaughan Voller and Fabio D.A. Aarão Reis
%Submitted to Int. J. Heat Mass Transfer 
%This work was supported as part of the 
%Center on Geo-processes in Mineral Carbon Storage, 
%an Energy Frontier Research Center 
%funded by the U.S. Department of Energy, Office of Science, 
%Basic Energy Sciences at the University of Minnesota under
%award #DE-SC0023429.

%FDAAR acknowledges support from the Brazilian agencies 
%CAPES (PrInt 88887.310427/2018-00) 
%CNPq (442233/2023-0, 305391/2018-6), and 
%FAPERJ (E-26/210.040/2020, E-26/201.050/2022).

%The objective of the paper and this code is to investigate anomalous 
%infiltration in partially saturated porous media. The domain consists 
%of a one dimensional tube, on the order of 1 or more meters in length, 
%with a background saturated hydraulic conductivity K. Within this tube 
%there is fractal (Cantor set) distribution of 1cm thick lenses with a much 
%lower saturated hydraulic conductivity. Initially the tube has a constant   
%partially saturated moisture content, at a uniform pressure head h0<0.
%At time t=0, the entrance of the tube is saturated and maintained at a
%pressure head of h=0. The aim of this code is to determine how the 
%distribution of lenses effects the infiltration of moisture into the tube. 

%This code is set to predict the advance of the infiltration length,
%with time and the evolution of the moisture profile. 
%The default settings are set for an NV=3 Cantor set with a domain of size 
%~100cm. This is sufficient to contain the 4th order fractal.
%The predictions obtained match the prediction of the log-log plot of
%infiltration length with time (Figure (2) N=3 in the paper) and the 
%moisture profile when the infiltration has reached 81 cm, the end of the 4th 
%order fractal. 

%Note, by changing the appropriate default settings, this code can be 
%used to obtain the data for all of the results in the body of the paper. 

%The algorithm for solving the Richards equation in this code is 
%detailed in appendix of the paper. This algorithm, based directly 
%on the numerical method presented in Celia et al, 
%A general mass-conservative numerical solution for the unsaturated flow 
%equation. Water Resources Research, 26(7):1483-1496, 1990.

%NOTE TO AVOID CONFUSION IN UNIT CONVERSION IN CONSITITIVE EQUs
%CALCULATIONS ARE DONE IN cm 




clear all
NV=3;  %Cantor set type, Other choices are NV=2 (homogeneous), NV=4, NV=5   
posmax=4 ; %number of orders of the fractal in the domain
GS=16; %number of numerical control volumes in 1cm lens.

zend=-(NV^(posmax)); % length of the posmax order fractal
nfrac=-GS*zend-1; % nodes in fractal of order posmax 
delz=-zend/(nfrac+1); %grid size
delt=.5; %time step   defulat 0.5
n=nfrac+300; % Increase nodes beyound end of fractal order length, 
             % this alliviates the effect of the downstream boundary.    

theta=zeros(n,2); %moisture content (*,1) old, (*,2) new
h=zeros(n,1);     % initial head
z=zeros(n,1);     %locations of nodes

%constitutive paramters in cm  (Suggested in Celia)  
alpha=1.611*1e6; 
thetas=0.287;
thetar=0.075;
beta=3.96;   
gamma=4.74;
Abase=1.175*1e6;


%initial and boundary conditions
hzero=0; %Boundary head values ay z=0
hint=-100; %Initial condition
hend=hint; %Boundary value at end of domain
%Boundary moisture value at z=0 and z=n*delz
thzero=alpha*(thetas-thetar)./(alpha+abs(hzero).^beta)+thetar;    
thend=alpha*(thetas-thetar)./(alpha+abs(hend).^beta)+thetar;

%Variable values of Saturated Conductivity and 
%constitutive values A and alp--these can change between background and lens 
%NOTE, to alleviate effects of down stream end condition,  
%Values are calculated one extra fractal order

%Settings for fractal layer 
 K=ones(GS,1); %GS is the number of CV's in unit lens
 A=ones(GS,1);
 alp=ones(GS,1);
 %Values in spaces between lenses
 for or=1:posmax+1
      Kmid=500*ones(GS*(NV-2)*(NV^(or-1)),1); %500 x lens value
      K=[K;Kmid;K];
      Amid=1*ones(GS*(NV-2)*(NV^(or-1)),1);   % currently 1 x lens value
      A=[A;Amid;A];
      alpmid=1*ones(GS*(NV-2)*(NV^(or-1)),1); %currently 1 x lens value
      alp=[alp;alpmid;alp];
 end

%Set values of arrays holding nodal values of Ks, A, and alpha. 
Ks=0.00944*K;%cm/s  
A=Abase*A;
alp=alpha*alp;

%For Homogeneous Case Overwrite: 
%Ks=0.00944*ones(n+1,1);  %Require one extra storage
%A=Abase*ones(n+1,1);
%alp=alpha*ones(n+1,1);

%set moisture
for ii=1:n
    h(ii,1)=hint; 
    theta(ii,1)=alp(ii)*(thetas-thetar)./(alp(ii)+abs(h(ii,1)).^beta)+thetar;
    z(ii,1)=(-ii)*delz; %locations of nodes with unknowns 
end
theta(:,2)=theta(:,1);
thint=theta(:,1); %initial theta 

%Coefficients for linear solver (Ax=b)
%diffusion
aw=zeros(n,1);
ap=zeros(n,1);
ae=zeros(n,1);
%gravity NOT USED IN THIS CODE
awz=zeros(n,1);
aez=zeros(n,1);

tstep=1;  %time step counter

Fluxin=0;  %Moisture infiltrated in
Fluxout=0; %out


%moisture value that defines position of infiltration front  
dryval=thint(1,1)+(0.005)*(thetas-thint(1,1)); 

%Loop for calculating infiltration to end of fractal order NV
%Variable ending locations at EndLocation (in cm) 
%can be imposed by replacing (nfrac+1,2) with (EndLocation*GS,1)
while theta(nfrac+1,2)< dryval
    thediff=1; %convergence parameter
    iter=1; %iteration counter 
    while thediff>5e-8 
        
        for ii=1:n
            if ii==1
                hw=(hzero+h(ii,1))/2;
                kw=Ks(ii)*A(ii)/(A(ii)+(-hw)^gamma);
            else
                hw=(h(ii-1,1)+h(ii,1))/2;
                kw=Ks(ii)*A(ii)/(A(ii)+(-hw)^gamma);
            end
            
            if ii==n
                he=(h(ii,1)+hend)/2;
                ke=Ks(ii+1)*A(ii+1)/(A(ii+1)+(-he)^gamma);
            else
                he=(h(ii,1)+h(ii+1,1))/2;
                ke=Ks(ii+1)*A(ii+1)/(A(ii+1)+(-he)^gamma);
            end
            
            
            %diffusion coefficients
            
            aw(ii,1)=kw/delz;
            ae(ii,1)=ke/delz;
            ap(ii,1)=aw(ii,1)+ae(ii,1);
            
            %gravity coefficients NOT USED in this Code
            awz(ii,1)=kw;
            aez(ii,1)=ke;
            
            %derivative dtheta/dh
            cval(ii,1)=beta*alp(ii)*(thetas-thetar)...
                      *abs(h(ii,1))^(beta-1)/(alp(ii)+abs(h(ii,1))^beta)^2;
        end
        
        
        %Prepare TriDiagonal Matix Algorithm  
        cT(1:n)=aw(1:n);
        cT(1)=0; %breaks link with boundary
        bT(1:n)=ae(1:n);
        bT(n)=0;
        aT(1:n)=ap(1:n)+(delz/delt)*cval(1:n);
        dT(1:n)=(delz/delt)*(cval(1:n).*h(1:n));
% NOTE with gravity add the term         +(awz(1:n)-aez(1:n));
        ex(1:n)=(theta(1:n,1)-theta(1:n,2));
        dT(1:n)=dT(1:n)+(delz/delt)*ex(1:n);
        dT(1)=dT(1)+aw(1)*hzero; %fixed boundary at y=0
        dT(n)=dT(n)+ae(n)*hend; %fixed boundary at y=ell
        
        %TDMA SOLVER
        cTD(1)=bT(1)/aT(1);
        for ii=2:n-1
            cTD(ii)=bT(ii)/(aT(ii)-cT(ii)*cTD(ii-1));
        end
        dTD(1)=dT(1)/aT(1);
        for ii=2:n
            dTD(ii)=(dT(ii)+cT(ii)*dTD(ii-1))/(aT(ii)-cT(ii)*cTD(ii-1));
        end
        
        %Note unknowns are on node 2 to n-1 h(2:n-1)
        h(n,1)=dTD(n);  
        for ii=n-1:-1:1
            h(ii,1)=dTD(ii)+cTD(ii)*h(ii+1,1);
        end
        h(:,1)=min(h(:,1),0); %Limiter on head value
        
        %update moisture -- with under-relaxtion of 0.95
        thetapre=theta(:,2);
        for ii=1:n
            theta(ii,2)=theta(ii,2)+.95*(alp(ii)*...
                         (thetas-thetar)./(alp(ii)+abs(h(ii,1)).^beta)+...
                          thetar-theta(ii,2)); 
        end
        
        %convergence measure (maximum nodal difference of moisture update)            
        thediff=max(abs(theta(2:n-1,2)-thetapre(2:n-1)));
        
        %iteration count 
        iter=iter+1;
       
    end  %convergence in time step 
    
    
   
    
    t(tstep)=tstep*delt; %store time value
    
    %calculate moisture into  and out of domain
    Fluxin=Fluxin+delt*(aw(1,1)*(hzero-h(1,1)));
%Note for gravity add the term     +awz(1,1));
    Fluxout=Fluxout+delt*(ae(n,1)*(h(n,1)-hend));
%Note for gravity add the term     +aez(n,1)));
    %Infiltration Length 
    FF(tstep)=Fluxin;
    %Mass Balance infiltration by flux - infiltration by balance
    %This value is a measure of mass errors, should be small
    MBal(tstep)=sum(theta(:,2)-thint)*delz-Fluxin+Fluxout;
    
    
    theta(:,1)=theta(:,2);
    tstep=tstep+1;
    
end


%Mass error expressed as length 
disp('Mass Error ( should be small 1e-4 or less)')
disp(delz*sum(theta(:,2)-thint)-Fluxin+Fluxout)

%slope of log-log of infiltration length 
pv=polyfit(log(t),log(FF),1);
disp('slope of log-log infiltration length vs time')
disp(pv(1))
disp('If value is > 0.5 super diffusive')

pfit=pv(1)*log(t)+pv(2); %linear fit to log-log plot 

figure (1)
title('Moisture Profile N=3, 4th order fractal',...
       'Interpreter','latex','FontSize', 16)
hold on
plot([0;-1*z/100;-1*zend/100],[thzero;theta(:,1);thend],'k-','LineWidth',2)
xlabel('position $x$ (m)','Interpreter','latex')
ylabel('mositure $\theta$','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter='latex';
ax.FontSize = 18;


figure(2)

title('Moisture infiltration length, N=3, 4th order fractal',...
       'Interpreter','latex','FontSize', 16)
hold on
plot(log(t),log(FF), '-b','LineWidth', 1)
plot(log(t),pfit,'--k','LineWidth', 1 )
xlabel('$\log(t)$','Interpreter','latex')
ylabel('$\log(F)$','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter='latex';
 ax.FontSize = 18;
  legend('prediction','linear fit',...
      'Location','southeast', 'Interpreter','Latex','FontSize', 16)
  legend boxoff
   


        
