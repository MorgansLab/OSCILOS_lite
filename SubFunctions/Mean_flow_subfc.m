%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MEAN FLOW
%
% This subroutine loads the mean flow parameters at the inlet of the geometry 
% from the file 'Mean_flow.txt'. The mean conservation equations are then 
% solved and the mean flow variables are obtained throughout the geometry.
% The mean velocity and mean temperature across the geometry are then
% plotted and saved in the output folder.
% 
% Last update : 07/01/2021
%
% Gamma is constant
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieving the data from the input file

fprintf("Loading the mean flow parameters... ");

filename='./Inputs/Mean_flow.txt';
fid=fopen(filename);

C_title         = textscan(fid, '%s', 6);           % read title
C_cell          = textscan(fid, '%f %f %f %f %f %f');        % read numeric data
fclose(fid);

P1      = C_cell{1}; % Mean pressure at the inlet        
T1      = C_cell{2}; % Mean temperature at the inlet           
index_M1_u1  = C_cell{3}; % Choice between Mach number (value: 1) or mean flow velocity (value: 2) at the inlet
M1_u1     = C_cell{4}; % Mach number or mean flow velocity at the inlet
index_gamma  = C_cell{5}; % Choice between a constant adiabatic index (value: 1) or changing with temperature (value: 2)
Delta_T_HS = C_cell{6}; % Temperature jump across the heat source (to be set to one if there is no heat source)

fprintf("Done.\n ");
fprintf("Computing the mean flow throughout the geometry... ");

%% Constant physical parameters

Ru       = 8.3145;
W_air    = 28.96512;
R_air    = Ru/W_air.*1000;
DeltaHr=0; % Initialisation of the heat source

%% Calculation of the mean flow - Inlet conditions

% The inlet gas is considered as air, since the gamma and W of the mixture
% is nearly the same as those of air, with an error of 5%.
% The burned gases are processed by the program to calculate the mass
% fraction of the mixture

gamma_cst = 1.4;
Cp_cst = gamma_cst ./ (gamma_cst - 1) .* R_air;

numSection    = length(x_sample)-1;                 % number of sections

% It is necessary to indicate that T_mean(a1:a2,b), herein, a1 denotes the 
% final state in the current interface, and a2 denotes the intermediate state
% b means the index number of the section

p_mean(1:2,1)     = P1;
T_mean(1:2,1)     = T1;

if index_gamma==1
    gamma(1:2,1)      = gamma_cst;
    Cp(1:2,1)         = Cp_cst;
    c_mean(1:2,1)     = (gamma(1,1)*R_air*T_mean(1,1)).^0.5;     
end

rho_mean(1:2,1)           = p_mean(1,1)./(R_air.*T_mean(1,1)); % mean density

% Zero-Mach warning - AMacL 03/08/2020

M_min = 1e-6;

if index_M1_u1==1
    if abs(M1_u1) < M_min
        warning("Near-zero inlet Mach number encountered ("+sprintf("%.1e",M1_u1)+") - M1 will be set to " + sprintf("%.1e",M_min));
        M_mean(1) = M_min;                % mean Mach number
    else
        M_mean(1) = M1_u1;                % mean Mach number
    end
    u_mean(1)     = M_mean(1).*c_mean(1); % mean velocity
        
elseif index_M1_u1==2
    if abs(M1_u1) < M_min * c_mean
        warning("Near-zero inlet velocity encountered ("+sprintf("%.1e",M1_u1)+") - u1 will be set to " + sprintf("%.1e",M_min * c_mean(1)));
        u_mean(1) = M_min * c_mean(1);    % mean velocity
    else
        u_mean(1) = M1_u1;                % mean velocity
    end
    M_mean(1)     = u_mean(1)./c_mean(1); % mean Mach number
        
else
    error("Value of index_M1_u1 ("+num2str(index_M1_u1)+") in "+filename+" unrecognised - supported values are 1 and 2");
end

%% Calculation of the mean flow - Different sections

for ss = 1:numSection-1 
    Theta(ss)  = (r_sample(ss+1)./r_sample(ss)).^2;       % Surface area ratio S2/S1
end

indexHA_num = 0;      % set the initial value to zero and it will be increased by 1 after every HA interface
for ss = 1:numSection-1 
    % At every interface, the changes are split into two steps:
    % 1. Cross sectional surface area change
    % 2. Heat addition
    %
    % --------------step 1-------------------------
    %
    [p_mean(1:2,ss+1),rho_mean(1:2,ss+1),u_mean(1:2,ss+1) ] = Fcn_calculation_TP_mean_WO_HeatAddition(p_mean(1,ss),rho_mean(1,ss),u_mean(1,ss),Theta(ss),gamma(1,ss));
    % ----------
    gamma(1:2,ss+1)   = gamma(1,ss);
    Cp(1:2,ss+1)      = Cp(1,ss);
    T_mean(1:2,ss+1)  = gamma(1,ss+1)/(gamma(1,ss+1)-1)*p_mean(1,ss+1)./(Cp(1,ss+1).*rho_mean(1,ss+1));                     
    c_mean(1:2,ss+1)  = ((gamma(1,ss+1) - 1).*Cp(1,ss+1).*T_mean(1,ss+1)).^0.5;
    M_mean(1:2,ss+1)  = u_mean(1,ss+1)./c_mean(1,ss+1);

    % --------------step 2-------------------------
    %
    if (SectionIndex(ss+1)==10) || (SectionIndex(ss+1)==11)             % with HA
            indexHA_num = indexHA_num + 1;              % this number is increased by 1
            % ---------first calculate the temperature, Cp, Rg after the heat addition interface -------------------------------------
            [TRatio(indexHA_num),c_mean(1,ss+1),DeltaHr(indexHA_num),...          % heat release rate per mass flow rate
                gamma(1,ss+1),Cp(1,ss+1)] =Fcn_GUI_INI_TP_calculation_products_after_HA(T_mean(2,ss+1),Delta_T_HS);
            T_mean(1,ss+1)    = TRatio(indexHA_num).*T_mean(2,ss+1);
            
            if index_gamma==1
                    gamma(1,ss+1)     = gamma_cst;
                    Cp(1,ss+1)        = Cp_cst;
                    c_mean(1,ss+1)    =  (gamma(1, 1)*R_air*T_mean(1,ss+1)).^0.5;        
            end
            
            Rg2 = Cp(1,ss+1)./(gamma(1,ss+1)./(gamma(1,ss+1)-1));
                    % ---------then, use the resolved temperature, Rg and the mean properties after the area changes to calculate the mean properties after HA -------------------
                    [p_mean(1,ss+1),rho_mean(1,ss+1),u_mean(1,ss+1)] = Fcn_calculation_TP_mean_W_HeatAddition(p_mean(2,ss+1),rho_mean(2,ss+1),u_mean(2,ss+1),Rg2,T_mean(1,ss+1));
                    M_mean(1,ss+1)  = u_mean(1,ss+1)./c_mean(1,ss+1);
            if index_gamma==1
                    DeltaHr(indexHA_num) = Cp(1,ss+1).*(T_mean(1,ss+1) - T_mean(2,ss+1)) + 0.5*(u_mean(1,ss+1).^2 - u_mean(1,ss).^2);              
            end
                    mass = rho_mean(2,ss+1).*u_mean(2,ss+1).*r_sample(ss+1).^2.*pi; % mass flow rate before HA
                    Q(indexHA_num)  = DeltaHr(indexHA_num).*mass;       % heat release rate                 
    end
end

%% Figures for the mean flow velocity and mean temperature

if DISP_FIGS
    fig2=figure('Name','Mean Flow Properties');
    set(fig2, 'Position', [150 75 1200 800])

    ha = tight_subplot(2,1,[0.15 0],[.10 .05],[.10 .05]);

    N = length(SectionIndex);
    
    x_plots(1,1:N-1) = x_sample(1:N-1);
    x_plots(2,1:N-1) = x_sample(2:N);

    for ss = 1:N-1
        ColorUDF{ss} = 'b';  % color of the line
    end
    
    if n_HS == 1
        indexHA=find(SectionIndex,11||10);
        for ss = indexHA(1):N-1
            ColorUDF{ss} = 'r';         % after the first heat addition interface, the color of plotted lines are set to red
        end
    end
  
    
    % Figure 1 : Mean flow velocity

    axes(ha(1)); 

    for ss = 1:N-1
        if TubeIndex(ss)==0
            y_plots(1:2,ss) = u_mean(1,ss);
        elseif TubeIndex(ss)==1 || TubeIndex(ss)==2
            y_plots(1:2,ss) = NaN;
        end   
    end

    hold on
    
    if n_HS == 1
        for ss = 1:length(indexHA)
            plot([x_plots(1,indexHA(ss)),x_plots(1,indexHA(ss))],[y_plots(1,indexHA(ss)-1),  y_plots(1,indexHA(ss))],'color','g','linewidth',2,'linestyle','--');
        end
    end 

    for ss = 1:N-1
        plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF{ss},'linewidth',2,'linestyle','-');
    end

    for ss = 1:N-1
        a = ss; 
        if TubeIndex(ss)==0
        elseif TubeIndex(ss)==1 || TubeIndex(ss)==2
            if ss == N-1
                y_plots(1:2,ss) = u_mean(1,ss);
            else
                y_plots(1:2,ss) = u_mean(1,ss:ss+1);
            end
            plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF{ss},'linewidth',2,'linestyle','-');
        end   
    end

    set(gca,'xlim',[x_sample(1), x_sample(end)],'YAxisLocation','left','Color','w');

    yvalue_max  = max(max(y_plots));
    yvalue_min  = min(min(y_plots));

    % This is my fix AMacL 07/2020
    ndps = 4; % Max. number of decimal places
    ymax = yvalue_max + 0.1*(yvalue_max-yvalue_min);
    ymin = yvalue_min - 0.1*(yvalue_max-yvalue_min);
    if ymax-ymin < 0.1^ndps
        % Increments graph at min. resolution
        ymin = round(yvalue_min - 2*0.1^ndps,ndps);
        ymax = round(yvalue_max + 2*0.1^ndps,ndps);
        ytickformat(sprintf('%%0.%df',ndps))
    end

    set(gca,'ylim',[ymin ymax])
    set(gca,'Position', [0.13 0.6 0.82 0.35]) % Leaves more room for yticks

    xlabel('$x$ [m]','Color','k','Interpreter','LaTex');
    xticklabels('auto')

    ylabel('$\overline{u}$ [m/s]','Color','k','Interpreter','LaTex');
    yticklabels('auto')

    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2)

    % Figure 2 : Mean temperature

    axes(ha(2));

    for ss = 1:N-1
        if TubeIndex(ss)==0
            y_plots(1:2,ss) = T_mean(1,ss);
        elseif TubeIndex(ss)==1 || TubeIndex(ss)==2
            y_plots(1:2,ss) = NaN;
        end  
    end

    hold on
    
    if n_HS == 1
        for ss = 1:length(indexHA)
            plot([x_plots(1,indexHA(ss)),x_plots(1,indexHA(ss))],[y_plots(1,indexHA(ss)-1),  y_plots(1,indexHA(ss))],'color','g','linewidth',2,'linestyle','--');
        end
    end 

    for ss = 1:N-1
        plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF{ss},'linewidth',2,'linestyle','-');
    end

    for ss = 1:N-1
        a = ss;
        if TubeIndex(ss)==0
        elseif TubeIndex(ss)==1 || TubeIndex(ss)==2
            if ss == N-1
                y_plots(1:2,ss) = T_mean(1,ss);
            else
                y_plots(1:2,ss) = T_mean(1,ss:ss+1);
            end
            plot([x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],'color',ColorUDF{ss},'linewidth',2,'linestyle','-');
        end   
    end  

    set(gca,'xlim',[x_sample(1), x_sample(end)],'YAxisLocation','left','Color','w');

    yvalue_max  = max(max(y_plots));
    yvalue_min  = min(min(y_plots)); 

    % This is my fix
    ndps = 4; % Max. number of decimal places
    ymax = yvalue_max + 0.1*(yvalue_max-yvalue_min);
    ymin = yvalue_min - 0.1*(yvalue_max-yvalue_min);
    if ymax-ymin < 0.1^ndps
        % Increments graph at min. resolution
        ymin = round(yvalue_min - 2*0.1^ndps,ndps);
        ymax = round(yvalue_max + 2*0.1^ndps,ndps);
        ytickformat(sprintf('%%0.%df',ndps))
    end

    set(gca,'ylim',[ymin ymax])
    set(gca,'Position', [0.13 0.1 0.82 0.35]) % Leaves more room for yticks

    xlabel('$x$ [m]','Color','k','Interpreter','LaTex');
    xticklabels('auto')

    ylabel('$\overline{T}$ [K]','Color','k','Interpreter','Latex');
    yticklabels('auto')

    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2);

    % Saving the figure

    if SAVE_FIGS
        saveas(fig2,'./Outputs/Initialisation/Mean_flow','fig')
    end
    if SAVE_PDFS
        save2pdf('./Outputs/Initialisation/Mean_flow',fig2,300)
    end
    if SMALL_PLOTS
        set(fig2, 'Position', [150 75 800 500])
        % Restore default font sizes
        set(ha(1), 'fontsize', get(gca,'DefaultAxesFontSize'));
        set(ha(2), 'fontsize', get(gca,'DefaultAxesFontSize'));
    end
end

fprintf("Done.\n ");

%% Functions used to compute the mean flow parameters throughout the domain

function [p_mean2,rho_mean2,u_mean2] =Fcn_calculation_TP_mean_WO_HeatAddition(p1,rho1,u1,Theta,gamma)
% This function is used to calculate the mean thermoproperties changing
% across the interface between two sections with different section surface
% area. 

options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
x0      = [1,1,1];  
F       = @(x)myfun_p_rho_u(x, p1,rho1,u1,Theta,gamma);
x2      = fsolve(F, x0, options);   % solve equation

iternum     = 0;
exitflag    = 0;
while exitflag > 4 || exitflag < 1 && iternum <20
    [x2, ~, exitflag]   = fsolve(F, x0, options);   % solve equation
    iternum = iternum + 1;
    x0      = x2;
end
if iternum>19
    h = msgbox('The mean flow calculation is not converged!'); 
end

p_mean2     = x2(1)*p1;
rho_mean2   = x2(2)*rho1;
u_mean2     = x2(3)*u1;

end

function F = myfun_p_rho_u(x, p1,rho1,u1,Theta,gamma)
% See previous function

p2      = x(1)*p1;
rho2    = x(2)*rho1;
u2      = x(3)*u1;

F(1) = Theta*rho2*u2 - rho1*u1;
if Theta>=1              % Area increasing
    F(2) = p2+rho2*u2^2 - (p1 + rho1*u1^2/Theta);
elseif Theta<1          % Area decreasing
    F(2) = p2*rho1^gamma - p1*rho2^gamma;
end
F(3) = gamma/(gamma-1)*(p2/rho2 - p1/rho1) + 0.5*(u2^2 - u1^2);
end

function [TRatio, cMean2, DeltaHr, gamma2, Cp2] = Fcn_GUI_INI_TP_calculation_products_after_HA(TMean1,Delta_T_HS)

% This function is used to calculate the mean properties after the heat
% addition interface
%
% indexHA_num is the order of heat addition interface from the inlet to the
% outlet 
% TMean1 denotes the incident mean temperature, 
% pMean1 denotes the incident mean pressure
% TRatio is the temperature jump ratio
% cMean2 denotes the speed of sound after the interface
% DeltaHr denotes the heat addition, or heat release rate from flame
% gamma2 denotes the specific heat ratio
% Cp2 denotes the heat capacity at constant pressure after the interface
%
% first created: 2014-12-03
% last modified: 2015-06-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)

TRatio  = Delta_T_HS;
TMean2  = TMean1*TRatio;
[cMean1,gamma1,Delta_hsMass1,cpMole1,cpMass,Rg1, meanMW]    = FcnCalHeatPropertiesAir(TMean1);
[cMean2,gamma2,Delta_hsMass2,cpMole2,Cp2,Rg2, meanMW]       = FcnCalHeatPropertiesAir(TMean2);
DeltaHr = Delta_hsMass2 - Delta_hsMass1;  

end

function [p_mean2,rho_mean2,u_mean2]= Fcn_calculation_TP_mean_W_HeatAddition(p_mean1,rho_mean1,u_mean1,Rg2,T_mean2)

% This function is used to calculate the mean thermoproperties changing
% across the flame front which is considered as compact
% The calculation  consists two steps
% calculate the properties after the flame, herein, considere there is
% section areas are constant which equals to that of section 2
% first created: 2014-04-01
% last modified: 2014-12-04
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
%-----------------------------step  ---------------------------------------

x1   = [p_mean1,rho_mean1,u_mean1];
T2   = T_mean2;
%--------------------------------------------------------------------------
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
x0 = x1;  % Make a starting guess at the solution
F       = @(x)myfun_W_HA(x,x1,Rg2,T2);
[x2,fval] = fsolve(F,x0,options); % Call solver
%
p_mean2     = x2(1);               
rho_mean2   = x2(2);
u_mean2     = x2(3);

end

function F = myfun_W_HA(x2,x1,Rg2,T2)
% p2    = x2(1);
% rho2  = x2(2);
% u2    = x2(3);
F(1) = x2(2)*x2(3) - x1(2)*x1(3);                               % mass conservation equation
F(2) = (x2(1)+x2(2)*x2(3)^2) - (x1(1)+x1(2)*x1(3)^2);           % momentum conservation equation
F(3) = x2(1) - x2(2)*Rg2*T2;                                    % gas law equation
end

function [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(T)
% FcnCalHeatPropertiesAir
% This function is used to calculate the thermodynamic properties of air
% [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(T)
% Example :
% [c,gamma,Delta_hsMass,cpMole,cpMass,Rg, meanMW] = FcnCalHeatPropertiesAir(300)
%
% Author: Jingxuan LI jingxuan.li@imperial.ac.uk
% last revision: 2015-06-02
%

load 'THERMO.mat' 
%
R       = THERMO.R;
W       = THERMO.Prod.W(7);
Rg      = R./W.*1000;
[Delta_hsMole, Delta_hsMass] = FcnCalDelta_hs(T,1,7,THERMO);
[cpMole,cpMass] = FcnCalCp(T,1,7,THERMO);
gamma   = cpMole./(cpMole - R);
c       = (gamma.*Rg.*T).^0.5;
meanMW  = W;
end

function [Delta_hsMole, Delta_hsMass, hsMole, hsMass] = FcnCalDelta_hs(T,flagGasType,index,THERMO)

% FcnCalDelta_hs  --- Calculation of sensible enthalpy
% T:            temperature
% flagGasType:  type of gas. 
%               1: combustion products. 2: fuel
% index:        index of selected type of fuel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 fuels can be precribed:
% CH4,
% C2H4,
% C2H6,
% C3H8,
% 1-butene (C4H8),
% n-butane (C4H10),
% isobutane (C4H10),
% Jet-A(C12H23, gas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 products can be prescribed:
% CO2
% CO
% H2O
% H2
% O2
% N2
% Air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMO: data for calculation
% 
%
% >>  [Delta_hsMole, Delta_hsMass, hsMole, hsMass] = FcnCalDelta_hs(2000,1,1,THERMO)
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk) and Aimee S.Morgans (a.morgans@imperial.ac.uk)
%
% last revision: 2015-06-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R   = THERMO.R;           % R = 8.3145
T0  = THERMO.T0;          % T0 = 298.15 
if T<=1000
    indexN=2*index-1;
elseif 1000<T
    indexN=2*index;
end
if flagGasType==1
        aT0 = THERMO.Prod.a(2*index-1,:);
        aT  = THERMO.Prod.a(indexN,:);
        bT0 = THERMO.Prod.b(2*index-1,:);
        bT  = THERMO.Prod.b(indexN,:);
        W   = THERMO.Prod.W(index);
elseif flagGasType==2
        aT0 = THERMO.Fuel.a(2*index-1,:);
        aT  = THERMO.Fuel.a(indexN,:);
        bT0 = THERMO.Fuel.a(2*index-1,:);
        bT  = THERMO.Fuel.a(indexN,:);
        W   = THERMO.Fuel.W(index);
end
%
% Delta_hs/RT = -a1 T^-2 + a2 lnT T^-1 + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hs at input temperature
temp = 0;
for ss=1:5
    temp = temp + aT(ss+2).*T.^ss./ss;
end
hs = (temp - aT(1)./T + aT(2).*log(T) + bT(1)).*R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hs at standard temperature
temp=0;
for ss=1:5
    temp=temp+aT0(ss+2).*T0.^ss./ss;
end
hs0             = (temp - aT0(1)./T0 + aT0(2).*log(T0) + bT0(1)).*R;
Delta_hsMole    = hs - hs0;
hsMole          = hs0;
Delta_hsMass    = Delta_hsMole/W*1e3;
hsMass          = hsMole/W*1e3;
end

function [cpMole,cpMass] = FcnCalCp(T,flagGasType,index,THERMO)
% FcnCalCp  --- Calculation of heat capacity at constant pressure
% T:            temperature
% flagGasType:  type of gas. 
%               1: combustion products. 2: fuel
% index:        index of selected type of fuel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 fuels can be precribed:
% CH4,
% C2H4,
% C2H6,
% C3H8,
% 1-butene (C4H8),
% n-butane (C4H10),
% isobutane (C4H10),
% Jet-A(C12H23, gas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 products can be prescribed:
% CO2
% CO
% H2O
% H2
% O2
% N2
% Air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMO:           data for calculation
% 
%
% >>  [cpMole,cpMass] = FcnCalCp(2000,1,1,THERMO)
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk) and Aimee S.Morgans (a.morgans@imperial.ac.uk)
%
% last revision: 2015-06-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = THERMO.R;

if      T <= 1000
    indexN=2*index-1;
elseif  T > 1000
    indexN=2*index;
end

if flagGasType== 1
        aT  = THERMO.Prod.a(indexN,:);
        W   = THERMO.Prod.W(index);
elseif flagGasType==2
        aT  = THERMO.Fuel.a(indexN,:);
        W   = THERMO.Fuel.W(index);
end
%         
% Cp/Ru = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
%
cpMole = R./T.^2.*polyval(aT(7:-1:1),T);
cpMass = cpMole./W.*1e3;
end
