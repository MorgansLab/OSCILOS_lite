%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FREQUENCY-DOMAIN SOLVER
%
% This subroutine loads the scan range and number of initialisations - fed to
% the frequency-domain solver - from the file 'Scan_range.txt'. The solver
% then finds the eigenvalues of the system, and the corresponding
% eigenmodes. The eigenvalues are saved as an output file called
% 'Eigenvalues.txt'. A map of the eigenvalues in the range of interest is
% then plotted and saved in the output folder. A number of eigenmodes,
% represented as the modulus of the acoustic pressure and velocity, are
% also plotted and saved in the output folder.
% 
% Last update : 18/01/2020 
%
% No indirect noise
% Linear solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieving the data from the input file - Scan range

fprintf("Loading the solver parameters... ");

filename='./Inputs/Scan_range.txt';
fid=fopen(filename);

C_title         = textscan(fid, '%s', 6);           % read title
C_cell          = textscan(fid, '%f %f %f %f %f %f');        % read numeric data
fclose(fid);

FreqMin      = C_cell{1}; % Minimum frequency, in Hertz              
FreqMax      = C_cell{2};  % Maximum frequency, in Hertz 
FreqNum = C_cell{3};  % Number of frequencies

GRMin      = C_cell{4}; % Minimum growth rate, in Hertz              
GRMax      = C_cell{5};  % Maximum growth rate, in Hertz 
GRNum = C_cell{6};  % Number of growth rates

fprintf("Done.\n ");
fprintf("Initialising the solver... ");

%% Pre-processing

C1 =  [1 1 0;1 -1 0;1 1 0];
C2 =  [1 1 0;1 -1 0;1 1 -1];

index1=1;
index2=1;

% Definition of the acoustic transfer matrices for each section

for ss = 1:numSection-1 % Loop on the sections
    if SectionIndex(ss+1)==0
        % Definition of the physical parameters
        cRatio  = c_mean(index2,ss+1)./c_mean(index1,ss);
        M1      = M_mean(index1,ss);
        M2      = M_mean(index2,ss+1);
        gamma1  = gamma(index1,ss);
        gamma2  = gamma(index2,ss+1);
        Thetaf   = Theta(ss);
        rho1    = rho_mean(index1,ss);
        rho2    = rho_mean(index2,ss+1);

        % First matrix
        B1_temp(1,1) =   0;
        B1_temp(1,2) =   cRatio;
        B1_temp(1,3) =   M1.*cRatio;
        B1_temp(3,1) =   M1.*gamma1./(gamma1-1);
        B1_temp(3,2) =   M1.^2;
        B1_temp(3,3) =  -M1./(gamma1-1);

        % Second matrix
        B2_temp(1,1) =   0;
        B2_temp(1,2) =   Thetaf;
        B2_temp(1,3) =   Thetaf.*M2;
        B2_temp(3,1) =   cRatio.* Thetaf.*M2.*gamma2./(gamma2-1);
        B2_temp(3,2) =   cRatio.* Thetaf.*M2.^2;
        B2_temp(3,3) =  -cRatio.* Thetaf.*M2./(gamma2-1); 

        if Thetaf>=1 % Area expansion
            B1_temp(2,1) =   Thetaf;
            B1_temp(2,2) = 2*M1;
            B1_temp(2,3) =   M1.^2;
            B2_temp(2,1) =   Thetaf;
            B2_temp(2,2) = 2*Thetaf.*M2;
            B2_temp(2,3) =   Thetaf.*M2.^2;
        elseif Thetaf<1 % Area contraction
            B1_temp(2,1) = 1/rho1.^gamma1;
            B1_temp(2,2) = 0;
            B1_temp(2,3) =-1/rho1.^gamma1;
            B2_temp(2,1) = 1/rho2.^gamma2;
            B2_temp(2,2) = 0;
            B2_temp(2,3) =-1/rho2.^gamma2;
        end  

        B1{1,ss}         = B1_temp;          % first step
        B2{1,ss}         = B2_temp;
        B1{2,ss}         = eye(3);           % second step
        B2{2,ss}         = eye(3);
        BC{ss}           = (B2_temp*C2)\(B1_temp*C1);
        
    elseif SectionIndex(ss+1)==10||11
%        indexHA = indexHA + 1;
        
        % First step
        
        % Definition of the physical parameters       
        cRatio  = c_mean(2,ss+1)./c_mean(1,ss);
        M1      = M_mean(1,ss);
        M2      = M_mean(2,ss+1);
        gamma1  = gamma(1,ss);
        gamma2  = gamma(2,ss+1);
        Thetaf   = Theta(ss);
        rho1    = rho_mean(1,ss);
        rho2    = rho_mean(2,ss+1);

        % First matrix        
        B1_temp(1,1) =   0;
        B1_temp(1,2) =   cRatio;
        B1_temp(1,3) =   M1.*cRatio;
        B1_temp(3,1) =   M1.*gamma1./(gamma1-1);
        B1_temp(3,2) =   M1.^2;
        B1_temp(3,3) =  -M1./(gamma1-1);

        % Second matrix
        B2_temp(1,1) =   0;
        B2_temp(1,2) =   Thetaf;
        B2_temp(1,3) =   Thetaf.*M2;
        B2_temp(3,1) =   cRatio.* Thetaf.*M2.*gamma2./(gamma2-1);
        B2_temp(3,2) =   cRatio.* Thetaf.*M2.^2;
        B2_temp(3,3) =  -cRatio.* Thetaf.*M2./(gamma2-1); 

        if Thetaf>=1
            B1_temp(2,1) =   Thetaf;
            B1_temp(2,2) = 2*M1;
            B1_temp(2,3) =   M1.^2;
            B2_temp(2,1) =   Thetaf;
            B2_temp(2,2) = 2*Thetaf.*M2;
            B2_temp(2,3) =   Thetaf.*M2.^2;
         elseif Thetaf<1
            B1_temp(2,1) = 1/rho1.^gamma1;
            B1_temp(2,2) = 0;
            B1_temp(2,3) =-1/rho1.^gamma1;
            B2_temp(2,1) = 1/rho2.^gamma2;
            B2_temp(2,2) = 0;
            B2_temp(2,3) =-1/rho2.^gamma2;
        end

        B1{1,ss} = B1_temp;              
        B2{1,ss} = B2_temp;
        
        % Second step
        
        cRatio  = c_mean(1,ss+1)./c_mean(2,ss+1);
        M1      = M_mean(2,ss+1);
        M2      = M_mean(1,ss+1);
        gamma1  = gamma(2,ss+1);
        gamma2  = gamma(1,ss+1);
        HA      = 1.00*DeltaHr./c_mean(2,ss+1).^2;

        B1a_temp(1,1) =   0;
        B1a_temp(1,2) =   cRatio;
        B1a_temp(1,3) =   M1.*cRatio;
        B1a_temp(2,1) =   1;
        B1a_temp(2,2) = 2*M1;
        B1a_temp(2,3) =   M1.^2;
        B1a_temp(3,1) =   M1.*gamma1./(gamma1-1);
        B1a_temp(3,2) =   M1.^2 - HA;
        B1a_temp(3,3) =  -M1.*(1/(gamma1-1) + HA);
        
        B2_temp(1,1) =   0;
        B2_temp(1,2) =   1;
        B2_temp(1,3) =   M2;
        B2_temp(2,1) =   1;
        B2_temp(2,2) = 2*M2;
        B2_temp(2,3) =   M2.^2;
        B2_temp(3,1) =   cRatio.*M2.*gamma2./(gamma2-1);
        B2_temp(3,2) =   cRatio.*M2.^2;
        B2_temp(3,3) =  -cRatio.*M2./(gamma2-1); 
        
        B1{2,ss} = B1a_temp;              % second step
        B2{2,ss} = B2_temp;
        BC{ss} = eye(3);
    end  
end

x_diff = transpose(diff(x_sample));
tau_plus      = x_diff./(c_mean(1,:)+u_mean(1,:));
tau_minus     = x_diff./(c_mean(1,:)-u_mean(1,:));
tau_c         = x_diff./ u_mean(1,:);       % convection time delay

fprintf("Done.\n ");
fprintf("Computing the eigenvalues... ");

%% Eigenvalues calculation

% Split the scan domain

Scan_GRSp        = linspace( GRMin, GRMax, GRNum);
Scan_FreqSp      = linspace( FreqMin, FreqMax, FreqNum);

% Split the contour domain

Cont_GRSp        = linspace( GRMin, GRMax, 10*GRNum);
Cont_FreqSp      = linspace( FreqMin, FreqMax, 10*FreqNum);

% Linear solver

% Initialisating the waitbar     
wx = 0;
wb = waitbar(wx,sprintf("%d%%",wx),'Name','Calculating eigenvalues');
nss = length(Scan_GRSp);
nkk = length(Scan_FreqSp);
tnum = nss*nkk;
nxi = 100; % Number of 'increments' on waitbar
nxx = 0;

eigen_num=1;
options = optimset('Display','off');        % the calculation results by fsolve will not be shown in the workspace
for ss = 1:length(Scan_GRSp)
    GR = Scan_GRSp(ss); % Current growth rate
    for kk = 1:length(Scan_FreqSp)
        omega = 2*pi*Scan_FreqSp(kk); % Current frequency
        if GR == 0
            GR = 1; % Minimal growth rate
        end
        if omega == 0
            omega =10; % Minimal angular frequency
        end
        s0 = GR+1i*omega;                                              % initial value                                                     
        [x,fval,exitflag] = fsolve(@(s)Fcn_DetEqn_Linear(s,R1h,R2h,M_mean,gamma,numSection,SectionIndex,tau_plus,tau_minus,tau_c,Flame_param1,Flame_param2,DeltaHr,c_mean,Theta,B1,B2,C1,C2,BC),s0,options);   % linear                          
        if exitflag==1
            eigenvalue(eigen_num) = x;
            eigenvalue_prec(eigen_num) = floor(eigenvalue(eigen_num)./10).*10;      % this is used to set the precision
            eigen_num = eigen_num+1;
        end
        wx = (nkk*(ss-1) + kk - 1)/(tnum-1);
        if (wx >= nxx/nxi)
            waitbar(wx,wb,sprintf("%d%%",floor(wx*100)));
            nxx = nxx + 1;
        end
    end
end
close(wb); % Close waitbar

fprintf("Done.\n ");
fprintf("Processing the eigenvalues... ");

% AMacL 12/08/2020 re-expressed 3 operations to use logical indexing (previously loops)
% Substituted EigValCol{1} for eigval_unique

% Parameters (AMacL 12/08/2020 - previously hard-coded)
eig_prec = 10;     % Precision to which eigenvalues are evaluated [Re 1/s, Imag rad/s]
eig_min_diff = 1; % Minimum spacing (absval) between consecutive eigenvalues before they are considered identical
eig_min_freq = 1;  % Minimum eigenvalue frequency considered [rad/s]

% Getting rid of duplicate values
[eigval_unique_prec,m,n] = unique(eigenvalue_prec); 
eigval_unique = eigenvalue(m); % these are the eigenvalues
eigval_unique = sort(eigval_unique);

% Getting rid of near-zero values
eigval_unique = eigval_unique( abs(imag(eigval_unique)) >= eig_min_freq );

% Getting rid of redundant values
eigenvalue_unique_diff = diff(eigval_unique);
eigval_unique = [eigval_unique( abs(eigenvalue_unique_diff) >= eig_min_diff ), eigval_unique(end)]; % this leaves the eigenvalues we want

% Getting rid of values out of the range of interest
eigval_unique = eigval_unique( ~( real(eigval_unique)<GRMin | real(eigval_unique)>GRMax...
                                  | abs(imag(eigval_unique))<(FreqMin*2*pi) | abs(imag(eigval_unique))>(FreqMax*2*pi) ) );

if isempty(eigval_unique)
    error("Solver error: No meaningful eigenvalues were found."...
        + newline + "This can happen if:"...
        + newline + "No modes lie between min_freq and max_freq (try increasing max_freq)"...
        + newline + "No modes lie between min_GR and max_GR (try widening GR range)"...
        + newline + "The phase of a boundary reflectance is too strong a function of frequency (check BCs)");
end

%% Store the data in the output file called 'Eigenvalues.txt'

data_num(:,1)   = abs(imag(eigval_unique)./2./pi);
data_num(:,2)   = real(eigval_unique);
%data_cell       = num2cell(data_num);
number=transpose([1:1:length(data_num(:,1))]);
data=transpose([number data_num]);

if SAVE_EIGS
    filename2='./Outputs/Results/Eigenvalues.txt';
    fid2=fopen(filename2,'w');
    fprintf(fid2, 'Mode number \t Frequency [Hz] \t Growth rate [1/s] \n');
    fprintf(fid2, '%1.0f \t \t %f \t \t %f \n', data);
    fclose(fid2);
end

fprintf("Done.\n ");

%% Contour plot

if DISP_FIGS
    
    fprintf("Plotting the eigenvalues... ");
    
    % Initialising the waitbar
    
    wx = 0;
    wb = waitbar(wx,sprintf("%d%%",wx),'Name','Calculating contour plot');
    nuu = length(Cont_GRSp);
    nkk = length(Cont_FreqSp);
    tnum = nuu*nkk;
    nxi = 100; % Number of 'increments' on waitbar
    nxx = 0;
    
    F   = zeros(length(Cont_GRSp),length(Cont_FreqSp)); % Used for the contour plot
    for uu = 1:length(Cont_GRSp)
        GR = Cont_GRSp(uu);
        for kk  = 1:length(Cont_FreqSp)
            omega       = 2*pi*Cont_FreqSp(kk);
            s           = GR+1i*omega;
            F(uu,kk) = Fcn_DetEqn_Linear(s,R1h,R2h,M_mean,gamma,numSection,SectionIndex,tau_plus,tau_minus,tau_c,Flame_param1,Flame_param2,DeltaHr,c_mean,Theta,B1,B2,C1,C2,BC);                            
        end
        wx = (nkk*(uu-1) + kk - 1)/(tnum-1);
        if (wx >= nxx/nxi)
            waitbar(wx,wb,sprintf("%d%%",floor(wx*100)));
            nxx = nxx + 1;
        end
    end
    close(wb); % Close waitbar
    
    ValCol{1}=F; % This matrix contains the contour plot data

    fig4=figure('Name','Eigenvalue Map');
    set(fig4, 'Position', [250 125 800 800])
    position_hAxes1=get(gca,'position');

    hold on
    contourf(Cont_GRSp./100,Cont_FreqSp,20*log10(abs(ValCol{1}'))); % B.B. 05/07/2019 Suppressed output
    plot(real(eigval_unique)./100,abs(imag(eigval_unique))./2./pi,'p','markersize',8,'color','k','markerfacecolor',[1,1,1],'LineWidth',0.5);
    hold off
    % AMacL 15/07/2020 have added this in but eigenvalues are also missed if
    % scan numbers are too low
    grid on

    set(gca,'YColor','k','Box','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2)
    xlabel(gca,  '$ \textrm{Growth rate}/100~~$ [$1/s$] ','Color','k','Interpreter','LaTex','FontSize',16);
    ylabel(gca,'$ \textrm{Frequency}~~$ [Hz]','Color','k','Interpreter','LaTex','FontSize',16);
    set(gca,'ylim',[FreqMin FreqMax],'YAxisLocation','left','Color','w');
    set(gca,'xlim',[GRMin   GRMax]./100);
    set(gca,'layer','top');
    colorbar 
    colormap(hot);
    hcb=colorbar;
    set(hcb,'Fontsize',16,'box','on','Unit','points','location','eastoutside')

    hTitle = title('Eigenvalues are located at minima');
    set(hTitle, 'interpreter','latex', 'fontunits','points','fontsize',20)

    set(hcb,'Fontsize',16,'box','on','Unit','points','location','eastoutside') % Workaround to reinitialize the colorbar

    % Saving the figure
    if SAVE_FIGS
        saveas(fig4,'./Outputs/Results/Eigenvalues_map','fig')
    end
    if SAVE_PDFS
        save2pdf('./Outputs/Results/Eigenvalues_map',fig4,300)
    end
    if SMALL_PLOTS
        set(fig4, 'Position', [250 125 600 600])
        set(hTitle, 'fontsize', 11)
        set(hcb, 'fontsize', get(gca, 'DefaultTextFontSize'))
        set(gca, 'fontsize', get(gca,'DefaultAxesFontSize'));
    end
    
    fprintf("Done.\n ");

end

%% Figures corresponding to the different eigenmodes

if DISP_FIGS
    
    fprintf("Plotting the mode shapes... ");

    n_plot_mod=min(PLOT_MODES,length(number)); 
    for yy=1:n_plot_mod
        numMode=number(yy);
        s_star = eigval_unique(numMode);             % eigenvalue
        
        [x_resample,p,u] = Fcn_calculation_eigenmode_Linear(s_star,R1h,tau_plus,tau_minus,tau_c,BC,c_mean,u_mean,rho_mean,numSection,SectionIndex,Flame_param1,Flame_param2,DeltaHr,Theta,B1,B2,C1,C2,x_sample);

        fig5(yy)=figure('Name',"Mode Shape "+num2str(numMode)+" ("+num2str(sprintf('%.2f',abs(imag(eigval_unique(numMode))./2./pi)))+" Hz)",'NumberTitle','off');
        set(fig5(yy), 'Position', [280+20*yy 150 1200 800])
        ha = tight_subplot(2,1,[0.05 0],[.10 .10],[.10 .05]);

        axes(ha(2)); 

        hold on

        for lll=1:length(x_sample)-1
            plot(x_resample(lll,:),abs(p(lll,:)),'-','color','k','Linewidth',2)
        end

        valueLevel = round(log10(max(max(abs(p)))));
        ymax1=ceil(max(max(abs(p)))./10^valueLevel).*10^valueLevel;
        ymin1=floor(min(min(abs(p)))./10^valueLevel).*10^valueLevel;
        ylimitUD=[ymin1 ymax1+0.25*(ymax1-ymin1)];
        ytickUD=linspace(ylimitUD(1),ylimitUD(2),6);

        for zz=1:length(ytickUD)
            yticklabelUD{zz}=num2str(ytickUD(zz));
        end

        yticklabelUD{end}='';
        xmax1=max(max(x_resample));
        xmin1=min(min(x_resample));
        xlimitUD=[xmin1 xmax1];
        xtickUD=linspace(xlimitUD(1),xlimitUD(2),6);

        set(gca,'YColor','k','Box','on');
        set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2)
        xlabel('$x $ [m]','Color','k','Interpreter','LaTex','FontSize',16);
        ylabel('$|~\hat{p}~|$ [Pa] ','Color','k','Interpreter','LaTex','FontSize',16);
        set(gca,'ylim', ylimitUD,'yTick',ytickUD,'yticklabel',yticklabelUD,'YAxisLocation','left');
        set(gca,'xlim',xlimitUD,'xtick',xtickUD);
        ylimit=get(gca,'ylim');
        yticklabels('auto');
        xticklabels('auto');
        set(gca,'layer','top');

        NSp = length(x_sample);
        if NSp < 10 %&& Outlet_type ~= 8 % B.B. 07/07/2019 - added test to avoid section lines when using HX (reserved to indicate HX extent)
            for k=2:length(x_sample)-1
               plot([x_sample(k),x_sample(k)],ylimit,'--','linewidth',0.5,'color','k') 
            end
        % NOTE : k=1:length(x_sample) was changed to k=2:length(x_sample)-1 to remove the dotted lines placed on the y-axis     
        end % B.B. 05/07/2019 STOP
        grid on
        hold off

        axes(ha(1)); 

        hold on
        for k=1:length(x_sample)-1
            plot(gca,x_resample(k,:),abs(u(k,:)),'-','color','k','Linewidth',2)
        end
        set(gca,'YColor','k','Box','on');
        set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',2)
        set(gca,'xlim',xlimitUD,'xtick',xtickUD,'xticklabel',[]);
        xlabel(gca,'','Color','k','Interpreter','LaTex','FontSize',16);
        ylabel(gca,'$|~\hat{u}~|$ [m/s] ','Color','k','Interpreter','LaTex','FontSize',16);
        ylimit=get(gca,'ylim');
        yticklabels('auto');
        set(gca,'layer','top');

        if NSp < 10 && Outlet_type ~= 8 % B.B. 07/07/2019 - added test to avoid section lines when using HX (reserved to indicate HX extent)
            for k=2:length(x_sample)-1
                plot([x_sample(k),x_sample(k)],ylimit,'--','linewidth',1,'color','k') 
            end
        % NOTE : k=1:length(x_sample) was changed to k=2:length(x_sample)-1 to remove the dotted lines placed on the y-axis        
        end   % B.B. 05/07/2019 STOP
        grid on
        hold off

        str = [['Mode ',num2str(numMode),' - Frequency = ',num2str(sprintf('%.2f',abs(imag(eigval_unique(numMode))./2./pi))), ' Hz'];""];
        hTitle_bis = title(str);
        set(hTitle_bis, 'interpreter','latex', 'fontunits','points','fontsize',20)

        % Saving the figure

        path_solver=['./Outputs/Results/Mode_',num2str(numMode)];
        if SAVE_FIGS
            saveas(fig5(yy),path_solver,'fig')
        end
        if SAVE_PDFS
            save2pdf(path_solver,fig5(yy),300)
        end
        if SMALL_PLOTS
            set(fig5(yy), 'Position', [280+20*yy 150 800 500])
            set(ha(1), 'fontsize', get(gca,'DefaultAxesFontSize'));
            set(ha(2), 'fontsize', get(gca,'DefaultAxesFontSize'));
        end

    end
    
    fprintf("Done.\n ");

end

%% Functions

function F = Fcn_DetEqn_Linear(s,R1h,R2h,M_mean,gamma,numSection,SectionIndex,tau_plus,tau_minus,tau_c,Flame_param1,Flame_param2,DeltaHr,c_mean,Theta,B1,B2,C1,C2,BC)
% Defines the equation to be solved by the linear solver

    % Boundary conditions - Pressure acoustic reflection coefficients
    R1_temp      = R1h(s);
    R2_temp      = R2h(s);
    Rs      = -0.5*M_mean(end)./(1 + 0.5*(gamma(end) - 1 ).*M_mean(end)); % B.B. 10/07/2019

    Te = 0; % No indirect noise
    
%    indexHA = 0;            % index of heat addition
%    indexHP = 0;            % index of heat perturbation

    G = eye(3);
    
    for ll = 1:numSection-1 
        D1 = diag([ exp(-s*tau_plus(ll)), exp( s*tau_minus(ll)), exp(-s*tau_c(ll))]);
        
%        if Te==0 % AMacL 07/2020 - remove -> Inf contributions from exp(-s*tau_c) for low mean flow
%            D1(3,3) = 0;
%        end
        
        if SectionIndex(ll+1)==0
            Z = BC{ll}*D1;
        elseif SectionIndex(ll+1)==10
%            indexHA = indexHA + 1;
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = B1{2,ll}*(B2{1,ll}\B1{1,ll}) + B2b;
            BC1     = Bsum*C1;
            BC2     = B2{2,ll}*C2;
            Z       = (BC2\BC1)*D1;            
         elseif SectionIndex(ll+1)==11
%            indexHA = indexHA + 1;
%            indexHP = indexHP + 1;
            FT = Flame_param1.*exp(-s.*Flame_param2.*1e-3);
            B2b     = zeros(3);
            % in case there are two heat addition, but the first one is a
            % steady one and the second one is unsteady
            B2b(3,2)= DeltaHr./c_mean(2,ll+1)./c_mean(1,ll)./Theta(ll).*FT;
            Bsum    = B1{2,ll}*(B2{1,ll}\B1{1,ll}) + B2b;
            BC1     = Bsum*C1;
            BC2     = B2{2,ll}*C2;
            Z       = (BC2\BC1)*D1;
        end   
        
        G = Z*G;
    end

    A1_minus        = 1000;
    A1_plus         = R1_temp.*A1_minus;
    E1              = 0;
    Array_LeftBD    = [A1_plus, A1_minus, E1].'; % Bug fix 31/07/2020 .' is elementwise transpose not complex transpose

    D1End           = diag([exp(-s*tau_plus(end)), exp( s*tau_minus(end)), Te.*exp(-s*tau_c(end))]);
    
%    if Te==0 % AMacL 07/2020
%        D1End(3,3) = 0; % Prevent Te.*exp(-s*tau_c) going -> 0*Inf (=NaN) limit for very large tau_c caused by low mean flow
%    end

    Array_RightBD   = D1End*G*Array_LeftBD;
    AN_plus         = Array_RightBD(1);
    AN_minus        = Array_RightBD(2);
    EN_plus         = Array_RightBD(3);

    F = (R2_temp.*AN_plus + Rs.*EN_plus) - AN_minus;  

end

function [x_resample,p,u]=Fcn_calculation_eigenmode_Linear(s_star,R1h,tau_plus,tau_minus,tau_c,BC,c_mean,u_mean,rho_mean,numSection,SectionIndex,Flame_param1,Flame_param2,DeltaHr,Theta,B1,B2,C1,C2,x_sample)
    % This function is used to plot the modeshape of selected mode, 
    % Linear

    R1          = R1h(s_star);
    Te          = 0;

    A_minus(1)  = 1;
    E(1)        = 0;
    A_plus(1)   = R1.*A_minus(1);
    Array(:,1)  = [A_plus(1),A_minus(1),E(1)].'; % B.B. 05/07/2019 - Suppressed output

    for pp = 1:numSection-1
        D1 = diag([ exp(-s_star*tau_plus(pp)), exp( s_star*tau_minus(pp)), exp(-s_star*tau_c(pp))]);
        %if E(1)==0 % AMacL 07/2020 - remove -> Inf contributions from exp(-s*tau_c) for low mean flow
        %    D1(3,3) = 0;
        %end
        if SectionIndex(pp+1)==0
            Z{pp} = BC{pp}*D1;
        elseif SectionIndex(pp+1)==10
            %indexHA = indexHA + 1;
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = B1{2,pp}*(B2{1,pp}\B1{1,pp}) + B2b;
            BC1     = Bsum*C1;
            BC2     = B2{2,pp}*C2;
            Z{pp}       = (BC2\BC1)*D1;            
        elseif SectionIndex(pp+1)==11
            %indexHA = indexHA + 1;
            %indexHP = indexHP + 1;
            FT = Flame_param1.*exp(-s_star.*Flame_param2.*1e-3);
            B2b     = zeros(3);
            temp    = D1*Array(:,pp);
            uf      = abs((temp(1) - temp(2))./(c_mean(1,pp).*rho_mean(1,pp)));         % velocity before the flame
            B2b(3,2)= DeltaHr./c_mean(2,pp+1)./c_mean(1,pp)./Theta(pp).*FT;
            Bsum    = B1{2,pp}*(B2{1,pp}\B1{1,pp}) + B2b;
            BC1     = Bsum*C1;
            BC2     = B2{2,pp}*C2;
            Z{pp}   = (BC2\BC1)*D1;            
        end
        Array(:,pp+1) = Z{pp}*Array(:,pp); % B.B. 05/07/2019 - output suppressed
    end

    for k = 1:numSection-1
        A_plus(k+1)     = Array(1,k+1);
        A_minus(k+1)    = Array(2,k+1);
        E(k+1)          = Array(3,k+1);
    end

    for k=1:length(x_sample)-1
        x_resample(k,:)=linspace(x_sample(k),x_sample(k+1),201);    % resample of x-coordinate
    end

    c_mean      = c_mean(1,:);
    u_mean      = u_mean(1,:);
    rho_mean    = rho_mean(1,:);

    for k = 1:length(x_sample)-1
        kw1_plus(k)  = s_star./(c_mean(k)+u_mean(k));
        kw1_minus(k) = s_star./(c_mean(k)-u_mean(k));

        p(k,:) =    A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))+...
                    A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1)));
        u(k,:) =    (A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))-...
                    A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1))))./rho_mean(k)./c_mean(k);
    end
end
