%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOUNDARY CONDITIONS
%
% This subroutine loads the acoustic boundary conditions at the inlet and 
% outlet of the geometry from files 'Inlet.txt' and 'Outlet.txt'
% respectively. The acoustic reflection coefficient at the inlet and outlet
% are then plotted and saved in the output folder.
% 
% Last update : 24/11/2020
% 
% Supported Boundary Conditions:
%
% Type 1 (Open end)
% Type 2 (Closed end) 
% Type 4 (Given amplitude and time delay) 
% Type 5 (Given amplitude and phase lag)
% Type 11 (Unflanged Levine-Schwinger by Norris-Sheng polynomial)
% Type 12 (Flanged Levine-Schwinger by Norris-Sheng polynomial)
%
% Unsupported Boundary Conditions (as of v2.0):
%
% Type 3 (Choked boundary condition)
% Type 6 (Polynomial with time delay from coefficients)
% Type 7 (Polynomial with time delay from user data)
% Type 8 (Heat Exchangers boundary condition)
% Type 9 (Gain and phase interpolated from data)
% Type 10 (Arbitrary function for gain and phase)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieving the data from the input file - Inlet

fprintf("Loading the acoustic boundary conditions... ");

filename1='./Inputs/Inlet.txt';
fid1=fopen(filename1);

C_title1= textscan(fid1, '%s', 4);             % read title
C_cell1  = textscan(fid1, '%f %f %f %f');      % read numeric data
fclose(fid1);

Inlet_type    = C_cell1{1}; % Inlet type              
Inlet_param1  = C_cell1{2}; % Parameter 1        
Inlet_param2  = C_cell1{3}; % Parameter 2
Inlet_param3  = C_cell1{4}; % Parameter 3

%% Retrieving the data from the input file - Outlet

filename2='./Inputs/Outlet.txt';
fid2=fopen(filename2);

C_title2 = textscan(fid2, '%s', 4);           % read title
C_cell2  = textscan(fid2, '%f %f %f %f');     % read numeric data
fclose(fid2);

Outlet_type   = C_cell2{1}; % Outlet type
Outlet_param1 = C_cell2{2}; % Parameter 1        
Outlet_param2 = C_cell2{3}; % Parameter 2
Outlet_param3 = C_cell2{4}; % Parameter 3

%% Retrieving the data from the input file - Scan range
% Used in the figures

filename='./Inputs/Scan_range.txt';
fid=fopen(filename);

C_title = textscan(fid, '%s', 6);             % read title
C_cell  = textscan(fid, '%f %f %f %f %f %f'); % read numeric data
fclose(fid);

FreqMin = C_cell{1}; % Minimum frequency, in Hertz              
FreqMax = C_cell{2};  % Maximum frequency, in Hertz 
FreqNum = C_cell{3};  % Number of frequencies

f_sample    = linspace(FreqMin, FreqMax, FreqNum);  % sampling frequency
s           = 2*pi*f_sample*1i;                     % s = i omega

%% Inlet setup

switch Inlet_type
    case 1                                      % open end, R=-1;
        R1h = @(sf) -1 .* ones(1,length(sf));
        
    case 2                                      % closed end, R=1; 
        R1h = @(sf) 1 .* ones(1,length(sf));
% 	 case 3                                     % choked (not implemented)
%        gamma           = gamma(1,1);
%        M1              = M_mean(1,1);
%        TEMP            = gamma*M1/(1+(gamma-1)*M1^2);
%        num1      = (1-TEMP)/(1+TEMP);
%        R1h = @(sf) num1 .* ones(1,length(sf));
    case 4                                      % time lag (in ms)
        tau_d1 = Inlet_param2./1000;
        R1h = @(sf) Inlet_param1.*exp(-sf.*tau_d1);
        
    case 5                                      % phase lag
        Phi = Inlet_param2;
        if Inlet_param3 ~= 0 % Use param3 == 1 to indicate phase in degrees
            Phi = Phi*pi/180;
        end
        R1h = @(sf) Inlet_param1.*exp(-Phi.*1i).*ones(1,length(sf));
        
    case 6                                      % Polynomial with time delay from coefficients
        C_inlet_BCpoly = getBCpoly(filename1, Inlet_param2, Inlet_param3, "POLYNOMIAL_COEFFICIENTS");
        
        num1 = cell2mat(C_inlet_BCpoly{1});
        den1 = cell2mat(C_inlet_BCpoly{2});
        tau_d1 = Inlet_param1./1000;
        
        R1h = @(sf) polyval(num1,sf)./polyval(den1,sf).*exp(-sf*tau_d1);

%     case 7                                    % Polynomial with time delay from user data (not yet implemented)

%     case 8                                    % Boundary condition for heat exchangers

    case 9                                      % Gain and phase interpolated from data
        C_inlet_BCdata = getBCdata(filename1, "GAIN_PHASE_DATA");

        In_freq      = C_inlet_BCdata(:,1); % Frequency
        In_gain      = C_inlet_BCdata(:,2); % Gain        
        In_phase     = C_inlet_BCdata(:,3); % Phase
        
        if Inlet_param1 ~= 0 % Use param1 == 1 to indicate phase in degrees
            In_phase = In_phase.*pi./180;
        end
        
        interp_type = 'makima';
        R1h = @(sf) interp1(In_freq, In_gain, imag(sf)./(2*pi), interp_type).*exp(-1i.*interp1(In_freq, In_phase, imag(sf)./(2*pi), interp_type));
        
    case 10                                     % Arbitrary function for gain and phase
        C_inlet_funcs = getBCfuncs(filename1, "GAIN_PHASE_FUNCTIONS");
        
        InletBCfunc_Gain = eval(C_inlet_funcs{1});
        InletBCfunc_Phase = eval(C_inlet_funcs{2});
        
        if Inlet_param1 ~= 0 % Use param1 == 1 to indicate phase in degrees
            Inlet_phase_mod = pi/180;
        else
            Inlet_phase_mod = 1;
        end
        
        R1h = @(sf) InletBCfunc_Gain(imag(sf)./(2*pi)).*exp(-1i.*InletBCfunc_Phase(imag(sf)./(2*pi)).*Inlet_phase_mod);
        
    case 11                                     % Unflanged Levine-Schwinger by Norris-Sheng polynomial
        a_LS = r_sample(1); % inlet radius
        R1h = @(sf) polyLevineSchwingerBC(sf, a_LS, c_mean(1,1));
        
    case 12                                     % Flanged Levine-Schwinger by Norris-Sheng polynomial
        a_LS = r_sample(1); % inlet radius
        R1h = @(sf) polyFlangedLevineSchwingerBC(sf, a_LS, c_mean(1,1));
        
end

%% Outlet setup

switch Outlet_type
    case 1                            % open end, R=-1;
        R2h = @(sf) -1 .* ones(1,length(sf));

    case 2                            % closed end, R=1; 
        R2h = @(sf) 1 .* ones(1,length(sf));
%     case 3                          % choked (not implemented)
%        gamma   = gamma(1,end);
%        M2      = M_mean(1,end);
%        TEMP    = (gamma-1)*M2/2;
%        num2  = (1-TEMP)/(1+TEMP);
%        R2 = @(sf) num2 .* ones(1,length(sf));
    case 4                            % time lag
        tau_d2 = Outlet_param2./1000;
        R2h = @(sf) Outlet_param1.*exp(-sf.*tau_d2);

    case 5                            % phase lag
        Phi = Outlet_param2;
        if Outlet_param3 ~= 0 % Use param3 == 1 to indicate phase in degrees
            Phi = Phi*pi/180;
        end
        R2h = @(sf) Outlet_param1.*exp(-Phi.*1i).*ones(1,length(sf));

    case 6                            % Polynomial with time delay from coefficients
        C_outlet_BCpoly = getBCpoly(filename2, Outlet_param2, Outlet_param3, "POLYNOMIAL_COEFFICIENTS");
        
        num2 = cell2mat(C_outlet_BCpoly{1});
        den2 = cell2mat(C_outlet_BCpoly{2});
        tau_d2 = Outlet_param1./1000;
        
        R2h = @(sf) polyval(num2,sf)./polyval(den2,sf).*exp(-sf*tau_d2);

%     case 7                          % Polynomial with time delay from user data (not yet implemented)

%     case 8                                    % Boundary condition for heat exchangers

    case 9                                      % Gain and phase interpolated from data
        C_outlet_BCdata = getBCdata(filename2, "GAIN_PHASE_DATA");

        Out_freq      = C_outlet_BCdata(:,1); % Frequency
        Out_gain      = C_outlet_BCdata(:,2); % Gain
        Out_phase     = C_outlet_BCdata(:,3); % Phase
        
        if Outlet_param1 ~= 0 % Use param1 == 1 to indicate phase in degrees
            Out_phase = Out_phase.*pi./180;
        end
        
        interp_type = 'makima';
        R2h = @(sf) interp1(Out_freq, Out_gain, imag(sf)./(2*pi), interp_type).*exp(-1i.*interp1(Out_freq, Out_phase, imag(sf)./(2*pi), interp_type));
        
    case 10                           % Arbitrary function for gain and phase
        C_outlet_funcs = getBCfuncs(filename2, "GAIN_PHASE_FUNCTIONS");
        
        OutletBCfunc_Gain = eval(C_outlet_funcs{1});
        OutletBCfunc_Phase = eval(C_outlet_funcs{1});
        
        fclose(fid);
        
        if Outlet_param1 ~= 0 % Use param1 == 1 to indicate phase in degrees
            Outlet_phase_mod = pi/180;
        else
            Outlet_phase_mod = 1;
        end
        
        R2h = @(sf) OutletBCfunc_Gain(imag(sf)./(2*pi)).*exp(-1i.*OutletBCfunc_Phase(imag(sf)./(2*pi)).*Outlet_phase_mod);

    case 11                           % Levine-Schwinger by Norris-Sheng polynomial
        a_LS = r_sample(end); % outlet radius
        R2h = @(sf) polyLevineSchwingerBC(sf, a_LS, c_mean(1,end));
            
    case 12                           % Flanged Levine-Schwinger by Norris-Sheng polynomial
        a_LS = r_sample(end); % outlet radius
        R2h = @(sf) polyFlangedLevineSchwingerBC(sf, a_LS, c_mean(1,end));

end

%% Plotting the reflection coefficients

% Setup
if DISP_FIGS
    fig4=figure('Name','Boundary Conditions - Reflection Coefficients');
    set(fig4, 'Position', [200 100 1200 800])

    ha = tight_subplot(2,2,[0.06 0.12],[.10 .10],[.10 .05]);

    R1 = R1h(s);
    R2 = R2h(s);

    % Modulus of the upstream reflection coefficient

    axes(ha(1)); 

    hold on
    hLine1=plot(f_sample,abs(R1),'-','color','b','Linewidth',2); 
    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
    xlabel(gca,'','Color','k','Interpreter','LaTex','FontSize',16);
    ylabel(gca,'Gain [-]','Color','k','Interpreter','LaTex','FontSize',16)
    %set(gca,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},'YAxisLocation','left','Color','w');
    %set(gca,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})
    if (Inlet_type==9) % plot interpolation points
        pts1=scatter(In_freq, In_gain, 80, 'g', 'LineWidth', 2);
        legend("Interpolated", "Data Points", 'Location', 'best');
    end

    % AMacL 07/2020
    Ylim=get(gca,'ylim');
    ymin = Ylim(1);
    ymax = Ylim(2);
    ndps = 2; % Max. number of decimal places
    if ymax-ymin < 0.1^ndps
        % Increments graph at min. resolution
        ymin = round(ymin - 2*0.1^ndps,ndps);
        ymax = round(ymax + 2*0.1^ndps,ndps);
        ytickformat(sprintf('%%0.%df',ndps))
    end
    ylim([ymin,ymax]);

    xticklabels('auto')
    yticklabels('auto')
    set(gca,'box','on','linewidth',2,'FontSize',16);
    title({'$\mathcal{R}_{inlet}$',''},'fontsize',20,'interpreter','latex')
    hold off

    % Phase of the upstream reflection coefficient

    axes(ha(3)); 

    hold on
    hLine2=plot(f_sample,unwrap(angle(R1),1.9*pi)./pi,'-','color', 'b','Linewidth',2);
    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
    xlabel(gca,'$f$ [Hz]','Color','k','Interpreter','LaTex','FontSize',16);
    ylabel(gca,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex','FontSize',16)
    %set(gca,'xlim',get(ha(1),'xlim'),'xTick',get(ha(1),'xTick'),'YAxisLocation','left','Color','w');
    if (Inlet_type==9) % plot interpolation points
        pts2=scatter(In_freq, unwrap(In_phase./-pi,1.9*pi), 80, 'g', 'LineWidth', 2);
        legend("Interpolated", "Data Points", 'Location', 'best');
    end
    xticklabels('auto')
    yticklabels('auto')
    set(gca,'box','on','linewidth',2,'FontSize',16);
    hold off

    % Modulus of the dowstream reflection coefficient

    axes(ha(2)); 

    hold on
    hLine1=plot(f_sample,abs(R2),'-','color','r','Linewidth',2); 
    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
    xlabel(gca,'','Color','k','Interpreter','LaTex','FontSize',16);
    ylabel(gca,'Gain [-]','Color','k','Interpreter','LaTex','FontSize',16)
    %set(gca,'xlim',[0 1000],'xTick',200:200:1000,'xticklabel',{},'YAxisLocation','left','Color','w');
    %set(gca,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})
    if (Outlet_type==9) % plot interpolation points
        pts3=scatter(Out_freq, Out_gain, 80, 'g', 'LineWidth', 2);
        legend("Interpolated", "Data Points", 'Location', 'best');
    end

    % AMacL 07/2020
    Ylim=get(gca,'ylim');
    ymin = Ylim(1);
    ymax = Ylim(2);
    ndps = 2; % Max. number of decimal places
    if ymax-ymin < 0.1^ndps
        % Increments graph at min. resolution
        ymin = round(ymin - 2*0.1^ndps,ndps);
        ymax = round(ymax + 2*0.1^ndps,ndps);
        ytickformat(sprintf('%%0.%df',ndps))
    end
    ylim([ymin,ymax]);
    xticklabels('auto')
    yticklabels('auto')
    set(gca,'box','on','linewidth',2,'FontSize',16);
    title({'$\mathcal{R}_{outlet}$',''},'fontsize',20,'interpreter','latex')
    hold off

    % Phase of the downstream reflection coefficient

    axes(ha(4)); 

    hold on
    hLine2=plot(f_sample,unwrap(angle(R2),1.9*pi)./pi,'-','color', 'r','Linewidth',2);
    set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
    set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
    xlabel(gca,'$f$ [Hz]','Color','k','Interpreter','LaTex','FontSize',16);
    ylabel(gca,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex','FontSize',16)
    %set(gca,'xlim',get(ha(1),'xlim'),'xTick',get(ha(1),'xTick'),'YAxisLocation','left','Color','w');
    if (Outlet_type==9) % plot interpolation points
        pts2=scatter(Out_freq, unwrap(Out_phase./-pi,1.9*pi), 80, 'g', 'LineWidth', 2);
        legend("Interpolated", "Data Points", 'Location', 'best');
    end
    xticklabels('auto')
    yticklabels('auto')
    set(gca,'box','on','linewidth',2,'FontSize',16);
    hold off

    % Saving the figure

    if SAVE_FIGS
        saveas(fig4,'./Outputs/Initialisation/Boundary_Conditions','fig')
    end
    if SAVE_PDFS
        save2pdf('./Outputs/Initialisation/Boundary_Conditions',fig4,300)
    end
    if SMALL_PLOTS
        for uu=1:4
            set(ha(uu), 'fontsize', get(gca,'DefaultAxesFontSize'));
        end
        set(fig4, 'Position', [200 100 800 600])
    end
end

fprintf("Done.\n ");

%% BC File IO Functions
function [D_mat] = getBCdata(fname, header_string)
    C     = {[]};
    D_mat = [];
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning("BC data could not be read from file: "+header_string+" not found");
    else
        while isempty(C{1}) && ~feof(fid)
            C = textscan(fid,'%f %f %f'); % Get data
            if isempty(C{1})
                fgets(fid); % skip line, try again
            end
        end
        if ~isempty(C{1})
            D_mat=cell2mat(C);
        else
            warning("BC data could not be read from file: data not found");
        end
    end
    fclose(fid);
end

function [C] = getBCfuncs(fname, header_string)
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning("BC functions could not be read from file: " + header_string + " not found");
    else
        C = {fgetl(fid),fgetl(fid)};
        if ~(C{1}(1)=='@') || ~(C{2}(1)=='@')
            warning("BC functions are not of expected format");
        end
    end
    fclose(fid);
end

function [C] = getBCpoly(fname, Nnum, Nden, header_string)
    C = {{[]},{[]}};
    L = "";
    pos = 0;
    fid = fopen(fname);
    while ~strcmp(fgetl(fid),header_string) && ~feof(fid)
    end
    if feof(fid)
        warning("BC polynomial coefficients could not be read from file: "+header_string+" not found");
    else
        % Numerator
        while isempty(C{1}{1}) && ~feof(fid)
            L = fgetl(fid);
            [C{1},pos] = textscan(L,'%f', Nnum+1); % Get data
        end
        if length(C{1}{1})<Nnum+1
            warning("Expected BC polynomial numerator coefficient(s) not found: order "+num2str(Nnum)+" specified so "+num2str(Nnum+1)+" coefficients expected");
        end
        if pos+1<length(L)
            D = textscan(L(pos+1:end),'%f');
            if ~isempty(D{1})
                warning(num2str(length(D{1}))+" unexpected BC polynomial numerator coefficient(s) found: order "+num2str(Nnum)+" specified so "+num2str(Nnum+1)+" coefficients expected");
            end
        end
        
        %Denominator
        while isempty(C{2}{1}) && ~feof(fid)
            L = fgetl(fid);
            [C{2},pos] = textscan(L,'%f', Nden+1); % Get data
        end
        if length(C{2}{1})<Nden+1
            warning("Expected BC polynomial denominator coefficient(s) not found: order "+num2str(Nden)+" specified so "+num2str(Nden+1)+" coefficients expected");
        end
        if pos+1<length(L)
            D = textscan(L(pos+1:end),'%f');
            if ~isempty(D{1})
                warning(num2str(length(D{1}))+" unexpected BC polynomial denominator coefficient(s) found: order "+num2str(Nden)+" specified so "+num2str(Nden+1)+" coefficients expected");
            end
        end
    end
    fclose(fid);
end

%% Levine-Schwinger functions

function [R] = polyLevineSchwingerBC(sf_LS, a_LS, c) % (sf = i omega values, a_LS = radius, c = soundspeed)
    % See Norris and Sheng 1989
    k_LS = imag(sf_LS) / c; % wavenumber
    He_LS = a_LS.*k_LS;
    alpha_LS = 0.2;
    beta_LS = -0.084;
    Rmod_LS = (1 + alpha_LS.*He_LS + beta_LS.*(He_LS.^2)) ./ (1 + alpha_LS.*He_LS + (0.5 + beta_LS).*(He_LS.^2));
    l_a_LS = (0.6133 + 0.027.*(He_LS.^2)) ./ (1 + 0.19.*(He_LS.^2));
    R = -1 .* Rmod_LS .* exp(-2i .* He_LS .* l_a_LS); % -ve exponent harmonic convention
end

function [R] = polyFlangedLevineSchwingerBC(sf_LS, a_LS, c) % (sf = i omega values, a_LS = radius, c = soundspeed)
    % See Norris and Sheng 1989
    k_LS = imag(sf_LS) / c; % wavenumber
    He_LS = a_LS.*k_LS;
    alpha_LS = 0.323;
    beta_LS = -0.077;
    Rmod_LS = (1 + alpha_LS.*He_LS + beta_LS.*(He_LS.^2)) ./ (1 + alpha_LS.*He_LS + (1 + beta_LS).*(He_LS.^2));
    l_a_LS = (0.82159 - 0.49.*(He_LS.^2)) ./ (1 - 0.46.*(He_LS.^3));
    R = -1 .* Rmod_LS .* exp(-2i .* He_LS .* l_a_LS); % -ve exponent harmonic convention
end

