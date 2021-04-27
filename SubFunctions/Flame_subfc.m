%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FLAME MODEL
%
% This subroutine loads the flame model from the file 'Flame.txt'. The gain
% and phase of the corresponding Flame Model are then plotted and saved in
% the output folder. If heat perturbations are not present inside the
% domain, this step has no impact on the final result.
% 
% Last update : 08/01/2021
% 
% Supported Flame models:
%
% n-tau model (Type 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If heat perturbations are absent, the n-tau model is used with a zero gain and phase
% This does not affect the final result

Flame_type    = 1; % n-tau model           
Flame_param1  = 0; % Gain set to zero        
Flame_param2  = 0; % Phase set to zero

%% Checking if there is a source of heat perturbations

if sum(~(SectionIndex-11))
    
    %% Retrieving the data from the input file - Inlet

    fprintf("Loading the flame model... ");

    filename1='./Inputs/Flame.txt';
    fid1=fopen(filename1);

    C_title1= textscan(fid1, '%s', 3);             % read title
    C_cell1  = textscan(fid1, '%f %f %f');      % read numeric data
    fclose(fid1);

    Flame_type    = C_cell1{1}; % Inlet type              
    Flame_param1  = C_cell1{2}; % Parameter 1        
    Flame_param2  = C_cell1{3}; % Parameter 2

    %% Retrieving the data from the input file - Scan range

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
    FTF = polyval(Flame_param1,s).*exp(-s.*Flame_param2.*1e-3);

    %% Plotting the FTF

    % Setup
    if DISP_FIGS
        fig3=figure('Name','Flame Transfer Function');
        set(fig3, 'Position', [200 100 1200 800])

        ha = tight_subplot(2,1,[0.06 0.12],[.10 .10],[.10 .05]);

        % Modulus of the FTF

        axes(ha(1)); 

        hold on
        hLine1=plot(f_sample,abs(FTF),'-','color','b','Linewidth',2); 
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
        xlabel(gca,'','Color','k','Interpreter','LaTex','FontSize',16);
        ylabel(gca,'Gain [-]','Color','k','Interpreter','LaTex','FontSize',16)
        %set(gca,'ylim',[0 1.2],'yTick',0:0.2:1.2,'yticklabel',{'',0.2:0.2:1})

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

        set(gca,'xticklabel',{});
        yticklabels('auto')
        set(gca,'box','on','linewidth',2,'FontSize',16);
        %title({'Flame Transfer Function',''},'fontsize',20,'interpreter','latex')
        hold off

        % Phase of the FTF

        axes(ha(2)); 

        hold on
        hLine2=plot(f_sample,unwrap(angle(FTF),1.9*pi)./pi,'-','color', 'b','Linewidth',2);
        set(gca,'YColor','k','Box','on','ygrid','on','xgrid','on');
        set(gca,'FontName','Helvetica','FontSize',16,'LineWidth',1)
        xlabel(gca,'$f$ [Hz]','Color','k','Interpreter','LaTex','FontSize',16);
        ylabel(gca,'Phase/$\pi$ [-]','Color','k','Interpreter','LaTex','FontSize',16)
        %set(gca,'xlim',get(ha(1),'xlim'),'xTick',get(ha(1),'xTick'),'YAxisLocation','left','Color','w');
        
        xticklabels('auto')
        yticklabels('auto')
        set(gca,'box','on','linewidth',2,'FontSize',16);
        hold off

        % Saving the figure

        if SAVE_FIGS
            saveas(fig3,'./Outputs/Initialisation/Flame_Model','fig')
        end
        if SAVE_PDFS
            save2pdf('./Outputs/Initialisation/Flame_Model',fig3,300)
        end
        if SMALL_PLOTS
            for uu=1:2
                set(ha(uu), 'fontsize', get(gca,'DefaultAxesFontSize'));
            end
            set(fig3, 'Position', [200 100 800 600])
        end
    end

    fprintf("Done.\n ");
end