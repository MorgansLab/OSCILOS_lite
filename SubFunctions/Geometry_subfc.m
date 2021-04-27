%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GEOMETRY
%
% This subroutine loads the geometry from the file 'Geometry.txt', plots it
% and saves it in the output folder.
% 
% Last update : 07/01/2021
% 
% No Liners, no Helmholtz resonators, A single heat source only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieving the data from the input file

fprintf("Loading the geometry... ");

filename='./Inputs/Geometry.txt';
fid=fopen(filename);

C_title         = textscan(fid, '%s', 4);           % read title
C_cell          = textscan(fid, '%f %f %f %f');        % read numeric data
fclose(fid);

x_sample      = C_cell{1}; % Axial references                   
r_sample      = C_cell{2}; % Radius              
SectionIndex  = C_cell{3}; % Section Index     
TubeIndex     = C_cell{4}; % Tube Index 

%% Updating the geometry for a gradually changing cross section area

% Error message if there is more than one heat source
n_HS=sum(~(SectionIndex-11)+~(SectionIndex-10));
n_FHS=sum(~(SectionIndex-11));

if n_HS > 1
    error("Geometry file: Too many heat sources (currently " ...
    + num2str(n_HS) ...
    + ") - Only one heat source is allowed in this version of OSCILOS_lite.")
end

% Warning if final TubeIndex is nonzero
if TubeIndex(end)~=0
    warning("Geometry file: Final section should have TubeIndex 0 (currently " ...
    + num2str(TubeIndex(end)) ...
    + ") - this will be taken as zero.")
    TubeIndex(end)=0;
end

% Warning of unsupported section indices
indexUS = find(SectionIndex~=0 & SectionIndex~=10 & SectionIndex~=11);
indexSHS = find(SectionIndex==10); % Is there a steady Heat Source?
indexFHS = find(SectionIndex==11); % Is there a fluctuating Heat Source?

if ~isempty(indexUS)
    warning("Geometry file: SectionIndex(es) " + num2str(SectionIndex(indexUS)) ...
        + " unsupported in this version and will be taken as zero.")
    SectionIndex(indexUS) = 0;
end

% Warning if only 2 sections provided (3 currently required)
if (length(x_sample)<2)
    error("Geometry file: insufficient number of sections to define geometry"...
        + newline + "Only " + num2str(length(x_sample)) + " provided.")
elseif (length(x_sample)==2 && TubeIndex(1)==0)
    warning("Geometry file: only 2 sections provided - solver requires 3."...
        + newline + "A 3rd section will be interpolated at the midpoint with TubeIndex matching the inlet.")
    x_sample(3) = x_sample(2);
    r_sample(3) = r_sample(2);
    SectionIndex(3) = SectionIndex(2);
    TubeIndex(3) = TubeIndex(2);
    
    x_sample(2) = 0.5*(x_sample(1)+x_sample(3));
    r_sample(2) = 0.5*(r_sample(1)+r_sample(3));
    SectionIndex(2) = 0;
    TubeIndex(2) = TubeIndex(1);
end

indexGC   = find(TubeIndex >= 1);
isNoGC    = isempty(indexGC);

if isNoGC==0
NumSplit  = 50;       % default splitting number
NumSplit = NumSplit.*ones(1,length(indexGC));
for ss = 1:length(indexGC)
    k = indexGC(ss);
    N = length(x_sample); % get the current length of x_sample
    x_sample((k+1:N) + NumSplit(ss)-2)        = x_sample(k+1:N); % move forward by NumSplit(ss)-1
    r_sample((k+1:N) + NumSplit(ss)-2)        = r_sample(k+1:N); % move forward by NumSplit(ss)-1
    SectionIndex((k+1:N) + NumSplit(ss)-2)    = SectionIndex(k+1:N);
    TubeIndex((k+1:N) + NumSplit(ss)-2)       = TubeIndex(k+1:N);
    
    x_sample(k:k+NumSplit(ss)-1) = linspace(x_sample(k),x_sample(k+NumSplit(ss)-1), NumSplit(ss));
    r_sample(k:k+NumSplit(ss)-1) = linspace(r_sample(k),r_sample(k+NumSplit(ss)-1), NumSplit(ss));
    SectionIndex(k+1:k+NumSplit(ss)-2) = 0;
    TubeIndex(k:k+NumSplit(ss)-2) = TubeIndex(indexGC(ss));
    
    if ss < length(indexGC)
        indexGC(ss+1:end) = indexGC(ss+1:end) + NumSplit(ss)-2;
    end    
end
end

%% Plotting the geometry

if DISP_FIGS
    fig1=figure('Name','Geometry');
    set(fig1, 'Position', [100 50 1200 400])
    hold on

    WallLineWidth = 2;
    W           = abs(x_sample(end) - x_sample(1));             % Length of the combustor
    H           = 3.5*max(r_sample);                              % Diameter of the combustor
    plot_ratio  = 1.1;                                          % Ratio of the axes limit to the combustor dimension
    axes_W      = plot_ratio*W;                                 % axes x width
    axes_H      = 2*plot_ratio*H;                               % axes y width        
    x_min       = x_sample(1) - (axes_W-W)./2;                  % axes x min
    y_min       = -axes_H./2;                                   % axes y min
    
    % These are used for the legend
    xMax = max(get(gca,'xlim'));
    yMax = max(get(gca,'ylim'));
    
    h1=plot(gca,xMax,yMax,'-b','linewidth',3);
    h2=plot(gca,xMax,yMax,'-g','linewidth',3);
    h3=plot(gca,xMax,yMax,'-r','linewidth',3);
    h4=plot(gca,xMax,yMax,'-m','linewidth',3);   
    
    % Combustor outline
    for s=1:length(x_sample)-1

        x_plot1=1000*[x_sample(s),x_sample(s+1)];
        x_plot5=1000*[x_sample(s+1),x_sample(s+1)];
        y_plot1=1000*[r_sample(s),r_sample(s)];
        y_plot2=1000*[-r_sample(s),-r_sample(s)];
        y_plot3=1000*[r_sample(s),r_sample(s+1)];
        y_plot4=1000*[-r_sample(s),-r_sample(s+1)];

        plot(x_plot1, y_plot1,'-k','linewidth',WallLineWidth);
        plot(x_plot1, y_plot2,'-k','linewidth',WallLineWidth);

        if r_sample(s) ~=r_sample(s+1)
            plot(x_plot5, y_plot3,'-k','linewidth',WallLineWidth);
            plot(x_plot5, y_plot4,'-k','linewidth',WallLineWidth);
        end

    end

    % Combustor's sections  
    for s = 1:length(SectionIndex)

            if SectionIndex(s)==0          % just interface
                indexColor = 'k';
                indexLinestyle  = '-';
                indexLineWidth  = 1;
            elseif SectionIndex(s)==10        % with heat addition but perturbation
                indexColor = 'm';
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            elseif SectionIndex(s)==11         % with heat addition and heat perturbations
                indexColor = 'r';
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            end 
            if s == 1
                indexColor = 'b';   % inlet
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            elseif s == length(SectionIndex)
                indexColor = 'g';   % outlet
                indexLinestyle  = '-';
                indexLineWidth  = 3;
            end
        if TubeIndex(s) >= 1 && s > 1
            indexLinestyle  = 'none';
        end
        plot(1e3*[x_sample(s),x_sample(s)],-1e3*[-r_sample(s),r_sample(s)],...
            'linestyle',indexLinestyle,...
            'color',indexColor,'linewidth',indexLineWidth);
    end

    % For a varying cross section area
    diffTubeIndex = diff(TubeIndex);
    indexVarTubeIndex = find(diffTubeIndex~=0);
    if ~isempty(indexVarTubeIndex)
        for k = 1:length(indexVarTubeIndex)-1 % -1 to prevent black line at outlet
            s = indexVarTubeIndex(k)+1;

            plot(1e3*[x_sample(s),x_sample(s)],-1e3*[-r_sample(s),r_sample(s)],...
            'linestyle','-',...
            'color','k','linewidth',1);
        end
    end

    % Handle properties
    set(gca,'xlim',1000*[x_min, x_min+axes_W],'fontsize',16);
    set(gca,'box','on','linewidth',0.5,'gridlinestyle','-.');
    set(gca,'xgrid','on','ygrid','on');
    set(gca,'ylim',1000*[y_min, y_min+axes_H]);
    xlabel('x~ [mm]','Color','k','Interpreter','LaTex');
    ylabel('r~ [mm]','Color','k','Interpreter','LaTex');   

    % Legend properties
    
    legend1 = ['Steady HS']; 
    legend2 = ['Fluctuating HS'];
    
    if isempty(indexSHS) && isempty(indexFHS)
        hlegend = legend([h1,h2],'Inlet','Outlet');
    elseif isempty(indexSHS) && ~isempty(indexFHS)
        hlegend = legend([h1,h2,h3],'Inlet','Outlet',legend2);
    elseif ~isempty(indexSHS) && isempty(indexFHS)
        hlegend = legend([h1,h2,h4],'Inlet','Outlet',legend1);
%    else
%         hlegend = legend([h1,h2,h3,h4],'Inlet','Outlet',legend2,legend1);
    end          
    
    set(hlegend,'fontsize',14,'location','eastoutside');

    % Saving the figure
    if SAVE_FIGS
        saveas(fig1,'./Outputs/Initialisation/Geometry','fig')
    end
    if SAVE_PDFS
        save2pdf('./Outputs/Initialisation/Geometry',fig1,300)
    end
    if SMALL_PLOTS
        pos = [100 50 800 300];
        set(fig1, 'Position', pos)
        % Restore default font sizes
        set(gca, 'fontsize', get(gca,'DefaultAxesFontSize'));
        set(hlegend,'fontsize',get(gca,'DefaultLegendFontSize'));
    end
end

fprintf("Done.\n ");