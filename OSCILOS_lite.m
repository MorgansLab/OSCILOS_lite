%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OSCILOS_lite
%
% This program is the 'lite' version of OSCILOS_long - without a 
% Graphical User Interface. It runs much faster than OSCILOS_long and can
% be used to carry out parametric investigations.
%
% Please update the input files before running this script. More
% information about the required inputs is included in the User Guide.
%
% Last update : 08/01/2021
%
% Authors: R. Gaudron, A. MacLaren, and previous contributors (see full
% list on the OSCILOS website).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALISATION
%
% This subroutine loads the 'Init.txt' file and initialises OSCILOS_lite

addpath('./SubFunctions')
run Init_subfc

%% GEOMETRY
%
% This subroutine loads the geometry from the file 'Geometry.txt', plots it
% and saves it in the output folder.

run Geometry_subfc

%% MEAN FLOW
%
% This subroutine loads the mean flow parameters at the inlet of the geometry 
% from the file 'Mean_flow.txt'. The mean conservation equations are then 
% solved and the mean flow variables are obtained throughout the geometry.
% The mean velocity and mean temperature across the geometry are then
% plotted and saved in the output folder.

run Mean_flow_subfc

%% FLAME MODEL
%
% This subroutine loads the flame model from the file 'Flame.txt'. The gain
% and phase of the corresponding Flame Model are then plotted and saved in
% the output folder. If heat perturbations are not present inside the
% domain, this step has no impact on the final result.

run Flame_subfc

%% BOUNDARY CONDITIONS
%
% This subroutine loads the acoustic boundary conditions at the inlet and 
% outlet of the geometry from files 'Inlet.txt' and 'Outlet.txt'
% respectively. The acoustic reflection coefficient at the inlet and outlet
% are then plotted and saved in the output folder.

run BC_subfc

%% FREQUENCY-DOMAIN SOLVER
%
% This subroutine loads the scan range and number of initialisations - fed to
% the frequency-domain solver - from the file 'Scan_range.txt'. The solver
% then finds the eigenvalues of the system, and the corresponding
% eigenmodes. The eigenvalues are saved as an output file called
% 'Eigenvalues.txt'. A map of the eigenvalues in the range of interest is
% then plotted and saved in the output folder. A number of eigenmodes,
% represented as the modulus of the acoustic pressure and velocity, are
% also plotted and saved in the output folder.

run Solver_subfc

%% MESSAGE BOX

time=toc;
time_final=num2str(sprintf('%.2f',time));
text=['Calculation completed in ',time_final,' seconds'];
cltext = " Close Plots ";
f = questdlg(text,'OSCILOS_lite','OK',	cltext, 'OK');
if f==cltext
    close all
end
fprintf("\n OSCILOS_lite ran successfully in %s seconds.\n",time_final)
