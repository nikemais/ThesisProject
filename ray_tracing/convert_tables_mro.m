
clear variables
close all
clc


%--------------------------------------------------------------------------
% load ray tracing model
%--------------------------------------------------------------------------

vis.c = load('test/table_c_op.txt'); % Plug your table 
ir.c = load('test/table_c_ir.txt'); % Plug your table 

% the angles are the same for both optical and infrared
attack = load('test/table_attack_angles.txt');
sideslip = load('test/table_sideslip_angles.txt')+90;


%--------------------------------------------------------------------------
% convert from grid frame to spacecraft frame
%--------------------------------------------------------------------------

N = length(attack);

R_GRD2SC = [ cosd(attack).*cosd(sideslip) ...
            -cosd(attack).*sind(sideslip) ...
            -sind(attack) ...
             sind(sideslip) ...
             cosd(sideslip) ...
             zeros(N,1) ...
             sind(attack).*cosd(sideslip) ...
            -sind(attack).*sind(sideslip) ...
             cosd(attack)];
         
vis.c = rot_mat_times_vec_vec(R_GRD2SC,vis.c);
ir.c = rot_mat_times_vec_vec(R_GRD2SC,ir.c);


%--------------------------------------------------------------------------
% save lookup tables
%--------------------------------------------------------------------------

table = [vis.c ir.c]; % This is the total RP force coefficient that you need

fid = fopen('table_v2.txt','w');

fprintf(fid,'  %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e \n',table');
fclose(fid);



% table = [attack sideslip vis.c];
% 
% fid = fopen('../GRACE_radiation_pressure_lookup_table_visible.txt','w');
% fprintf(fid,'# Lookup table for radiation pressure force coefficients for the GRACE satellites\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# The table was produced by TU Delft for the Swarm DISC TOLEOS project. The table\n');
% fprintf(fid,'# and information on the  TOLEOS project are available on the project''s website:\n');
% fprintf(fid,'# http://thermosphere.tudelft.nl/toleos.html\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# The force coefficients are provided in the satellite body frame. In the nominal\n');
% fprintf(fid,'# flight situation, the x-axis is pointing towards the other satellite, the z-axis\n');
% fprintf(fid,'# is pointing towards Earth, and the y-axis completes the right-hand system. The\n');
% fprintf(fid,'# coefficients are provided for the visible (vis) wavelength.\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# Lookup columns (first number), data columns (second number), and table rows ():\n');
% fprintf(fid,'2 6 %d\n',size(table,1));
% fprintf(fid,'# Attack angle   Sideslip angle Cx_vis         Cy_vis         Cz_vis\n');
% fprintf(fid,'# (°)            (°)            (m2)           (m2)           (m2)\n');
% fprintf(fid,'  %14.7e %14.7e %14.7e %14.7e %14.7e\n',table');
% fclose(fid);
% 
% 
% 
% 
% table = [attack sideslip ir.c];
% 
% fid = fopen('../GRACE_radiation_pressure_lookup_table_infrared.txt','w');
% fprintf(fid,'# Lookup table for radiation pressure force coefficients for the GRACE satellites\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# The table was produced by TU Delft for the Swarm DISC TOLEOS project. The table\n');
% fprintf(fid,'# and information on the  TOLEOS project are available on the project''s website:\n');
% fprintf(fid,'# http://thermosphere.tudelft.nl/toleos.html\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# The force coefficients are provided in the satellite body frame. In the nominal\n');
% fprintf(fid,'# flight situation, the x-axis is pointing towards the other satellite, the z-axis\n');
% fprintf(fid,'# is pointing towards Earth, and the y-axis completes the right-hand system. The\n');
% fprintf(fid,'# coefficients are provided for the infrared (ir) wavelength.\n');
% fprintf(fid,'#\n');
% fprintf(fid,'# Lookup columns (first number), data columns (second number), and table rows ():\n');
% fprintf(fid,'2 6 %d\n',size(table,1));
% fprintf(fid,'# Attack angle   Sideslip angle Cx_ir          Cy_ir          Cz_ir\n');
% fprintf(fid,'# (°)            (°)            (m2)           (m2)           (m2)\n');
% fprintf(fid,'  %14.7e %14.7e %14.7e %14.7e %14.7e\n',table');
% fclose(fid);




