close all
clear variables
clc


% pth = 'C:/Users/nike/Documents/ThesisProject/matlab_scratchpad/symmetry_test/';
% fn = 'test_2';
pth = 'C:/Users/nike/Documents/ThesisProject/matlab_scratchpad/MRO_lowfidelity/';
fn = 'MRO_lowfidelity';
M = read_collada_file([pth fn '.dae']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert units from mm to m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.pos = M.pos/1000;
M.v0 = M.v0/1000;
M.v1 = M.v1/1000;
M.v2 = M.v2/1000;

% note: no rotation needed for GRACE-FO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's check if triangulation is watertight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result =  check_watertightness(M.v0,M.v1,M.v2,M.nrl);
if result ~= true
    error('Mesh is not watertight!')
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a plot to check if assign materials to triangles worked
% (this includes exporting as Collada file and importing to Matlab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cols = [166,206,227
31,120,180
178,223,138
51,160,44
251,154,153
227,26,28
253,191,111
255,127,0
202,178,214
106,61,154
255,255,153
177,89,40
90,90,90
190,190,190]/255;

for k = 1:length(M.materials)
    disp(k)

    x = [M.v0(M.id==k,1) M.v1(M.id==k,1) M.v2(M.id==k,1)];
    y = [M.v0(M.id==k,2) M.v1(M.id==k,2) M.v2(M.id==k,2)];
    z = [M.v0(M.id==k,3) M.v1(M.id==k,3) M.v2(M.id==k,3)];
    for n = 1:size(x,1)
        fill3(x(n,:),y(n,:),z(n,:),cols(mod(k-1,size(cols,1))+1,:))
        hold on
    end
    
    % show normals
    for n = 1:size(x,1)
        xm = mean(x(n,:));
        ym = mean(y(n,:));
        zm = mean(z(n,:));
        
        v1 = [x(n,1),y(n,1),z(n,1)];
        v2 = [x(n,2),y(n,2),z(n,2)];
        v3 = [x(n,3),y(n,3),z(n,3)];
        nrl = cross(v2-v1,v3-v1);
        nrl = nrl/sqrt(sum(nrl.^2));
        len = 0.2;
        plot3(xm + [0, len*nrl(1)],...
              ym + [0, len*nrl(2)],...
              zm + [0, len*nrl(3)],'k','LineWidth',2)
    end
end
hold off
axis equal
xlabel('x')
zlabel('z')
ylabel('y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print list materials with names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(M.materials)
    fprintf('%d: %s\n',k,M.materials(k).name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define thermo-optical surface properties of materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1: SolarPanelSun: Si_Glass_Solar_Array (from product spec doc of Swarm)
M.materials(1).optical_absorption           = 0.90; % use Swarm
M.materials(1).optical_diffuse_reflection   = 0.07;
M.materials(1).optical_specular_reflection  = 0.03;
M.materials(1).infrared_absorption          = 0.81;
M.materials(1).infrared_diffuse_reflection  = 0.16;
M.materials(1).infrared_specular_reflection = 0.03;

% % 2: SolarPanelBM : Si_Glass_Solar_Array (from product spec doc of Swarm)
% M.materials(2).optical_absorption           = 0.90; % use Swarm
% M.materials(2).optical_diffuse_reflection   = 0.07;
% M.materials(2).optical_specular_reflection  = 0.03;
% M.materials(2).infrared_absorption          = 0.81;
% M.materials(2).infrared_diffuse_reflection  = 0.16;
% M.materials(2).infrared_specular_reflection = 0.03;
% 
% % 3: GPSAntenna: white paint, silicate --> from Space Systems Engineering book
% M.materials(3).optical_absorption           = 0.14;
% M.materials(3).optical_diffuse_reflection   = 0.80; % my own guess
% M.materials(3).optical_specular_reflection  = 0.06; % my own guess
% M.materials(3).infrared_absorption          = 0.90;
% M.materials(3).infrared_diffuse_reflection  = 0.05; % my own guess
% M.materials(3).infrared_specular_reflection = 0.05; % my own guess
% 
% % 4: SBandAntenna: white paint, silicate --> from Space Systems Engineering book
% M.materials(4).optical_absorption           = 0.14;
% M.materials(4).optical_diffuse_reflection   = 0.80; % my own guess
% M.materials(4).optical_specular_reflection  = 0.06; % my own guess
% M.materials(4).infrared_absorption          = 0.90;
% M.materials(4).infrared_diffuse_reflection  = 0.05; % my own guess
% M.materials(4).infrared_specular_reflection = 0.05; % my own guess
% 
% % 5: CESSWhite: white paint, silicate --> from Space Systems Engineering book
% M.materials(5).optical_absorption           = 0.14;
% M.materials(5).optical_diffuse_reflection   = 0.80; % my own guess
% M.materials(5).optical_specular_reflection  = 0.06; % my own guess
% M.materials(5).infrared_absorption          = 0.90;
% M.materials(5).infrared_diffuse_reflection  = 0.05; % my own guess
% M.materials(5).infrared_specular_reflection = 0.05; % my own guess
% 
% % 6: CESSBlack: black paint, polyurethane --> from Space Systems Engineering book
% M.materials(6).optical_absorption           = 0.90;
% M.materials(6).optical_diffuse_reflection   = 0.05; % my own guess
% M.materials(6).optical_specular_reflection  = 0.05; % my own guess
% M.materials(6).infrared_absorption          = 0.90;
% M.materials(6).infrared_diffuse_reflection  = 0.05; % my own guess
% M.materials(6).infrared_specular_reflection = 0.05; % my own guess


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check if coefficient add up to one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(M.materials)
    if M.materials(n).optical_absorption + M.materials(n).optical_diffuse_reflection + M.materials(n).optical_specular_reflection ~= 1
        error('optical')
    end
    if M.materials(n).infrared_absorption + M.materials(n).infrared_diffuse_reflection + M.materials(n).infrared_specular_reflection ~= 1
        error('infrared')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Convert into my own format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v0 = M.v0;
v1 = M.v1;
v2 = M.v2;
id = M.id;
pa_optical = zeros(size(id));
pd_optical = zeros(size(id));
ps_optical = zeros(size(id));
pa_infrared = zeros(size(id));
pd_infrared = zeros(size(id));
ps_infrared = zeros(size(id));
for k = 1:length(M.materials)
    pa_optical(id==k) = M.materials(k).optical_absorption;
    pd_optical(id==k) = M.materials(k).optical_diffuse_reflection;
    ps_optical(id==k) = M.materials(k).optical_specular_reflection;
    pa_infrared(id==k) = M.materials(k).infrared_absorption;
    pd_infrared(id==k) = M.materials(k).infrared_diffuse_reflection;
    ps_infrared(id==k) = M.materials(k).infrared_specular_reflection;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final check plot before saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrl = cross(v1-v0,v2-v0);
nrl = nrl./repmat(sqrt(sum(nrl.^2,2)),1,3);
for k = 1:size(v0,1)
    % show faces
    fill3([v0(k,1),v1(k,1),v2(k,1)],...
          [v0(k,2),v1(k,2),v2(k,2)],...
          [v0(k,3),v1(k,3),v2(k,3)],cols(mod(id(k)-1,size(cols,1))+1,:))
    hold on
    
    % show normals
    xm = mean([v0(k,1),v1(k,1),v2(k,1)]);
    ym = mean([v0(k,2),v1(k,2),v2(k,2)]);
    zm = mean([v0(k,3),v1(k,3),v2(k,3)]);
    len = 0.2;
    plot3(xm + [0, len*nrl(k,1)],...
          ym + [0, len*nrl(k,2)],...
          zm + [0, len*nrl(k,3)],'k','LineWidth',2)
end
hold off
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_txt([pth fn '_nrl.txt'],nrl')
save_txt([pth fn '_v0.txt'],v0')
save_txt([pth fn '_v1.txt'],v1')
save_txt([pth fn '_v2.txt'],v2')
save_txt([pth fn '_pa_optical.txt'],pa_optical')
save_txt([pth fn '_pd_optical.txt'],pd_optical')
save_txt([pth fn '_ps_optical.txt'],ps_optical')
save_txt([pth fn '_pa_infrared.txt'],pa_infrared')
save_txt([pth fn '_pd_infrared.txt'],pd_infrared')
save_txt([pth fn '_ps_infrared.txt'],ps_infrared')
save_txt_int([pth fn '_id.txt'],id')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




