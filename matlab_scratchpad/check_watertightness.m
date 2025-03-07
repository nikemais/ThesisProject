function [result] = check_watertightness(v1,v2,v3,nrl)
% [result] = check_watertightness(v1,v2,v3,nrl)

result = true;

N = size(v1,1);

for n = 1:N
    %----------------------------------------------------------------------
    % check v1/v2 edge
    %----------------------------------------------------------------------
    sa =   sum(all(v1(n,:)==v1,2) & all(v2(n,:)==v2,2)) ... 
         + sum(all(v1(n,:)==v2,2) & all(v2(n,:)==v3,2)) ... 
         + sum(all(v1(n,:)==v3,2) & all(v2(n,:)==v1,2));
    
    sb =   sum(all(v1(n,:)==v1,2) & all(v2(n,:)==v3,2)) ... 
         + sum(all(v1(n,:)==v3,2) & all(v2(n,:)==v2,2)) ... 
         + sum(all(v1(n,:)==v2,2) & all(v2(n,:)==v1,2));
     
    % this checks if an edge appears in 2 triangles, and if the normals
    % point into the same direction
    if sa ~= 1 || sb ~= 1
        s12 = false; % bad
    else
        s12 = true;
    end
    
    %----------------------------------------------------------------------
    % check v2/v3 edge
    %----------------------------------------------------------------------
    sa =   sum(all(v2(n,:)==v1,2) & all(v3(n,:)==v2,2)) ... 
         + sum(all(v2(n,:)==v2,2) & all(v3(n,:)==v3,2)) ... 
         + sum(all(v2(n,:)==v3,2) & all(v3(n,:)==v1,2));
    
    sb =   sum(all(v2(n,:)==v1,2) & all(v3(n,:)==v3,2)) ... 
         + sum(all(v2(n,:)==v3,2) & all(v3(n,:)==v2,2)) ... 
         + sum(all(v2(n,:)==v2,2) & all(v3(n,:)==v1,2));
     
    % this checks if an edge appears in 2 triangles, and if the normals
    % point into the same direction
    if sa ~= 1 || sb ~= 1
        s23 = false; % bad
    else
        s23 = true;
    end
    
    %----------------------------------------------------------------------
    % check v3/v1 edge
    %----------------------------------------------------------------------
    sa =   sum(all(v3(n,:)==v1,2) & all(v1(n,:)==v2,2)) ... 
         + sum(all(v3(n,:)==v2,2) & all(v1(n,:)==v3,2)) ... 
         + sum(all(v3(n,:)==v3,2) & all(v1(n,:)==v1,2));
    
    sb =   sum(all(v3(n,:)==v1,2) & all(v1(n,:)==v3,2)) ... 
         + sum(all(v3(n,:)==v3,2) & all(v1(n,:)==v2,2)) ... 
         + sum(all(v3(n,:)==v2,2) & all(v1(n,:)==v1,2));
     
    % this checks if an edge appears in 2 triangles, and if the normals
    % point into the same direction
    if sa ~= 1 || sb ~= 1
        s31 = false; % bad
    else
        s31 = true;
    end
        
    if s12 == false || s23 == false || s31 == false
        warning('Edge %d is not good.',n)
        result = false;
    end
end
fprintf("ok")
%----------------------------------------------------------------------
% check normals
%----------------------------------------------------------------------
n = find(sum(cross(v2-v1,v3-v1).*nrl,2)<0);
if ~isempty(n)
    warning('Normal %d is not good.',n)
end
    
    