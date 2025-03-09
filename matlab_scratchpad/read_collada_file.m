function [M] = read_collada_file(filename)


fid = fopen(filename);
s = fscanf(fid,'%c');
fclose(fid);

% search for the text segment that contains the materials
tag = '<library_materials>';
na = strfind(s,tag);
na = na + length(tag); % s(na) is the first character after the start tag

tag = '</library_materials>';
nb = strfind(s(na:end),tag);
nb = na + nb - 2; % s(nb) is the first character before the end tag

% string that contains the materials
sm = s(na:nb);
tag = '<material';
ind = strfind(sm,tag);
clear M
M.materials.id = [];
M.materials.name = [];
M.materials.parent = [];
M.materials = repmat(M.materials,length(ind),1);
for n = 1:length(ind)
    na = ind(n);                                    % first index of tag = character "<"
    nb = na - 1 + find(sm(na:end)=='>',1,'first');  % last index of tag = character ">"
    st = sm(na:nb);
    
    % get ID
    txt = 'id="';
    ka = strfind(st,txt) + length(txt);
    kb = ka - 2 + find(st(ka:end)=='"',1,'first');
    M.materials(n).id = st(ka:kb);
    
    % get name
    txt = 'name="';
    ka = strfind(st,txt) + length(txt);
    kb = ka - 2 + find(st(ka:end)=='"',1,'first');
    M.materials(n).name = st(ka:kb);
end

% search for the text segment that contains the materials
tag = '<library_geometries>';
na = strfind(s,tag);
na = na + length(tag); % s(na) is the first character after the start tag

tag = '</library_geometries>';
nb = strfind(s(na:end),tag);
nb = na + nb - 2; % s(nb) is the first character before the end tag

% string that contains the geometry
sg = s(na:nb);


tag = '<source';
ind = strfind(sg,tag);
for n = 1:length(ind)

    na = ind(n);
    nb = na - 1 + find(sg(na:end)=='>',1,'first');
    tmp = sg(na:nb);
    
    if contains(tmp,'positions">')
        na = nb + 1;
        nb = strfind(sg(na:end),'</source>');
        nb = na + nb(1) - 2;
        tmp = sg(na:nb); % contents between <source ...> and </source>
        
        na = strfind(tmp,'<float_array');
        na = na + find(tmp(na:end)=='>',1,'first');
        nb = strfind(tmp,'</float_array>') - 1;
        M.pos = sscanf(tmp(na:nb),'%f');
        M.pos = reshape(M.pos,3,length(M.pos)/3)';
    end

    if contains(tmp,'normals">')
        na = nb + 1;
        nb = strfind(sg(na:end),'</source>');
        nb = na + nb(1) - 2;
        tmp = sg(na:nb); % contents between <source ...> and </source>
        
        na = strfind(tmp,'<float_array');
        na = na + find(tmp(na:end)=='>',1,'first');
        nb = strfind(tmp,'</float_array>') - 1;
        M.nrl = sscanf(tmp(na:nb),'%f');
        M.nrl = reshape(M.nrl,3,length(M.nrl)/3)';
    end
end
   
tag = '<triangles';
ind = strfind(sg,tag);
for n = 1:length(ind)

    na = ind(n);
    nb = na - 1 + find(sg(na:end)=='>',1,'first');
    tmp = sg(na:nb);
    
    % get material ID
    txt = 'material="';
    ka = strfind(tmp,txt) + length(txt);
    kb = ka - 2 + find(tmp(ka:end)=='"',1,'first');
    this_material = tmp(ka:kb);
    
    na = nb + 1;
    nb = strfind(sg(na:end),'</triangles>');
    nb = na - 1 + nb(1);
    
    tmp = sg(na:nb); % contents between <triangles ...> and </triangles>
    
    ind_input = strfind(tmp,'<input');
    for k = 1:length(ind_input)
        na = ind_input(k);
        nb = na - 1 + find(tmp(na:end)=='>',1,'first');
        tmp_input = tmp(na:nb);
        
        % get semantic
        txt = 'semantic="';
        ka = strfind(tmp_input,txt) + length(txt);
        kb = ka - 2 + find(tmp_input(ka:end)=='"',1,'first');
        sem = tmp_input(ka:kb);
        
        % get offset
        txt = 'offset="';
        ka = strfind(tmp_input,txt) + length(txt);
        kb = ka - 2 + find(tmp_input(ka:end)=='"',1,'first');
        off = tmp_input(ka:kb);
        
        if strcmp(sem,'VERTEX')
            offset_vertex = str2double(off);
        elseif strcmp(sem,'NORMAL')
            offset_normal = str2double(off);
        end
    end
    
    na = strfind(tmp,'<p');
    na = na + find(tmp(na:end)=='>',1,'first');
    nb = strfind(tmp,'</p>') - 1;
    parent = sscanf(tmp(na:nb),'%d'); % index starts with 0, change that to 1
    
    for k = 1:length(M.materials)
        if strcmp(M.materials(k).id,this_material)
            M.materials(k).parent = parent + 1;
            M.materials(k).offset_vertex = offset_vertex;
            M.materials(k).offset_normal = offset_normal;
            M.materials(k).num_inputs = length(ind_input);
            break
        end
    end
end

c = 0;
for k = 1:length(M.materials)
    if ~isempty(M.materials(k).parent)
        c = c + length(M.materials(k).parent) / (M.materials(k).num_inputs*3);
    end
end
M.v0 = zeros(c,3);
M.v1 = zeros(c,3);
M.v2 = zeros(c,3);
M.id = zeros(c,1);
c = 0;
for k = 1:length(M.materials)
    ind1 = M.materials(k).parent(1+0*M.materials(k).num_inputs+M.materials(k).offset_vertex:3*M.materials(k).num_inputs:end); % vertex 1
    ind2 = M.materials(k).parent(1+1*M.materials(k).num_inputs+M.materials(k).offset_vertex:3*M.materials(k).num_inputs:end); % vertex 2
    ind3 = M.materials(k).parent(1+2*M.materials(k).num_inputs+M.materials(k).offset_vertex:3*M.materials(k).num_inputs:end); % vertex 3
    indn = M.materials(k).parent(1+M.materials(k).offset_normal:3*M.materials(k).num_inputs:end); % normal
    
    % check sequence of vertices
    for n = 1:length(ind1)
        %fprintf('%d-%d',k,n)
        v1 = M.pos(ind1(n),:);
        v2 = M.pos(ind2(n),:);
        v3 = M.pos(ind3(n),:);
        
        my_nrl = cross(v2-v1,v3-v1);
        my_nrl = my_nrl ./ repmat(sqrt(sum(my_nrl.^2)),1,3);
        collada_nrl = M.nrl(indn(n),:);
        
        if my_nrl * collada_nrl' < 0
            %warning('Swapping vertices of triangle %d of material %d to correct normal.',n,k)
            
            tmp = M.pos(ind2(n),:);
            M.pos(ind2(n),:) = M.pos(ind3(n),:);
            M.pos(ind3(n),:) = tmp;
            
            tmp = v2;
            v2 = v3;
            v3 = tmp;
        end
        
        c = c + 1;
        M.v0(c,:) = v1;
        M.v1(c,:) = v2;
        M.v2(c,:) = v3;
        M.id(c) = k;
    end
end

% recalculate normals (for indexing)
M.nrl = cross(M.v1-M.v0,M.v2-M.v0);
M.nrl = M.nrl ./repmat(sqrt(sum(M.nrl.^2,2)),1,3);