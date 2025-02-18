function save_txt(filename, array)
    % SAVE_TXT Saves an array to a text file
    %   save_txt(filename, array) saves the given array to the specified filename.
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Check if the file opened successfully
    if fid == -1
        error('Could not open file for writing.');
    end
    
    % Determine the number of columns
    [~, cols] = size(array);
    
    % Format string for writing rows
    formatSpec = [repmat('%g ', 1, cols-1), '%g\n']; 
    
    % Write the array to the file
    fprintf(fid, formatSpec, array');
    
    % Close the file
    fclose(fid);
    
    disp(['Array saved to ', filename]);
end