function save_txt_int(filename, array)
    % SAVE_TXT_INT Saves an integer array to a text file
    %   save_txt_int(filename, array) saves the given integer array to the specified filename.
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Check if the file opened successfully
    if fid == -1
        error('Could not open file for writing.');
    end
    
    % Determine the number of columns
    [~, cols] = size(array);
    
    % Format string for writing rows (integer format)
    formatSpec = [repmat('%d ', 1, cols-1), '%d\n']; 
    
    % Write the array to the file
    fprintf(fid, formatSpec, array');
    
    % Close the file
    fclose(fid);
    
    disp(['Integer array saved to ', filename]);
end