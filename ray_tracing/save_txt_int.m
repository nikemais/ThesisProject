function save_txt_int(filename,A)
% save_txt_int(filename,A)

fid = fopen(filename,'w');
fprintf(fid,'%% %d %d\n',size(A,1),size(A,2));
for m = 1:size(A,1)
    for n = 1:size(A,2)
        fprintf(fid,'%24d  ',A(m,n));
    end
    fprintf(fid,'\n');
end
fclose(fid);
    
