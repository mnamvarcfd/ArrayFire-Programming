function arrayToFile(variable_name, fileName, arrayType, filePath)

% Typical arrayType for matlab are: double, single, int64, ...
% "help format" to get more.

% Usage example : 
% arrayToFile(rho, 'rho_final', 'double', '/Users/jesusgp/Desktop/VM/Github/AF_LBM_JG/MAIN_AF_LBM_POISEUILLE2D/workingDirectory/');
% If file is in the same directory use something like:
% arrayToFile(rho, 'rho_final', 'double', './');


arrayName = '.array';
fileName_array = strcat(filePath,fileName,arrayName);

% Write ascii dimensions
fileID = fopen(fileName_array, 'w+', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', fileName_size); end

sizeVariable = size(variable_name);

if numel(sizeVariable)>4
    error('ArrayFire does not support array with dim greater than 4.') ;
end

try
    fprintf(fileID,'%s ',num2str(sizeVariable(1)));
catch
    error('Array is empty, not writing any data to file');
end    
try
    fprintf(fileID,'%s ',num2str(sizeVariable(2)));
catch
    fprintf(fileID,'%s ',num2str(1));
end    
try
    fprintf(fileID,'%s ',num2str(sizeVariable(3)));
catch
    fprintf(fileID,'%s ',num2str(1));
end    
try
    fprintf(fileID,'%s\n',num2str(sizeVariable(4)));
catch
    fprintf(fileID,'%s\n',num2str(1));
end    

% Write the data content of the array
fwrite(fileID, variable_name, arrayType);

% closing file
fclose(fileID);
