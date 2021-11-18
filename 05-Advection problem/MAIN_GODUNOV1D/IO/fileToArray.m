function variable_name = fileToArray(fileName, arrayType, filePath)

% Typical arrayType for matlab are: double, single, int64, ...
% "help format" to get more.
% Usage example: 
% rho = fileToArray('rho_40000', 'double', '/Users/jesusgp/Desktop/VM/Github/AF_LBM_JG/MAIN_AF_LBM_POISEUILLE2D/workingDirectory/');
% If file is in the same directory use something like:
% rho = fileToArray('rho_40000', 'double', './'); 


arrayName = '.array';
fileName_array = strcat(filePath,fileName,arrayName);

% Read ascii dimensions
fileID = fopen(fileName_array, 'r', 'ieee-le');
if fileID == -1, error('Cannot open file: %s', fileName_array); end
size_mat_string = fgets(fileID);
spaceLocations = find(size_mat_string == ' ');
size_mat(1) = str2double(size_mat_string(1:spaceLocations(1)-1));
size_mat(2) = str2double(size_mat_string(spaceLocations(1):spaceLocations(2)-1));
size_mat(3) = str2double(size_mat_string(spaceLocations(2):spaceLocations(3)-1));
size_mat(4) = str2double(size_mat_string(spaceLocations(3):end));

% Read the data content of the array and reshape
variable_name = fread(fileID, Inf, arrayType);
fclose(fileID);

variable_name = reshape(variable_name,size_mat);

