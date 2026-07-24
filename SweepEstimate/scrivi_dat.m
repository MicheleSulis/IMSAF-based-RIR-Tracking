function  scrivi_dat(filename,x)

% scrivi_dat(filename,x)
% scrive x in un file dat 'filename.dat'

fp = fopen(filename, 'wb');

%fwrite(fp, length(x), 'uint32');
fwrite(fp, x, 'float64');

%fwrite(fp, length(y), 'uint32');
%fwrite(fp, y, 'double');

fclose(fp);