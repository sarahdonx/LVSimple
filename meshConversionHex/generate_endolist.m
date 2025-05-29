clear all; close all; clc;

filename = 'endoNodes.txt';
fid = fopen(filename, 'r');
EndoNodes=[];
while ~feof(fid)
    tline=fgetl(fid);
    faces=sscanf(tline,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
    if ~isempty(faces)
       EndoNodes=[EndoNodes faces'];
    end
end
fclose(fid);

fid = fopen('endoList.point', 'w');
fprintf(fid, '%d\n', length(EndoNodes));
for i = 1 : length(EndoNodes)
    fprintf(fid, '%d\n', EndoNodes(i)); 
end
fclose(fid);