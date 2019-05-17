%You shouldn't need to touch this

function changeFile(newHeight,newPhi,newTheta,month,day,year,hour,minute,second,fileFolder)
fileLocation = strcat(fileFolder,'\NameRef.txt');
% Read txt into cell A
fid = fopen(fileLocation,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Change cell A
A{11} = sprintf('  h1 = %d', newHeight);
A{12} = sprintf('  phi1 = %d', newPhi);
A{13} = sprintf('  thet1 = %d', newTheta);
A{24} = sprintf('  mn = %d', month);
A{25} = sprintf('  ida = %d', day);
A{26} = sprintf('  iyr = %d', year);
A{27} = sprintf('  ihro = %d', hour);
A{28} = sprintf('  mino = %d', minute);
A{29} = sprintf('  seco = %d', second);


% Write cell A into txt
fid = fopen('My Test\NameRef.txt', 'w');
for i = 1:numel(A)
    if A{i} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose('all');
end
