for i =0:9    
    filename = 'ImageFileName00';
    add = num2str(i);
    suffix = '.dcm';
    filename=strcat(filename,add);
    filename=strcat(filename , suffix);
    a(:,:,i+1)=dicomread(filename);
end

for i =10:99    
    filename = 'ImageFileName0';
    add = num2str(i);
    suffix = '.dcm';
    filename=strcat(filename,add);
    filename=strcat(filename , suffix);
    a(:,:,i+1)=dicomread(filename);
end

for i =100:122   
    filename = 'ImageFileName';
    add = num2str(i);
    suffix = '.dcm';
    filename=strcat(filename,add);
    filename=strcat(filename , suffix);
    a(:,:,i+1)=dicomread(filename);
end