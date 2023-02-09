% written by M Rutledge, 2/9/2011

cdnow=cd;      cd(filenames{file}(1:end-4));
switch i  
    case 1
        f_to_1 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 2
        f_to_2 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 3
        f_to_3 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 4
        f_to_4 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 5
        f_to_5 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 6
        f_to_6 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 7
        f_to_7 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 8
        f_to_8 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 9
        f_to_9 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 10
        f_to_10 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 11
        f_to_11 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 12
        f_to_12 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 13
        f_to_13 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 14
        f_to_14 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 15
        f_to_15 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
    case 16
        f_to_16 =  fopen([filenames{file}(max(strfind(filenames{file},'\'))+1:end-4),'_ch',num2str(i),'.bin'],'a');
end

cd(cdnow);