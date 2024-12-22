function [yyyy] = getYear()

m = date;
tmp = strsplit(m,'-');

yyyy = tmp{3};

end