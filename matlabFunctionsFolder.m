function folder=matlabFunctionsFolder()
[pathstr,name]=fileparts(which('identifierOfMainDir4NSKToolBox.m'));
p=find(pathstr(1:end-1)==filesep);
folder=pathstr(1:p(end)-1);