function allFiles=myNSKFunctionList(NSKRoot)
%functions returns a structure of all the functions and mat files in the sub dirs of the NSK toolbox
allFiles=[];
if nargin==0
    NSKRoot=fileparts(which('identifierOfMainDir4NSKToolBox'));
end
D=dir(NSKRoot);
D = D(find(~cellfun(@(x) x(1)=='.',{D(:).name}))); %remove ./.. entries

for i=1:length(D)
    [~,origName,ext]=fileparts(D(i).name);
    file=origName;
    file(isspace(file))='_';
    file(file=='.')='_';
    file(file=='-')='_';
    file(file=='~')='_';
    file(file=='(')='_';
    file(file==')')='_';
    file(file=='@')=[];
    file(file=='+')=[];
    file(file=='´')=[];
    file(file(1)=='_')='x';
    if isstrprop(file(1), 'digit')
        disp([file ' changed to X' file(2:end)] );
        file(1)='X';
    end
    
    if D(i).isdir
        allFiles.(file)=myNSKFunctionList([NSKRoot filesep origName]);
    else
        allFiles.(file)=ext(2:end);
    end
end