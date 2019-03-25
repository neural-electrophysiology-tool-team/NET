function str = fromUpperToCamelCase(upperscore_compound)
% converts upperscore_compound to CamelCase
%
% Not always exactly invertible
%
% Examples:
%   fromUpperToCamelCase('ONE')            -->  'one'
%   fromUpperToCamelCase('ONE_TWO_THREE')  -->  'oneTwoThree'
%   fromUpperToCamelCase('#$ONE_TWO_THREE') --> 'oneTwoThree'
%   fromUpperToCamelCase('One_Two_Three')  --> !error! lower case only mixes with alphanumerics
%   fromUpperToCamelCase('5_TWO_THREE')    --> !error! cannot start with a digit
assert(isempty(regexp(upperscore_compound, '\s', 'once')), ... 
  'white space is not allowed' )
assert(~ismember(upperscore_compound(1), '0':'9'), 'string cannot begin with a digit')
str = lower(upperscore_compound);
str = regexprep(str, '(^|[_\W]+)([a-zA-Z])', '${upper($2)}');
str = [lower(str(1)) str(2:end)];
end