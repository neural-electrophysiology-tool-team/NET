function [y]=choose_real(x)

n=length(x);
i=1;

while ((isreal(x(i))==0 || x(i)<=0) && i<n)
    i=i+1;
end;

y=x(i);

