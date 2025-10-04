function f = attack_au(t,n)
%FAULT Summary of this function goes here
%   Detailed explanation goes here

f = zeros(n,length(t));

for i = 1:length(t)
    if (t(i) < 150)
        f(1:n,i) = 0;
    else
        if (t(i)>= 150 && t(i) < 500)
            f(1:n,i) = 2*1.5*t(i)/750;
        else
            if (t(i)>= 500 && t(i) < 750)
                f(1:n,i) = 0;
            else
                if (t(i)>= 750 && t(i) < 1200)
                    f(1:n,i) = 0.75;
                else
                    if (t(i)>= 1200 && t(i) < 1350)
                        f(1:n,i) = 0;
                    else
                        f(1:n,i) = 1*sin(5*pi*(t(i)-1350)/1000);
                    end 
                end 
            end 
        end 
    end
end

end