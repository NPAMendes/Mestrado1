function f = attack_ay(t,n)
%FAULT Summary of this function goes here
%   Detailed explanation goes here

f = zeros(n,(length(t)+1));

for i = 1:length(t)
    if (t(i) < 150)
        f(1:n,i) = 0;
    else
        if (t(i)>= 150 && t(i) < 550)
            f(1:n,i) = 1.5;
        else
            if (t(i)>= 550 && t(i) < 650)
                f(1:n,i) = 0;
            else
                if (t(i)>= 650 && t(i) < 1300)
                    f(1:n,i) = 1*sin(pi*(t(i)-1300)/130);
                else
                    if (t(i)>= 1300 && t(i) < 1400)
                        f(1:n,i) = 0;
                    else
                        if (t(i)>= 1400 && t(i) < 1900)
                            f(1:n,i) = t(i)/380-3;
                        else
                            f(1:n,i) = 0;
                        end
                    end 
                end 
            end 
        end 
    end
end

end