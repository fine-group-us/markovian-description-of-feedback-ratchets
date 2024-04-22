function theta_1x = Protocol(deltax)
% This fucnction computes the protocol with a measure uncertanity \Delta x
L=2;
a=1;
xdis=linspace(0,L,1000); % x is shifted to the interval (0,L=2)
theta_1x=zeros(1, length(xdis));
p=0;
for i=1:length(xdis)
    x=xdis(i);
    l_u= mod(x+deltax,L); % upper limit of the measure interval
    l_l=mod(x-deltax,L);  % lower limit of the measure interval
    if (deltax==0)
        p=(x<=a);
        theta_1x(i)=p;
    else
        if (l_l<=l_u)
            if (l_u<=a)
                p=(l_u-l_l)/(2*deltax);
            else %l_u >a
                if (l_l<=a)
                    p=(a-l_l)/(2*deltax);
                else
                    p=0;
                end
            end
        else %l_l>l_u
            if (l_l<=a)
                p=(a-l_l+l_u)/(2*deltax);
            else %l_l>a
                if (l_u<=a)
                    p=l_u/(2*deltax);
                else
                    p= a/(2*deltax);
                end
            end
        end
        theta_1x(i)=p;
    end
end