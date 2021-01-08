% bisection algorithm
% capacity under the joint TPC+PAC constraints for a massive MIMO channel under
% favorable propagation
% inputs: 
% W : channel gains
% alpha : coefficient for grade of service
% PT : total transmit power constraint
% P1 : per-antenna power constraints


function BA_TPC_PAC(W,alpha,P1,PT)
m=length(W);
W_alpha=W.*alpha;

% desired uncertainty interval
bisectionerror=1;

% upper bound and lower bound for mu
xl=0;
xu=max(W_alpha);

while(bisectionerror>1e-8)
    
	% midpoint
	midpoint=(1/2)*(xl+xu);
    
	% redefine upper bound or lower bound
	for k=1:m
        a=max(0,(midpoint)^(-1)-(W_alpha(k)^(-1)));
        r(k)=alpha(k)*min(((P1(k))/(alpha(k))),a);
	end
         
	f_value=sum(r)-PT;
	if f_value<0
        xu=midpoint;
	elseif f_value>0
        xl=midpoint;
	end

	bisectionerror=abs(f_value);
	%if f_value=0, then bisectionerror=0, hence bisectionerror>1e-8 is false.
 
end
save('PAC_TPC_BA.mat')
end
