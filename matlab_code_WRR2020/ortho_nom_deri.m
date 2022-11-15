

function [output] = ortho_nom_deri(order,x)

 load('orthopolycoeff.mat');


 

% switch order
%     case 1
%         output = sqrt(3)*(2*x-1);
%     case 2
%         output = 6*sqrt(5)*(x.^2-x+1/6);
%     case 3
%         output1 = 20*sqrt(7)*(x.^3-1.5*x.^2+0.6*x-1/20);
%     otherwise
%         warning('order must be integer between 1 and 3');
%         output = NaN;
 
 
switch order
    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
       
        order_ma = order+1; 
        nsamp = length(x);
        y = zeros(nsamp,order_ma);
        for i = 1:order
             y(:,i) = ORTHOPOLYCOEFLEG(order_ma,i)*(order_ma-i)*x.^(order-i);
        end
     
        output = sum(y,2);
        
    otherwise
         warning('order must be an integer between 1 and 15');
         output = NaN;
        
end


end

    