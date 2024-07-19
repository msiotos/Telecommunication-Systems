function [est_bit] = PAM_4_to_bits(X,A)
    i=1;
    for i=1:length(X)

        if(X(i)== -3*A)          
          est_bit(i)=0;
          est_bit(i+1)=0;

        elseif(X(i)== -A)
          est_bit(i)=0;
          est_bit(i+1)=1; 

        elseif(X(i)== A)
          est_bit(i)=1;
          est_bit(i+1)=1; 

        elseif(X(i)== 3*A)
          est_bit(i)=1;
          est_bit(i+1)=0; 
        end
        i=i+2;
    end
end
