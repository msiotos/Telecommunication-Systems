function [X] = bits_to_4_PAM(bit_seq,A)
    j=1;
    Z = [-3*A, -1*A, A , 3*A];
    X=zeros(1,length(bit_seq)/2);
    for i=1:2:length(bit_seq)
        if(bit_seq(i)==0 && bit_seq(i+1)==0)
            X(j) = Z(1);
        elseif(bit_seq(i)==0 && bit_seq(i+1)==1)
            X(j) = Z(2); 
        elseif(bit_seq(i)==1 && bit_seq(i+1)==1)
            X(j) = Z(3);
        elseif(bit_seq(i)==1 && bit_seq(i+1)==0)
            X(j) = Z(4);
        end
        j=j+1;
    end
end