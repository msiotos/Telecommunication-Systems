function X = bits_to_2PAM(b)
X = 1:length(b);
for i=1:length(b)
if b(i)==0
X(i) = +1;
elseif b(i)==1
X(i) = -1;
else
disp('Error');
return;
end
end
end