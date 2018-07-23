% numerical measure
% voiced accuracy
reference = [[1250:1:3000] [3800:1:5600] [7800:1:11200] [13300:1:18700] [20300:1:22600] [25000:1:27400] []];
detected = [];
for i=1:length(B1)
    detected = [detected B1(i)*w_length:1:E1(i)*w_length];
end
count = 0;
for j=1:length(reference)
    for k=1:length(detected)
        if reference(j)==detected(k);
            count=count+1;
        end
    end
end
accuracy_voiced = count/length(reference);
 
% unvoiced accuracy
if isempty(E2)||isempty(B2)
    break;
else
   un_reference = [[5800:1:6800] [7100:1:7750] [11300:1:13300] [18800:1:19700] [28800:1:31200]]; 
   un_detected = [];
   for i=1:length(B2)
    un_detected = [un_detected B2(i)*w_length:1:E2(i)*w_length];
   end
count = 0;
for j=1:length(un_reference)
    for k=1:length(un_detected)
        if un_reference(j)==un_detected(k);
            count=count+1;
        end
    end
end
accuracy_un = count/length(un_reference);
end