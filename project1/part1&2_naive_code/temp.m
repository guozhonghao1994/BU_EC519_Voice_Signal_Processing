
cankao = pitch(:,2);
cankao = cankao';
for i=1:length(key)
    if isnan(cankao(i))
        key(i)=NaN;
    end
end
for i=1:length(key)
    key(i)=key(i)-10;
end
sum = 0;
for i=1:length(key)
    if isnan(key(i))||isnan(cankao(i))
        continue
    else
        sum = sum+(key(i)-cankao(i))^2;
    end
end
result = sqrt(sum/length(key));







stem(cankao)
hold on
stem(key,'r')