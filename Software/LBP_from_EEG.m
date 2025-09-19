nb_seizures = [0,4,0,4,6,2,0,3,0,0,2,0,2,0,0,2];


for pat=[2,4,5,6,8,11,13,16]
for i=1:nb_seizures(pat)
load("no_backup\Pat"+string(pat)+"\Sz"+string(i)+".mat");
size(EEG)

l = 6; % 64 different LBP's if l is 6
BIG_LBP = zeros(size(EEG));
for col = 1:length(EEG(1,:))
    EEG_column = EEG(:,col);
    LBP = dec2bin(0,l);
    for k = 2:length(EEG_column)
        LBP = circshift(LBP,-1);
        if (EEG_column(k) > EEG_column(k-1))
            LBP(l) = '1';
        else
            LBP(l) = '0';
        end
        BIG_LBP(k-1,col) = bin2dec(LBP);
    end
end
LBP = BIG_LBP(1:end-1,:);

save("no_backup\Pat"+string(pat)+"\LBP_"+string(i)+".mat", "LBP")
end
end