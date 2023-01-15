mysource
myencoder
mydecoder
load("sourcevars.mat")
load('Lmin.mat')
fprintf("Entropy H[X/S] of dependent source=%f\n",entropy);
fprintf("Lmin with window elements included=%f\n",length(totalencodestr)/length(source));
fprintf("Lmin without window elements included=%f\n",Lmin);
