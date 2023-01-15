rng(1);
alphabetset = 'a':char('a'+3);%param to change
fixedlengthentropy=2;%param to change
%numOfStates=length(alphabetset)+1;
numOfStates=4;%param to change
%creates a markov source with random transition matrix
mc=mcmix(numOfStates,Zeros=2*numOfStates);
fprintf("Is my markov chain ergodic=%d\n",isergodic(mc));%tells whether the generated chain is ergodic or not
requiredmatchlength=9;%param to change%typical match length
windowsize=2^(requiredmatchlength*fixedlengthentropy);%approx typical set size
numofLetters=windowsize*8;%source length is 8 times the window size
X = simulate(mc,numofLetters);
mc.P;%(i,j) element specifies prob from ith state to jth state
xFix = asymptotics(mc);%each state probability q(s)

%H[X/S] calculation
entropy=0;
for j=1:numOfStates   
    H=0;
    symbols = (1:numOfStates); % Alphabet vector
    prob = mc.P(j,:); % alphabet probability vector at a given state
    prob=prob(prob~=0);
    %[dict,avglen] = huffmandict((1:length(prob)),prob);
    %entropyPerState=avglen;%E[L]
    entropyPerState=log2(prob)*transpose(-prob);%H[X/s]
    entropy=entropy+xFix(j)*entropyPerState;% H[X/S]
end
%disp(entropy);
%creates a random symbol for each transition of state
%T=alphabetset(randi([1 length(alphabetset)],numOfStates,numOfStates))
%creates a symbol for each transition of state
T=[ 'abcd'
    'abcd'
    'abcd'
    'abcd'];
source='';
for i=2:length(X)
 %source=strcat(source,T(X(i-1),X(i)));%alphabetset(X(i)));
 source(i-1)=T(X(i-1),X(i));
end
%source=num2str(transpose(X));
%source=source(~isspace(source));
%disp(source);
save('sourcevars.mat','source','entropy','alphabetset','mc','X','xFix','windowsize','T');
%simplot(mc,X);
%figure;
%graphplot(mc,'LabelEdges',true,'ColorEdges',true,'ColorNodes',true);
