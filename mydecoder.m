
clearvars -except  windowsize source torcver alphabetset Lmin
%filling rcver window
%rcverwindow=source(1:windowsize);
assert(strlength(source)~=0,'Generate source string')
assert(strlength(torcver)~=0,'Run encoder')
fprintf("LZ decoding started\n");
rcverwindow='';
for i=1:windowsize
    currentparselength=0;
    rcvedstr=alphabetset(bin2dec(torcver(1:ceil(log2(strlength(alphabetset)))))+1);
    currentparselength=currentparselength+ceil(log2(strlength(alphabetset)));
    rcverwindow=append(rcverwindow,rcvedstr);
    rcvedstr=[];
    torcver(1:currentparselength)=[];
end
totalrcvd=rcverwindow;
mlen=0;
totallength=0;
abc=0;
while ~isempty(torcver)
    currentparselength=0;
    k=strfind(torcver,'1');%finding the first 1 in unary binary code
    e=k(1);%consider the nearest match for parsing%also indicates no. of zeros present in unary binary code
    
    %parsing and decoding a match
    if e==1%rcved match of length 1
        currentparselength=currentparselength+1;%as n is 1
        %decode based on fixed length coding from alphabet set
        rcvedstr=alphabetset(bin2dec(torcver(2:ceil(log2(strlength(alphabetset)))+1))+1);
        %parsed this length from encoded string
        currentparselength=currentparselength+ceil(log2(strlength(alphabetset)));
    else%decode n and u based on the location first 1 of n
        %extending parse length till unary zeros
        currentparselength=currentparselength+e-1;
        %decoding n
        lenatrcver=bin2dec(torcver(e:e*2-1));
        %k=strfind(torcver(e*2:strlength(torcver)),'1');
        %r=k(1);
        %parsing and decoding u based on windowing size length
        distatrcver=bin2dec(torcver(e*2:e*2+ceil(log2(windowsize))-1));
        if distatrcver==0
            distatrcver=windowsize;
        end
        %decoding the orginal str based on decoded n and u
        if distatrcver<lenatrcver%if u < n implies that part of the string is sent and needs to be multiplied
            rcvedstrmain=rcverwindow(windowsize-distatrcver+1:windowsize);
            remlen=lenatrcver;
            while remlen~=0
                if remlen>=strlength(rcvedstrmain)
                    %add string present in window until the length to be
                    %added is less than u
                    rcvedstr=strcat(rcvedstr,rcvedstrmain);
                    remlen=remlen-strlength(rcvedstrmain);
                else
                    %if length to be added is less than u only part of the
                    %string found in window has to be added
                    rcvedstr=strcat(rcvedstr,rcvedstrmain(1:remlen));
                    remlen=0;
                end
            end
            rcvedstr;
        else
            rcvedstr=rcverwindow(windowsize-distatrcver+1:windowsize-distatrcver+lenatrcver);
        end
        currentparselength=currentparselength+ceil(log2(windowsize))+e;%extending parse to n,u pair
    end
    %adding the newly decoded str to the window
    rcverwindow=append(rcverwindow,rcvedstr);
    rcverwindow(1:strlength(rcvedstr))=[];
    %adding to output string
    totalrcvd=append(totalrcvd,rcvedstr);
    %totallength=totallength+strlength(rcvedstr);
    %abc=abc+1;
    %mlen=totallength/abc
    rcvedstr=[];
    %deleting the parsed part from the binary string
    torcver(1:currentparselength)=[];
end
totalrcvd;
isoutmatched=strcmp(totalrcvd,source);%check if decoded output matches with source
fprintf("\nIs decoded string matched with source=%d\n",isoutmatched);