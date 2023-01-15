clear
load('sourcevars.mat');
windowstartptr=1;
windowendptr=windowsize;
encodeptr=1;
totalencodestr='';

%source string has to be greater than window size
assert(strlength(source)>windowsize,'source string length has to be greater than window size')
fprintf("LZ encoding started\n");
%filling the initial window for rcver
for i=1:windowsize
        %matchedlength=1;
        %ntobin=dec2bin(matchedlength,2*floor(log2(matchedlength))+1);
        matchstr=source(encodeptr);
        utobin=dec2bin(strfind(alphabetset,source(encodeptr))-1,ceil(log2(strlength(alphabetset))));
        %torcver=strcat(ntobin,utobin);
        torcver=utobin;
        %totalencodestr=strcat(totalencodestr,torcver);
        totalencodestr=append(totalencodestr,torcver);
        encodeptr=encodeptr+1;
end
initialprefixlength=strlength(totalencodestr);
while encodeptr~=strlength(source)+1
    hasMatchfound=0;
    matchedlength=1;%length of the match
    matchindex=1;%where match has occured in the window
    matchstr='';%string being matched
    %search by increasing the length of match str till window size
    for complength=2:windowsize
        out=[];%set of match indexes 
        if encodeptr+complength-1<=strlength(source)
            %virtually window length will be till last but one index of match str
            out=strfind(source(windowstartptr:windowendptr+complength-1),source(encodeptr:encodeptr+complength-1));
        else
            break;%implies no source str is present beyond this length for matching
        end
        if isempty(out)==1
            break;%implies no match for this in window
        else
            hasMatchfound=1;%record the matched str details and proceed for a bigger length
            matchedlength=complength;
            if ismatrix(out)==true%implies multiple matches are found
                matchindex=out(length(out));%taking the length closer to the window's end
            else%implies single match is found 
                matchindex=out;
            end
            matchstr=source(encodeptr:encodeptr+complength-1);
        end
    
    end
    if hasMatchfound
        [torcver,windowstartptr,windowendptr]=encodeInp(matchedlength,matchindex,windowstartptr,windowsize,source,alphabetset,encodeptr);
    else
        [torcver,windowstartptr,windowendptr]=encodeInp(1,0,windowstartptr,windowsize,source,alphabetset,encodeptr);
    end
    %{
    if hasMatchfound==1
        matchstr;
        ntobin=dec2bin(matchedlength,2*(floor(log2(matchedlength)))+1);%matched length is n %unary binary code
        matchdistfromencptr=windowsize- matchindex+1;%calculating u
        if matchdistfromencptr==windowsize
            utobin=dec2bin(0,ceil(log2(windowsize)));%assigning zero to the highest value 
        else
            utobin=dec2bin(matchdistfromencptr,ceil(log2(windowsize)));%fixed length encoding
        end
        torcver=strcat(ntobin,utobin);
        %sliding the window
        windowstartptr=windowstartptr+matchedlength;
        if strlength(source)>=windowstartptr+windowsize-1
            windowendptr=windowstartptr+windowsize-1;
        else
            windowendptr=strlength(source);
        end
    else%encoding using fixed length coding
        complength=1;
        matchedlength=1;%setting n to 1
        matchstr=source(encodeptr);
        ntobin=dec2bin(matchedlength,2*(floor(log2(matchedlength)))+1);%fixed length encoding
        utobin=dec2bin(strfind(alphabetset,source(encodeptr))-1,ceil(log2(strlength(alphabetset))));%fixed length encoding
        torcver=strcat(ntobin,utobin);
        windowstartptr=windowstartptr+complength;
        if strlength(source)>=windowstartptr+windowsize-1
            windowendptr=windowstartptr+windowsize-1;
        else
            windowendptr=strlength(source);
        end
    end
    %}
    encodeptr=encodeptr+matchedlength;
    %totalencodestr=strcat(totalencodestr,torcver);
    totalencodestr=append(totalencodestr,torcver);
    %Lmin=(strlength(totalencodestr))/(encodeptr)
end
torcver=totalencodestr;
%disp(totalencodestr)
Lmin=(strlength(torcver)-initialprefixlength)/(strlength(source)-windowsize);
%disp(Lmin);
save('Lmin.mat','Lmin','totalencodestr');
%this function converts n and u values in binary and updates window
%parameters
function [torcver,windowstartptr,windowendptr]=encodeInp(matchedlength,matchindex,windowstartptr,windowsize,source,alphabetset,encodeptr)
        %matchedlength
        ntobin=dec2bin(matchedlength,2*(floor(log2(matchedlength)))+1);%matched length is n %unary binary code
        matchdistfromencptr=windowsize- matchindex+1;%calculating u
        if matchedlength==1
            utobin=dec2bin(strfind(alphabetset,source(encodeptr))-1,ceil(log2(strlength(alphabetset))));%fixed length encoding
        else
            if matchdistfromencptr==windowsize
                utobin=dec2bin(0,ceil(log2(windowsize)));%assigning zero to the highest value 
            else
                utobin=dec2bin(matchdistfromencptr,ceil(log2(windowsize)));%fixed length encoding
            end
        end
        %torcver=strcat(ntobin,utobin);
        torcver=append(ntobin,utobin);
        %sliding the window
        windowstartptr=windowstartptr+matchedlength;
        if strlength(source)>=windowstartptr+windowsize-1
            windowendptr=windowstartptr+windowsize-1;
        else
            windowendptr=strlength(source);
        end
end




