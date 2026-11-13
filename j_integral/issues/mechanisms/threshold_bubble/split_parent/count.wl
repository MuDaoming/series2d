$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
prop={l1^2-1,l2^2,(l1+l2)^2-5/4,(l1+l2+Q)^2-1,(l1+Q)^2};top={1,1,1,1,0};
BLFamilyDefine[thresholdsplit,26/7,prop,{l1,l2},{Q},{},{Q^2->2},top,{},top];
start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Export["result.json",<|"count"->Length[mi],"masters"->ToString[mi,InputForm],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];
