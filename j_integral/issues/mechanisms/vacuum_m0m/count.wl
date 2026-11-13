$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
BLFamilyDefine[vacm0m,26/7,{l1^2-1,l2^2,(l1+l2)^2-1},{l1,l2},{},{},{},{1,1,1},{},{1,1,1}];
start=AbsoluteTime[];mi=BLMaximalCutMasters[{1,1,1},{},"PreferredMI"->"Dot"];
Export["result.json",<|"masters"->ToString[mi,InputForm],"count"->Length[mi],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];
