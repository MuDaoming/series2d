$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
prop={(l1+P)^2-1,l1^2-1,l2^2-(5/4),(l1+l2+Q)^2-1,(l1+Q)^2,(l2+P)^2,(l2+Q)^2};sp={P^2->1,Q^2->1,P Q->5/24};top={1,1,1,1,0,0,0};
BLFamilyDefine[dpbubblesplit,26/7,prop,{l1,l2},{P,Q},{},sp,top,{},top];
start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Put[mi,"masters.wl"];Export["result.json",<|"masters"->ToString[mi,InputForm],"count"->Length[mi],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];
