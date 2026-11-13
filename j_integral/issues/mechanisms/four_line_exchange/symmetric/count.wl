$HistoryLength=0;out=DirectoryName[$InputFileName];SetDirectory[out];
Get["/home/mudaoming/Documents/Packages/blade/blade-main/Blade.wl"];
Unprotect[Blade`Private`CacheDirectory];Clear[Blade`Private`CacheDirectory];
Blade`Private`CacheDirectory[key_]:=FileNameJoin[{out,"cache",ToString[Blade`Private`BLFamily,InputForm]<>"_"<>key}];
BLNthreads=1;BLSSolveNthreads=1;BLSetReducerOptions["IntegralOrdering"->2];
prop={l1^2-(0),l2^2,(l2+Q)^2,(l1+l2+P)^2,(l1+P)^2,(l1+Q)^2,(l2+P)^2};
sp={P^2->149/60,Q^2->38113/4620,P Q->22483/4620};top={1,1,1,1,0,0,0};
BLFamilyDefine[foursymmetric,26/7,prop,{l1,l2},{P,Q},{},sp,top,{},top];
start=AbsoluteTime[];mi=BLMaximalCutMasters[top,{},"PreferredMI"->"Dot"];Put[mi,"masters.wl"];Export["result.json",<|"masters"->ToString[mi,InputForm],"count"->Length[mi],"seconds"->(AbsoluteTime[]-start)|>,"RawJSON"];Print["RESULT=",mi];Exit[];
