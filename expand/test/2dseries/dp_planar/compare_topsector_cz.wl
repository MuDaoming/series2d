SetDirectory[DirectoryName[$InputFileName]];

Get["myMod.wl"];
Get["topsector_C_z_data.wl"];
invS = Get["invS"];

B = 3;
Nprop = 8;
rhs = Join[ConstantArray[1, B], ConstantArray[0, Nprop]];

mmaCandZ = invS . rhs;
mmaCList = mmaCandZ[[1 ;; B]];
mmaC = Total[mmaCandZ[[1 ;; B]]];
mmaZ = mmaCandZ[[B + 1 ;; B + Nprop]];

Clear[toZpNumerator, sameZpQ, status];
toZpNumerator[expr_] := myMod[Expand[Numerator[Together[expr]]], p];
sameZpQ[a_, b_] := TrueQ[Expand[toZpNumerator[a - b]] === 0];
status[ok_] := If[ok, "OK", "DIFF"];

beforeCOK = sameZpQ[cppBeforeC, mmaC];
afterCOK = sameZpQ[cppAfterC, mmaC];
beforeCListOK = MapThread[sameZpQ, {cppBeforeCandZ[[1 ;; B]], mmaCList}];
afterCListOK = MapThread[sameZpQ, {cppAfterCandZ[[1 ;; B]], mmaCList}];
beforeZOK = MapThread[sameZpQ, {cppBeforeZ, mmaZ}];
afterZOK = MapThread[sameZpQ, {cppAfterZ, mmaZ}];
beforeCandZOK = MapThread[sameZpQ, {cppBeforeCandZ, mmaCandZ}];
afterCandZOK = MapThread[sameZpQ, {cppAfterCandZ, mmaCandZ}];

report = Join[
   {
    "p = " <> ToString[p],
    "mmaCandZLength = " <> ToString[Length[mmaCandZ]],
    "cppBeforeCandZLength = " <> ToString[Length[cppBeforeCandZ]],
    "cppAfterCandZLength = " <> ToString[Length[cppAfterCandZ]],
    "",
    "C before_convert vs MMA: " <> status[beforeCOK],
    "C after_convert  vs MMA: " <> status[afterCOK],
    ""
    },
   {"c_i before_convert vs MMA:"},
   Table[
    "c" <> ToString[i] <> ": " <> status[beforeCListOK[[i]]],
    {i, Length[beforeCListOK]}
    ],
   {""},
   {"c_i after_convert vs MMA:"},
   Table[
    "c" <> ToString[i] <> ": " <> status[afterCListOK[[i]]],
    {i, Length[afterCListOK]}
    ],
   {""},
   {"Z before_convert vs MMA:"},
   Table[
    "z" <> ToString[i] <> ": " <> status[beforeZOK[[i]]],
    {i, Length[beforeZOK]}
    ],
   {""},
   {"Z after_convert vs MMA:"},
   Table[
    "z" <> ToString[i] <> ": " <> status[afterZOK[[i]]],
    {i, Length[afterZOK]}
    ],
   {""},
   {
    "CandZ before total diff count = " <>
     ToString[Count[beforeCandZOK, False]],
    "CandZ before diff positions = " <>
     ToString[Flatten[Position[beforeCandZOK, False]]],
    "CandZ after  total diff count = " <>
     ToString[Count[afterCandZOK, False]],
    "CandZ after  diff positions = " <>
     ToString[Flatten[Position[afterCandZOK, False]]]
    }
   ];

Export["compare_topsector_cz_report.txt", StringRiffle[report, "\n"], "Text"];

Print[StringRiffle[report, "\n"]];
