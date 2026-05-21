If[Length[$ScriptCommandLine] =!= 3,
  Print["Usage: wolframscript -file de_to_wl.wl <dmaster_relations> <output_wl>"];
  Exit[1];
];

inputPath = $ScriptCommandLine[[2]];
outputPath = $ScriptCommandLine[[3]];

lines = StringTrim /@ Import[inputPath, "Lines"];
pLine = SelectFirst[lines, StringStartsQ["# p ="]];
If[MissingQ[pLine],
  Print["Error: input has no '# p =' header"];
  Exit[1];
];
p = ToExpression[StringTrim[StringDrop[pLine, StringLength["# p ="]]]];
ratBound = Floor[Sqrt[p/2]];

splitTopLevel[str_, sep_] := Module[
  {chars = Characters[str], depth = 0, cur = "", out = {}, ch},
  Do[
    ch = chars[[i]];
    Switch[ch,
      "(" | "[" | "{", depth++,
      ")" | "]" | "}", depth--
    ];
    If[ch === sep && depth === 0,
      AppendTo[out, StringTrim[cur]];
      cur = "",
      cur = cur <> ch
    ],
    {i, Length[chars]}
  ];
  If[StringTrim[cur] =!= "", AppendTo[out, StringTrim[cur]]];
  out
];

parsePoly[str_] := Module[{s, terms, coeffs = <||>, term, degree, coeff, parts},
  s = StringTrim[str];
  If[StringStartsQ[s, "("] && StringEndsQ[s, ")"],
    s = StringTake[s, {2, -2}]
  ];
  terms = splitTopLevel[s, "+"];
  Do[
    term = StringTrim[t];
    If[term === "", Continue[]];
    If[StringContainsQ[term, "delta"],
      parts = StringSplit[term, "*delta"];
      coeff = If[StringTrim[parts[[1]]] === "", 1, ToExpression[StringTrim[parts[[1]]]]];
      degree = If[StringContainsQ[term, "^"],
        ToExpression[StringTrim[StringSplit[term, "^"][[-1]]]],
        1
      ],
      coeff = ToExpression[term];
      degree = 0
    ];
    coeffs[degree] = Mod[Lookup[coeffs, degree, 0] + coeff, p],
    {t, terms}
  ];
  If[Length[coeffs] === 0, Return[0]];
  Sum[Lookup[coeffs, k, 0] delta^k, {k, 0, Max[Keys[coeffs]]}]
];

rationalReconstruct[a_Integer] := Module[
  {aa = Mod[a, p], u, v, q, den, num},
  If[aa == 0, Return[0]];
  u = {1, 0, p};
  v = {0, 1, aa};
  While[Abs[v[[3]]] > ratBound,
    q = Quotient[u[[3]], v[[3]]];
    {u, v} = {v, u - q v};
  ];
  num = v[[3]];
  den = v[[2]];
  If[den < 0, num = -num; den = -den];
  If[
    den == 0 || Abs[num] > ratBound || Abs[den] > ratBound ||
      GCD[num, den] =!= 1 || Mod[den aa - num, p] =!= 0,
    $Failed,
    num/den
  ]
];

reconstructPolynomial[expr_] := Module[{poly, coeffs, rec},
  poly = PolynomialMod[Expand[expr], p];
  coeffs = CoefficientList[poly, x];
  rec = rationalReconstruct /@ coeffs;
  If[MemberQ[rec, $Failed], Return[$Failed]];
  Sum[rec[[i + 1]] x^i, {i, 0, Length[rec] - 1}]
];

reconstructRational[expr_] := Module[{t, num, den, rnum, rden},
  If[expr === 0, Return[0]];
  t = Together[expr];
  num = Numerator[t];
  den = Denominator[t];
  rnum = reconstructPolynomial[num];
  rden = reconstructPolynomial[den];
  If[rnum === $Failed || rden === $Failed, $Failed, Together[rnum/rden]]
];

normalizeFiniteRational[expr_] := Module[{t, num, den},
  If[expr === 0, Return[0]];
  t = Together[expr];
  num = PolynomialMod[Expand[Numerator[t]], p];
  den = PolynomialMod[Expand[Denominator[t]], p];
  Together[num/den]
];

headerQ[line_] := StringStartsQ[line, "## d("];
masterFromHeader[line_] := StringTake[line, {StringLength["## d("] + 1, -2}];

blocks = {};
current = {};
Do[
  If[headerQ[line],
    If[current =!= {}, AppendTo[blocks, current]];
    current = {line},
    If[current =!= {} && line =!= "", AppendTo[current, line]]
  ],
  {line, lines}
];
If[current =!= {}, AppendTo[blocks, current]];
If[blocks === {},
  Print["Error: input has no derivative relation blocks"];
  Exit[1];
];

parseBlock[block_] := Module[
  {master, expr, terms, term, close, poly, label, assoc = <||>},
  master = masterFromHeader[First[block]];
  expr = StringRiffle[
    Select[Rest[block], !StringStartsQ[#, "#"] && # =!= "" &],
    " "
  ];
  expr = StringTrim[First[StringSplit[expr, "= 0"]]];
  If[StringEndsQ[expr, "+"], expr = StringDrop[expr, -1]];
  terms = splitTopLevel[expr, "+"];
  Do[
    term = StringTrim[t];
    If[term === "", Continue[]];
    close = StringPosition[term, ")*"];
    If[close === {},
      Print["Error: bad relation term: ", term];
      Exit[1];
    ];
    close = close[[1, 1]];
    poly = StringTake[term, close];
    label = StringTrim[StringDrop[term, close + 1]];
    assoc[label] = parsePoly[poly],
    {t, terms}
  ];
  <|"master" -> master, "terms" -> assoc|>
];

relations = parseBlock /@ blocks;
masters = relations[[All, "master"]];
masterIndex = AssociationThread[masters -> Range[Length[masters]]];

ADelta = ConstantArray[0, {Length[masters], Length[masters]}];
Do[
  dLabel = "d(" <> relations[[i, "master"]] <> ")";
  terms = relations[[i, "terms"]];
  If[!KeyExistsQ[terms, dLabel],
    Print["Error: missing derivative coefficient for ", dLabel];
    Exit[1];
  ];
  den = terms[dLabel];
  Do[
    label = key;
    If[label === dLabel, Continue[]];
    If[!KeyExistsQ[masterIndex, label],
      Print["Error: relation references non-master label: ", label];
      Exit[1];
    ];
    ADelta[[i, masterIndex[label]]] = -terms[label]/den,
    {key, Keys[terms]}
  ],
  {i, Length[relations]}
];

AXFinite = Together[ADelta /. delta -> 1 - x];
AX = Map[reconstructRational, AXFinite, {2}];
failures = Position[AX, $Failed];
outputComment = "";
If[failures =!= {},
  outputComment =
    "(* rational reconstruction failed at " <> ToString[failures, InputForm] <>
    "; output AX is over GF(p) with p = " <> ToString[p, InputForm] <>
    " and coefficients reduced modulo p; reconstruction bound was " <>
    ToString[ratBound, InputForm] <> " *)\n\n";
  AX = Map[normalizeFiniteRational, AXFinite, {2}],
  outputComment =
    "(* AX has been rationally reconstructed from GF(p), p = " <>
    ToString[p, InputForm] <> "; reconstruction bound was " <>
    ToString[ratBound, InputForm] <> " *)\n\n"
];

WriteString[
  outputPath,
  outputComment <>
  "masters = " <> ToString[masters, InputForm] <> ";\n\n" <>
  "AX = " <> ToString[AX, InputForm] <> ";\n"
];

Print["wrote ", outputPath, " with rational reconstruction bound ", ratBound];
