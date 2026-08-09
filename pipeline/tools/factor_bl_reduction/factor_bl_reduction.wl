(* Factor one explicitly selected BL sector reduction over GF(p). *)

environmentArgs = Environment /@ {
  "FACTOR_BL_REDUCTION_INPUT",
  "FACTOR_BL_REDUCTION_OBJECT",
  "FACTOR_BL_REDUCTION_SECTOR",
  "FACTOR_BL_REDUCTION_OUTPUT"
};
scriptArgs = Which[
  Length[$ScriptCommandLine] == 5, Rest[$ScriptCommandLine],
  Length[$ScriptCommandLine] == 4, $ScriptCommandLine,
  AllTrue[environmentArgs, StringQ[#] && # =!= "" &], environmentArgs,
  True, {}
];
If[Length[scriptArgs] != 4,
  Print["Usage: wolframscript -file factor_bl_reduction.wl " <>
    "<reduction_output> <object_label> <sector> <output_path>"];
  Exit[1]
];

inputPath = scriptArgs[[1]];
wantedObject = StringReplace[StringTrim[scriptArgs[[2]]], Whitespace -> ""];
wantedSector = StringTrim[scriptArgs[[3]]];
outputPath = scriptArgs[[4]];

ClearAll[fail];
fail[msg_] := (Print["[factor_bl_reduction] Error: " <> ToString[msg]]; Exit[1]);

If[!FileExistsQ[inputPath], fail["Cannot open reduction output: " <> inputPath]];
lines = StringSplit[
  Import[inputPath, "Text"],
  RegularExpression["\\r\\n|\\n|\\r"]
];

primeLines = Select[lines, StringMatchQ[StringTrim[#], "# p =" ~~ ___] &];
If[Length[primeLines] != 1, fail["Expected exactly one '# p = ...' header"]];
p = Quiet@Check[
  ToExpression@StringTrim@Last@StringSplit[First[primeLines], "="],
  $Failed
];
If[!IntegerQ[p] || p <= 1, fail["Invalid finite-field prime"]];

ClearAll[normalizeSector];
normalizeSector[text_] := Module[{v},
  v = Quiet@Check[ToExpression[StringTrim[text]], $Failed];
  If[!ListQ[v] || !VectorQ[v, IntegerQ] || !AllTrue[v, (# == 0 || # == 1) &],
    fail["Invalid sector: " <> text]
  ];
  StringReplace[ToString[v, InputForm], " " -> ""]
];

wantedSector = normalizeSector[wantedSector];

sectorStarts = Flatten@Position[
  lines,
  s_String /; StringStartsQ[StringTrim[s], "sector="],
  {1}
];
If[sectorStarts === {}, fail["No sector reduction blocks found"]];

blocks = Table[
  start = sectorStarts[[i]];
  stop = If[i < Length[sectorStarts], sectorStarts[[i + 1]] - 1, Length[lines]];
  Take[lines, {start, stop}],
  {i, Length[sectorStarts]}
];

ClearAll[fieldValue];
fieldValue[block_, key_] := Module[{hits},
  hits = Select[block, StringStartsQ[StringTrim[#], key <> "="] &];
  If[Length[hits] == 0, Missing["NotFound"],
    StringTrim@StringDrop[StringTrim@First[hits], StringLength[key] + 1]
  ]
];

matches = Select[blocks,
  Quiet@Check[normalizeSector[fieldValue[#, "sector"]] == wantedSector, False] &&
  StringReplace[fieldValue[#, "object"], Whitespace -> ""] == wantedObject &
];
If[Length[matches] == 0,
  fail["No block matches object=" <> wantedObject <> " sector=" <> wantedSector]
];
If[Length[matches] != 1,
  fail["Multiple blocks match object=" <> wantedObject <> " sector=" <> wantedSector]
];
block = First[matches];

status = fieldValue[block, "status"];
If[status =!= "success", fail["Selected block has status=" <> ToString[status]]];
If[MemberQ[StringTrim /@ block, "zero"], fail["Selected block is a zero contribution"]];
If[fieldValue[block, "free_master"] === "1", fail["Selected block is a free master"]];

ClearAll[parseCoeffList];
parseCoeffList[text_] := Module[{v},
  v = Quiet@Check[ToExpression[StringTrim[text]], $Failed];
  If[!ListQ[v] || !VectorQ[v, IntegerQ] || Length[v] == 0,
    fail["Invalid polynomial coefficient list"]
  ];
  Mod[v, p]
];

denText = fieldValue[block, "den"];
If[MissingQ[denText], fail["Selected block has no denominator"]];
denCoeffs = parseCoeffList[denText];

termLines = Select[block, StringStartsQ[StringTrim[#], "term "] &];
If[termLines === {}, fail["Selected block has no reduction terms"]];
terms = Map[
  Function[line,
    eq = StringPosition[line, "="];
    If[Length[eq] != 1, fail["Invalid term line"]];
    pos = First@First@eq;
    <|
      "label" -> StringTrim@StringTake[line, {6, pos - 1}],
      "coeffs" -> parseCoeffList@StringDrop[line, pos]
    |>
  ],
  termLines
];

ClearAll[polyFromCoeffs, coeffsFromPoly];
polyFromCoeffs[c_] := Sum[c[[i + 1]] delta^i, {i, 0, Length[c] - 1}];
coeffsFromPoly[poly_] := Module[{d},
  d = Exponent[poly, delta];
  If[d < 0, {}, Mod[CoefficientList[poly, delta], p]]
];

ClearAll[factorData];
factorData[coeffs_] := Module[{poly, fl, unit, pairs, normalizedPairs},
  poly = PolynomialMod[polyFromCoeffs[coeffs], p];
  If[TrueQ[poly == 0], fail["Cannot factor a zero polynomial"]];
  fl = FactorList[poly, Modulus -> p];
  unit = Mod[fl[[1, 1]], p];
  pairs = Rest[fl];
  normalizedPairs = Map[
    Function[pair,
      fac = PolynomialMod[pair[[1]], p];
      lc = Mod[Coefficient[fac, delta, Exponent[fac, delta]], p];
      monic = PolynomialMod[fac PowerMod[lc, -1, p], p];
      unit = Mod[unit PowerMod[lc, pair[[2]], p], p];
      {coeffsFromPoly[monic], pair[[2]]}
    ],
    pairs
  ];
  <|"unit" -> unit, "factors" -> normalizedPairs,
    "degree" -> Exponent[poly, delta], "poly" -> poly|>
];

denFact = factorData[denCoeffs];
termFacts = Map[Append[#, "factorization" -> factorData[#["coeffs"]]] &, terms];

ClearAll[factorKey];
factorKey[v_] := StringReplace[ToString[v, InputForm], " " -> ""];
allFactorCoeffs = Join[
  denFact["factors"][[All, 1]],
  Flatten[(#["factorization"]["factors"][[All, 1]]) & /@ termFacts, 1]
];
uniqueFactorCoeffs = DeleteDuplicates[allFactorCoeffs];
factorIds = Association@Table[factorKey[uniqueFactorCoeffs[[i]]] -> ("F" <> ToString[i]),
  {i, Length[uniqueFactorCoeffs]}];

ClearAll[toExponentAssociation];
toExponentAssociation[pairs_] := Association@Map[
  (factorIds[factorKey[#[[1]]]] -> #[[2]]) &,
  pairs
];
denExponents = toExponentAssociation[denFact["factors"]];

structuredTerms = Map[
  Function[t,
    numExponents = toExponentAssociation[t["factorization"]["factors"]];
    ids = Union[Keys[numExponents], Keys[denExponents]];
    net = Association@Table[
      id -> (Lookup[numExponents, id, 0] - Lookup[denExponents, id, 0]),
      {id, ids}
    ];
    <|
      "label" -> t["label"],
      "unit" -> Mod[t["factorization"]["unit"] PowerMod[denFact["unit"], -1, p], p],
      "numerator" -> Select[net, # > 0 &],
      "denominator" -> Association@KeyValueMap[#1 -> -#2 &, Select[net, # < 0 &]],
      "rawDegree" -> t["factorization"]["degree"]
    |>
  ],
  termFacts
];

usedIds = Union@Flatten[
  Join[Keys[#["numerator"]], Keys[#["denominator"]]] & /@ structuredTerms
];
usedFactors = Select[
  Table[
    <|"id" -> ("F" <> ToString[i]), "coeffs" -> uniqueFactorCoeffs[[i]],
      "degree" -> (Length[uniqueFactorCoeffs[[i]]] - 1)|>,
    {i, Length[uniqueFactorCoeffs]}
  ],
  MemberQ[usedIds, #["id"]] &
];

factorParameterCount = Total[Lookup[usedFactors, "degree", {}]];
termConstantCount = Length[structuredTerms];
structuredUnknownCount = factorParameterCount + termConstantCount;
rawUnknownCount = denFact["degree"] + Total[(#["rawDegree"] + 1) & /@ structuredTerms];
minimumSeriesCoefficients = structuredUnknownCount;
minimumExpansionDegree = Max[-1, minimumSeriesCoefficients - 1];

ClearAll[rebuildPolynomial, factorizationVerifiedQ];
rebuildPolynomial[fd_] := PolynomialMod[
  fd["unit"] Times @@ ((polyFromCoeffs[#[[1]]]^#[[2]]) & /@ fd["factors"]),
  p
];
factorizationVerifiedQ[fd_] := TrueQ[PolynomialMod[fd["poly"] - rebuildPolynomial[fd], p] == 0];
factorizationVerified = factorizationVerifiedQ[denFact] &&
  AllTrue[termFacts, factorizationVerifiedQ[#["factorization"]] &];
If[!factorizationVerified, fail["Internal factorization verification failed"]];

idToCoeffs = Association@Table[
  ("F" <> ToString[i]) -> uniqueFactorCoeffs[[i]],
  {i, Length[uniqueFactorCoeffs]}
];
ClearAll[polyFromExponents];
polyFromExponents[a_Association] := PolynomialMod[
  Times @@ KeyValueMap[(polyFromCoeffs[idToCoeffs[#1]]^#2) &, a],
  p
];
structureVerified = And @@ MapThread[
  Function[{raw, structured},
    structuredNumerator = PolynomialMod[
      structured["unit"] polyFromExponents[structured["numerator"]],
      p
    ];
    structuredDenominator = polyFromExponents[structured["denominator"]];
    TrueQ[PolynomialMod[
      raw["factorization"]["poly"] structuredDenominator -
        denFact["poly"] structuredNumerator,
      p
    ] == 0]
  ],
  {termFacts, structuredTerms}
];
If[!structureVerified, fail["Internal cancelled-structure verification failed"]];

ClearAll[listText, exponentText];
listText[v_] := StringReplace[ToString[v, InputForm], " " -> ""];
exponentText[a_Association] := If[Length[a] == 0, "1",
  StringRiffle[KeyValueMap[#1 <> "^" <> ToString[#2] &, a], " * "]
];

out = OpenWrite[outputPath];
If[out === $Failed, fail["Cannot open output path: " <> outputPath]];
WriteString[out,
  "# format = bl_factor_structure_v1\n",
  "# source = ", inputPath, "\n",
  "# p = ", ToString[p], "\n",
  "# object = ", wantedObject, "\n",
  "# sector = ", wantedSector, "\n\n",
  "[summary]\n",
  "terms=", ToString[Length[structuredTerms]], "\n",
  "raw_den_degree=", ToString[denFact["degree"]], "\n",
  "raw_unknown_coefficients=", ToString[rawUnknownCount], "\n",
  "unique_factors=", ToString[Length[usedFactors]], "\n",
  "factor_coefficients=", ToString[factorParameterCount], "\n",
  "term_constants=", ToString[termConstantCount], "\n",
  "structured_unknown_coefficients=", ToString[structuredUnknownCount], "\n",
  "minimum_series_coefficients=", ToString[minimumSeriesCoefficients], "\n",
  "minimum_expansion_degree=", ToString[minimumExpansionDegree], "\n",
  "factorization_verified=1\n",
  "structure_verified=1\n",
  "rank_check_required=1\n\n",
  "[factors]\n"
];
Do[
  WriteString[out,
    f["id"], " degree=", ToString[f["degree"]],
    " coeffs=", listText[f["coeffs"]], "\n"
  ],
  {f, usedFactors}
];
WriteString[out, "\n[terms]\n"];
Do[
  WriteString[out,
    "term=", t["label"], "\n",
    "unit=", ToString[t["unit"]], "\n",
    "numerator=", exponentText[t["numerator"]], "\n",
    "denominator=", exponentText[t["denominator"]], "\n\n"
  ],
  {t, structuredTerms}
];
Close[out];

Print["Wrote factor structure to: " <> outputPath];
Print["structured_unknown_coefficients = " <> ToString[structuredUnknownCount]];
Print["minimum_expansion_degree = " <> ToString[minimumExpansionDegree]];
