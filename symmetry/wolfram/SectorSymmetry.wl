(* ::Package:: *)

BeginPackage["SectorSymmetry`"];

FindSectorRelations::usage =
  "FindSectorRelations[{U,F}, x] finds sector relations induced by permutations of x[i].";

PakCanonicalizeUF::usage =
  "PakCanonicalizeUF[{U,F}, vars] returns the joint Pak canonical form and witnesses.";

RestrictToSector::usage =
  "RestrictToSector[{U,F}, sector, vars] sets parameters outside sector to zero.";

VerifySectorMap::usage =
  "VerifySectorMap[source, target, rules, {U,F}, vars] verifies a proposed parameter map.";

Options[FindSectorRelations] = {
  "IncludeDegenerateSectors" -> False,
  "VerifyMappings" -> True
};

Begin["`Private`"];

ClearAll[
  coefficientKey, enumerateSectors, extractParameterVariables,
  polynomialRows, prefixKey, pakSearch, canonicalizeSector,
  buildRelation
];

coefficientKey[coefficient_] :=
  ToString[InputForm[Together[FactorTerms[Expand[coefficient]]]]];

extractParameterVariables[polynomials_List, parameterHead_Symbol] := Module[
  {indices, expected},
  indices = DeleteDuplicates[
    Cases[
      polynomials,
      HoldPattern[parameterHead[index_Integer?Positive]] :> index,
      Infinity
    ]
  ];
  If[indices === {},
    Return[{}]
  ];
  indices = Sort[indices];
  expected = Range[Max[indices]];
  If[indices =!= expected,
    Message[FindSectorRelations::params, parameterHead, indices];
    Return[$Failed]
  ];
  parameterHead /@ expected
];

FindSectorRelations::params =
  "Parameters with head `1` must have consecutive positive indices starting at 1; found `2`.";

enumerateSectors[count_Integer?NonNegative] :=
  Rest[IntegerDigits[Range[0, 2^count - 1], 2, count]];

RestrictToSector[
  polynomials : {_, _},
  sector : {(_Integer)..},
  vars_List
] /; Length[sector] === Length[vars] := Module[
  {rules},
  rules = Thread[Pick[vars, sector, 0] -> 0];
  Expand[polynomials /. rules]
];

polynomialRows[polynomial_, vars_List, tag_Integer] := Module[
  {rules},
  rules = CoefficientRules[Expand[polynomial], vars];
  ({tag, coefficientKey[Last[#]], First[#]} &) /@ rules
];

prefixKey[rows_List, order_List] := ToString[
  InputForm[
    Sort[
      ({#[[1]], #[[2]], #[[3, order]]} &) /@ rows
    ]
  ]
];

pakSearch[rows_List, variableCount_Integer] := Module[
  {states, expanded, keys, bestKey, finalKey},
  states = {{}};
  Do[
    expanded = DeleteDuplicates[
      Flatten[
        Table[
          Append[state, candidate],
          {state, states},
          {candidate, Complement[Range[variableCount], state]}
        ],
        1
      ]
    ];
    keys = prefixKey[rows, #] & /@ expanded;
    bestKey = First[Sort[keys]];
    states = Pick[expanded, keys, bestKey],
    {variableCount}
  ];
  finalKey = If[states === {{}}, prefixKey[rows, {}], prefixKey[rows, First[states]]];
  <|
    "Key" -> finalKey,
    "Orders" -> states
  |>
];

PakCanonicalizeUF[
  polynomials : {_, _},
  vars_List
] := Module[
  {rows, result},
  rows = Join[
    polynomialRows[polynomials[[1]], vars, 0],
    polynomialRows[polynomials[[2]], vars, 1]
  ];
  result = pakSearch[rows, Length[vars]];
  Join[
    result,
    <|"Rows" -> rows|>
  ]
];

canonicalizeSector[
  polynomials : {_, _},
  sector_List,
  vars_List
] := Module[
  {restricted, activePositions, activeVars, canonical},
  restricted = RestrictToSector[polynomials, sector, vars];
  activePositions = Flatten[Position[sector, 1]];
  activeVars = vars[[activePositions]];
  canonical = PakCanonicalizeUF[restricted, activeVars];
  <|
    "Sector" -> sector,
    "RestrictedPolynomials" -> restricted,
    "ActivePositions" -> activePositions,
    "ActiveVariables" -> activeVars,
    "Degenerate" -> (
      TrueQ[Expand[restricted[[1]]] === 0] ||
      TrueQ[Expand[restricted[[2]]] === 0]
    ),
    "CanonicalKey" -> canonical["Key"],
    "CanonicalOrders" -> canonical["Orders"]
  |>
];

VerifySectorMap[
  sourceSector_List,
  targetSector_List,
  rules_List,
  polynomials : {_, _},
  vars_List
] := Module[
  {source, target, differences},
  source = RestrictToSector[polynomials, sourceSector, vars];
  target = RestrictToSector[polynomials, targetSector, vars];
  differences = Together[Expand[#]] & /@ ((source /. rules) - target);
  <|
    "UVerified" -> TrueQ[differences[[1]] === 0],
    "FVerified" -> TrueQ[differences[[2]] === 0],
    "Differences" -> differences
  |>
];

buildRelation[
  source_Association,
  target_Association,
  polynomials_List,
  vars_List,
  verifyQ_
] := Module[
  {sourceOrder, targetOrder, rules, verification},
  sourceOrder = First[source["CanonicalOrders"]];
  targetOrder = First[target["CanonicalOrders"]];
  rules = Thread[
    source["ActiveVariables"][[sourceOrder]] ->
    target["ActiveVariables"][[targetOrder]]
  ];
  verification = If[
    TrueQ[verifyQ],
    VerifySectorMap[
      source["Sector"], target["Sector"], rules, polynomials, vars
    ],
    <|"UVerified" -> Missing["NotChecked"], "FVerified" -> Missing["NotChecked"]|>
  ];
  <|
    "Sector" -> source["Sector"],
    "Representative" -> target["Sector"],
    "ParameterRules" -> rules,
    "UVerified" -> verification["UVerified"],
    "FVerified" -> verification["FVerified"]
  |>
];

FindSectorRelations[
  polynomials : {_, _},
  parameterHead_Symbol,
  OptionsPattern[]
] := Module[
  {
    vars, sectors, data, degenerate, retained, grouped, orbits,
    relations, includeDegenerateQ, verifyQ
  },
  vars = extractParameterVariables[polynomials, parameterHead];
  If[vars === $Failed, Return[$Failed]];

  sectors = enumerateSectors[Length[vars]];
  data = canonicalizeSector[polynomials, #, vars] & /@ sectors;
  degenerate = Select[data, TrueQ[#["Degenerate"]] &];

  includeDegenerateQ = TrueQ[OptionValue["IncludeDegenerateSectors"]];
  retained = If[
    includeDegenerateQ,
    data,
    Select[data, ! TrueQ[#["Degenerate"]] &]
  ];

  grouped = GatherBy[
    retained,
    {Total[#["Sector"]], #["CanonicalKey"]} &
  ];
  grouped = SortBy[grouped, First[Sort[#[[All, "Sector"]]]] &];

  verifyQ = OptionValue["VerifyMappings"];
  orbits = Map[
    Function[group,
      Module[{members, representative, target, orbitRelations},
        members = Sort[group[[All, "Sector"]]];
        representative = Last[members];
        target = First[Select[group, #["Sector"] === representative &]];
        orbitRelations = buildRelation[#, target, polynomials, vars, verifyQ] & /@ group;
        <|
          "Representative" -> representative,
          "Members" -> members,
          "Relations" -> SortBy[orbitRelations, #["Sector"] &]
        |>
      ]
    ],
    grouped
  ];

  relations = Flatten[orbits[[All, "Relations"]], 1];
  <|
    "Parameters" -> vars,
    "SectorCount" -> Length[sectors],
    "DegenerateSectors" -> Sort[degenerate[[All, "Sector"]]],
    "Representatives" -> orbits[[All, "Representative"]],
    "Orbits" -> orbits,
    "Relations" -> SortBy[relations, #["Sector"] &]
  |>
];

End[];
EndPackage[];
