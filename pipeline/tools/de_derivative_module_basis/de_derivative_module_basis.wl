BeginPackage["DEDerivativeModuleBasis`"];

BuildDerivativeMatrices::usage =
  "BuildDerivativeMatrices[A,x,r] returns B_0,...,B_r with M^(n)=B_n M.";
BuildDerivativeModuleBasis::usage =
  "BuildDerivativeModuleBasis[A,x,r,Modulus->p] clears the common " <>
  "denominator of the rows of B_0,...,B_r and computes a weak-Popov " <>
  "basis of their F_p[x]-row module.";
CoordinatesInDerivativeModule::usage =
  "CoordinatesInDerivativeModule[result,target,x] reduces one rational " <>
  "target row by a basis returned by BuildDerivativeModuleBasis and " <>
  "returns its polynomial coordinates.";

Options[BuildDerivativeModuleBasis] = {
  Modulus -> 2305843009213693951,
  "Verbose" -> True,
  "ColumnShifts" -> Automatic
};

Begin["`Private`"];

ClearAll[
  fieldScalar, toFieldPolynomial, rationalFieldData,
  polynomialLCMList, rowDegree, leadingPosition, leadingCoefficient,
  zeroRowQ, normalizeRow, weakPopovBasis
];

fieldScalar[value_, p_Integer] := Module[
  {rational = Together[value], numerator, denominator},
  numerator = Mod[Numerator[rational], p];
  denominator = Mod[Denominator[rational], p];
  If[denominator === 0, Return[$Failed]];
  Mod[numerator PowerMod[denominator, -1, p], p]
];

toFieldPolynomial[polynomial_, variable_Symbol, p_Integer] := Module[
  {expanded = Expand[polynomial], degree},
  If[expanded === 0, Return[0]];
  degree = Exponent[expanded, variable];
  PolynomialMod[
    Sum[
      fieldScalar[Coefficient[expanded, variable, power], p] variable^power,
      {power, 0, degree}
    ],
    p
  ]
];

rationalFieldData[expression_, variable_Symbol, p_Integer] := Module[
  {rational, numerator, denominator},
  rational = Together[expression];
  numerator = toFieldPolynomial[Numerator[rational], variable, p];
  denominator = toFieldPolynomial[Denominator[rational], variable, p];
  If[numerator === $Failed || denominator === $Failed || denominator === 0,
    Return[$Failed]
  ];
  <|"Numerator" -> numerator, "Denominator" -> denominator|>
];

polynomialLCMList[polynomials_List, variable_Symbol, p_Integer] :=
  Fold[
    PolynomialLCM[#1, #2, variable, Modulus -> p] &,
    1,
    polynomials
  ];

zeroRowQ[row_List] := And @@ (TrueQ[# === 0] & /@ row);

normalizedShifts[row_List, shifts_] := Replace[
  shifts,
  Automatic :> ConstantArray[0, Length[row]]
];

rowDegree[row_List, variable_Symbol, shifts_:Automatic] := If[
  zeroRowQ[row],
  -Infinity,
  Max[MapThread[
    If[#1 === 0, -Infinity, Exponent[#1, variable] + #2] &,
    {row, normalizedShifts[row, shifts]}
  ]]
];

leadingPosition[row_List, variable_Symbol, shifts_:Automatic] := Module[
  {degree, localShifts, shiftedDegrees},
  If[zeroRowQ[row], Return[0]];
  localShifts = normalizedShifts[row, shifts];
  shiftedDegrees = MapThread[
    If[#1 === 0, -Infinity, Exponent[#1, variable] + #2] &,
    {row, localShifts}
  ];
  degree = Max[shiftedDegrees];
  Last@Flatten@Position[
    shiftedDegrees,
    degree
  ]
];

leadingCoefficient[row_List, variable_Symbol, shifts_:Automatic] := Module[
  {position, degree, localShifts},
  localShifts = normalizedShifts[row, shifts];
  position = leadingPosition[row, variable, localShifts];
  If[position === 0, Return[0]];
  degree = rowDegree[row, variable, localShifts] - localShifts[[position]];
  Coefficient[row[[position]], variable, degree]
];

normalizeRow[
  row_List, variable_Symbol, p_Integer, shifts_:Automatic
] := Module[
  {coefficient},
  If[zeroRowQ[row], Return[row]];
  coefficient = Mod[leadingCoefficient[row, variable, shifts], p];
  PolynomialMod[
    PowerMod[coefficient, -1, p] row,
    p
  ]
];

weakPopovBasis[
  inputRows_List, variable_Symbol, p_Integer, verbose_,
  shifts_:Automatic
] := Module[
  {rows, transformation, iterations = 0, positions, duplicate, indices,
   low, high, lowDegree, highDegree, factor, order, coefficient,
   localShifts},
  rows = PolynomialMod[#, p] & /@ inputRows;
  localShifts = normalizedShifts[First[rows], shifts];
  transformation = IdentityMatrix[Length[rows]];
  While[True,
    positions = leadingPosition[#, variable, localShifts] & /@ rows;
    duplicate = Select[
      DeleteDuplicates[DeleteCases[positions, 0]],
      Count[positions, #] > 1 &
    ];
    If[duplicate === {}, Break[]];
    indices = Flatten@Position[positions, First[duplicate]];
    {low, high} = Take[indices, 2];
    lowDegree = rowDegree[rows[[low]], variable, localShifts];
    highDegree = rowDegree[rows[[high]], variable, localShifts];
    If[lowDegree > highDegree,
      {low, high} = {high, low};
      {lowDegree, highDegree} = {highDegree, lowDegree};
    ];
    factor = Mod[
      leadingCoefficient[rows[[high]], variable, localShifts] *
        PowerMod[
          Mod[
            leadingCoefficient[rows[[low]], variable, localShifts],
            p
          ],
          -1,
          p
        ],
      p
    ] variable^(highDegree - lowDegree);
    rows[[high]] = PolynomialMod[
      rows[[high]] - factor rows[[low]],
      p
    ];
    transformation[[high]] = PolynomialMod[
      transformation[[high]] - factor transformation[[low]],
      p
    ];
    iterations++;
    If[verbose && Mod[iterations, 1000] === 0,
      Print["  weak-Popov row reductions: ", iterations]
    ];
  ];
  order = Select[Range[Length[rows]], !zeroRowQ[rows[[#]]] &];
  order = SortBy[
    order,
    leadingPosition[rows[[#]], variable, localShifts] &
  ];
  rows = rows[[order]];
  transformation = transformation[[order]];
  Do[
    coefficient = Mod[
      leadingCoefficient[rows[[index]], variable, localShifts], p
    ];
    rows[[index]] = PolynomialMod[
      PowerMod[coefficient, -1, p] rows[[index]], p
    ];
    transformation[[index]] = PolynomialMod[
      PowerMod[coefficient, -1, p] transformation[[index]], p
    ],
    {index, Length[rows]}
  ];
  <|
    "Rows" -> rows,
    "GeneratorCoordinates" -> transformation,
    "Iterations" -> iterations,
    "LeadingPositions" ->
      (leadingPosition[#, variable, localShifts] & /@ rows),
    "RowDegrees" -> (rowDegree[#, variable, localShifts] & /@ rows),
    "ColumnShifts" -> localShifts
  |>
];

CoordinatesInDerivativeModule[
  result_Association, targetRow_List, variable_Symbol
] := Module[
  {p, denominator, basis, fieldData, polynomialTarget, remainder,
   coordinates, position, degree, basisIndex, basisDegree, factor,
   steps = 0, quotient, badDenominators, shifts},
  p = result["Modulus"];
  denominator = result["CommonDenominator"];
  basis = result["PolynomialBasis"];
  shifts = Lookup[
    result, "ColumnShifts", ConstantArray[0, Length[targetRow]]
  ];
  fieldData = rationalFieldData[#, variable, p] & /@ targetRow;
  If[!FreeQ[fieldData, $Failed],
    Return[Failure[
      "BadCoefficient",
      <|"Message" -> "A target coefficient is undefined modulo p."|>
    ]]
  ];
  badDenominators = Select[
    fieldData,
    PolynomialRemainder[
      denominator, #["Denominator"], variable, Modulus -> p
    ] =!= 0 &
  ];
  If[badDenominators =!= {},
    Return[Failure[
      "OutsideCommonDenominator",
      <|"Denominators" -> Lookup[badDenominators, "Denominator"]|>
    ]]
  ];
  polynomialTarget = Map[
    Function[entry,
      quotient = PolynomialQuotient[
        denominator, entry["Denominator"], variable, Modulus -> p
      ];
      PolynomialMod[entry["Numerator"] quotient, p]
    ],
    fieldData
  ];
  If[FailureQ[polynomialTarget], Return[polynomialTarget]];
  remainder = polynomialTarget;
  coordinates = ConstantArray[0, Length[basis]];
  While[!zeroRowQ[remainder],
    position = leadingPosition[remainder, variable, shifts];
    degree = rowDegree[remainder, variable, shifts];
    basisIndex = FirstPosition[
      leadingPosition[#, variable, shifts] & /@ basis,
      position,
      Missing["NotFound"]
    ];
    If[MissingQ[basisIndex],
      Return[<|
        "MemberQ" -> False,
        "Reason" -> "NoBasisLeadingPosition",
        "Remainder" -> remainder
      |>]
    ];
    basisIndex = First[basisIndex];
    basisDegree = rowDegree[basis[[basisIndex]], variable, shifts];
    If[basisDegree > degree,
      Return[<|
        "MemberQ" -> False,
        "Reason" -> "RemainderDegreeBelowBasisDegree",
        "Remainder" -> remainder
      |>]
    ];
    factor = Mod[
      leadingCoefficient[remainder, variable, shifts] *
        PowerMod[
          Mod[
            leadingCoefficient[
              basis[[basisIndex]], variable, shifts
            ],
            p
          ],
          -1,
          p
        ],
      p
    ] variable^(degree - basisDegree);
    coordinates[[basisIndex]] = PolynomialMod[
      coordinates[[basisIndex]] + factor,
      p
    ];
    remainder = PolynomialMod[
      remainder - factor basis[[basisIndex]],
      p
    ];
    steps++;
  ];
  <|
    "MemberQ" -> True,
    "Coordinates" -> coordinates,
    "CoordinateDegrees" ->
      (If[# === 0, -Infinity, Exponent[#, variable]] & /@ coordinates),
    "ScalarCoefficientCount" -> Total[
      If[# === 0, 0, Exponent[#, variable] + 1] & /@ coordinates
    ],
    "ReductionSteps" -> steps
  |>
];

BuildDerivativeMatrices[matrix_, variable_Symbol, maximumOrder_Integer] :=
  NestList[
    Map[Together, D[#, variable] + # . matrix, {2}] &,
    IdentityMatrix[Length[matrix]],
    maximumOrder
  ];

BuildDerivativeModuleBasis[
  matrix_, variable_Symbol, maximumOrder_Integer, OptionsPattern[]
] := Module[
  {p, verbose, derivativeMatrices, dimension, generatorRows,
   fieldData, commonDenominator, polynomialRows, popov, columnShifts},
  p = OptionValue[Modulus];
  verbose = TrueQ[OptionValue["Verbose"]];
  dimension = Length[matrix];
  columnShifts = Replace[
    OptionValue["ColumnShifts"],
    Automatic :> ConstantArray[0, dimension]
  ];
  If[
    Length[columnShifts] =!= dimension ||
      !VectorQ[columnShifts, IntegerQ],
    Return[Failure[
      "InvalidColumnShifts",
      <|"Message" -> "ColumnShifts must be an integer vector of length N."|>
    ]]
  ];
  If[verbose, Print["building B_0,...,B_", maximumOrder]];
  derivativeMatrices = BuildDerivativeMatrices[
    matrix, variable, maximumOrder
  ];
  generatorRows = Join @@ derivativeMatrices;
  If[verbose, Print["converting rational rows to F_p(x)"]];
  fieldData = Map[
    rationalFieldData[#, variable, p] &,
    generatorRows,
    {2}
  ];
  If[!FreeQ[fieldData, $Failed],
    Return[Failure[
      "BadCoefficient",
      <|"Message" -> "A coefficient denominator vanishes modulo p."|>
    ]]
  ];
  commonDenominator = polynomialLCMList[
    Flatten[fieldData[[All, All, "Denominator"]]],
    variable,
    p
  ];
  If[verbose,
    Print[
      "common denominator degree = ",
      Exponent[commonDenominator, variable]
    ]
  ];
  polynomialRows = Map[
    Function[row,
      Map[
        Function[entry,
          PolynomialMod[
            entry["Numerator"] *
              PolynomialQuotient[
                commonDenominator,
                entry["Denominator"],
                variable,
                Modulus -> p
              ],
            p
          ]
        ],
        row
      ]
    ],
    fieldData
  ];
  If[verbose, Print["reducing polynomial row module"]];
  popov = weakPopovBasis[
    polynomialRows, variable, p, verbose, columnShifts
  ];
  <|
    "Modulus" -> p,
    "Variable" -> variable,
    "DerivativeOrder" -> maximumOrder,
    "Dimension" -> dimension,
    "RawGeneratorCount" -> Length[generatorRows],
    "CommonDenominator" -> commonDenominator,
    "ColumnShifts" -> columnShifts,
    "PolynomialBasis" -> popov["Rows"],
    "RationalBasis" -> Together[popov["Rows"]/commonDenominator],
    "BasisInDerivativeGenerators" -> popov["GeneratorCoordinates"],
    "ModuleRank" -> Length[popov["Rows"]],
    "LeadingPositions" -> popov["LeadingPositions"],
    "RowDegrees" -> popov["RowDegrees"],
    "ReductionIterations" -> popov["Iterations"]
  |>
];

End[];
EndPackage[];
