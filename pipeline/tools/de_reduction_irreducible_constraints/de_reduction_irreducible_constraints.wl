BeginPackage["DEIrreduciblePoleConstraints`"];

ComputeIrreduciblePoleConstraintStructures::usage =
  "ComputeIrreduciblePoleConstraintStructures[A, x, d, maxDI] computes " <>
  "joint reduction-pole numerator constraints at an irreducible factor " <>
  "d(x).  Local coefficients are represented over F_p[x]/(d), while all " <>
  "linear systems are expanded exactly over F_p.";

ProperNumeratorTupleContainedQ::usage =
  "ProperNumeratorTupleContainedQ[tuple, structure] tests a tuple " <>
  "(n_dI,...,n_1), with one proper numerator polynomial per master and " <>
  "pole layer, against a computed irreducible-pole structure.";

Options[ComputeIrreduciblePoleConstraintStructures] = {
  "Modulus" -> 2305843009213693951,
  "EquationOrderStart" -> Automatic,
  "EquationOrderStep" -> 2,
  "EquationOrderLimit" -> 16,
  "StableRepeats" -> 2,
  "FuchsianResidueOnly" -> False
};

Begin["`Private`"];

ClearAll[
  fieldScalar, polynomialOverField, monicPolynomial, factorValuation,
  nonzeroRowQ, rowBasis, seriesAdd, seriesScale, seriesMultiply,
  seriesInverse, seriesComposePolynomial, buildDeltaSeries,
  rationalLocalSeries, buildTaylorSystem, projectedJetBasis,
  projectedTaylorBasisFuchsian,
  multiplicationMatrixFor, expandKMatrix, buildNumeratorTransform,
  buildConstraintMatrix, polynomialCoordinates, removeFactorPower,
  localFieldGeneratorBasis, localFieldRowSpaceFp
];

fieldScalar[value_, modulus_Integer] := Module[
  {rational = Together[value], numerator, denominator},
  numerator = Mod[Numerator[rational], modulus];
  denominator = Mod[Denominator[rational], modulus];
  Mod[numerator PowerMod[denominator, -1, modulus], modulus]
];

polynomialOverField[polynomial_, variable_, modulus_Integer] := Module[
  {degree, coefficients},
  If[TrueQ[PossibleZeroQ[polynomial]], Return[0]];
  degree = Exponent[polynomial, variable];
  coefficients = Table[
    fieldScalar[Coefficient[polynomial, variable, power], modulus],
    {power, 0, degree}
  ];
  Sum[coefficients[[power + 1]] variable^power, {power, 0, degree}]
];

monicPolynomial[polynomial_, variable_, modulus_Integer] := Module[
  {fieldPolynomial, degree, leading},
  fieldPolynomial = PolynomialMod[
    polynomialOverField[polynomial, variable, modulus],
    modulus
  ];
  degree = Exponent[fieldPolynomial, variable];
  leading = Coefficient[fieldPolynomial, variable, degree];
  PolynomialMod[
    Expand[fieldPolynomial PowerMod[leading, -1, modulus]],
    modulus
  ]
];

factorValuation[
  polynomial_,
  factor_,
  variable_,
  modulus_Integer
] := Module[{current, quotient, remainder, valuation = 0},
  If[TrueQ[PossibleZeroQ[polynomial]], Return[Infinity]];
  current = polynomialOverField[polynomial, variable, modulus];
  While[True,
    {quotient, remainder} = PolynomialQuotientRemainder[
      current, factor, variable, Modulus -> modulus
    ];
    If[remainder =!= 0, Break[]];
    valuation++;
    current = quotient;
  ];
  valuation
];

removeFactorPower[
  polynomial_,
  factor_,
  power_Integer,
  variable_,
  modulus_Integer
] := Nest[
  PolynomialQuotient[
    #,
    factor,
    variable,
    Modulus -> modulus
  ] &,
  polynomialOverField[polynomial, variable, modulus],
  power
];

nonzeroRowQ[row_List, modulus_Integer] :=
  AnyTrue[row, Mod[#, modulus] =!= 0 &];

rowBasis[matrix_List, modulus_Integer] := If[
  matrix === {} || Length[matrix] === 0,
  {},
  Select[
    RowReduce[matrix, Modulus -> modulus],
    nonzeroRowQ[#, modulus] &
  ]
];

seriesAdd[left_List, right_List, reduce_] :=
  MapThread[reduce[#1 + #2] &, {left, right}];

seriesScale[series_List, scalar_, reduce_] :=
  (reduce[scalar #] & /@ series);

seriesMultiply[
  left_List,
  right_List,
  maximumOrder_Integer,
  reduce_
] := Table[
  reduce[
    Sum[
      left[[index + 1]] right[[power - index + 1]],
      {index, 0, power}
    ]
  ],
  {power, 0, maximumOrder}
];

seriesInverse[
  series_List,
  maximumOrder_Integer,
  reduce_,
  inverse_
] := Module[{result, inverseConstant},
  inverseConstant = inverse[First[series]];
  result = ConstantArray[0, maximumOrder + 1];
  result[[1]] = inverseConstant;
  Do[
    result[[power + 1]] = reduce[
      -inverseConstant Sum[
        series[[index + 1]] result[[power - index + 1]],
        {index, 1, power}
      ]
    ],
    {power, 1, maximumOrder}
  ];
  result
];

seriesComposePolynomial[
  polynomial_,
  argumentSeries_List,
  variable_,
  maximumOrder_Integer,
  modulus_Integer,
  reduce_
] := Module[{degree, coefficients, result, power},
  degree = If[
    TrueQ[PossibleZeroQ[polynomial]],
    0,
    Exponent[polynomial, variable]
  ];
  coefficients = Table[
    fieldScalar[
      Coefficient[polynomial, variable, index],
      modulus
    ],
    {index, 0, degree}
  ];
  result = ConstantArray[0, maximumOrder + 1];
  power = ConstantArray[0, maximumOrder + 1];
  power[[1]] = 1;
  Do[
    result = seriesAdd[
      result,
      seriesScale[power, coefficients[[index + 1]], reduce],
      reduce
    ];
    power = seriesMultiply[
      power, argumentSeries, maximumOrder, reduce
    ],
    {index, 0, degree}
  ];
  result
];

buildDeltaSeries[
  factor_,
  variable_,
  generator_,
  maximumOrder_Integer,
  modulus_Integer,
  reduce_,
  inverse_
] := Module[
  {
    derivative, result, target, factorSeries, derivativeSeries,
    inverseDerivative, correction, iterations
  },
  derivative = D[factor, variable];
  result = ConstantArray[0, maximumOrder + 1];
  result[[1]] = generator;
  target = ConstantArray[0, maximumOrder + 1];
  If[maximumOrder >= 1, target[[2]] = 1];
  iterations = Ceiling[Log[2, maximumOrder + 1]] + 2;
  Do[
    factorSeries = seriesComposePolynomial[
      factor, result, variable, maximumOrder, modulus, reduce
    ];
    factorSeries = seriesAdd[
      factorSeries,
      seriesScale[target, -1, reduce],
      reduce
    ];
    derivativeSeries = seriesComposePolynomial[
      derivative, result, variable, maximumOrder, modulus, reduce
    ];
    inverseDerivative = seriesInverse[
      derivativeSeries, maximumOrder, reduce, inverse
    ];
    correction = seriesMultiply[
      factorSeries, inverseDerivative, maximumOrder, reduce
    ];
    result = seriesAdd[
      result,
      seriesScale[correction, -1, reduce],
      reduce
    ],
    {iterations}
  ];
  result
];

rationalLocalSeries[
  expression_,
  variable_,
  deltaSeries_List,
  minimumPower_Integer,
  maximumPower_Integer,
  modulus_Integer,
  reduce_,
  inverse_
] := Module[
  {
    together, numerator, denominator, seriesOrder,
    numeratorSeries, denominatorSeries, numeratorValuation,
    denominatorValuation, shift, numeratorRegular, denominatorRegular,
    maximumRegularOrder, quotient, inverseConstant
  },
  If[TrueQ[PossibleZeroQ[expression]],
    Return[Association@Table[power -> 0,
      {power, minimumPower, maximumPower}]]
  ];
  together = Together[expression];
  numerator = polynomialOverField[
    Numerator[together], variable, modulus
  ];
  denominator = polynomialOverField[
    Denominator[together], variable, modulus
  ];
  seriesOrder = Length[deltaSeries] - 1;
  numeratorSeries = seriesComposePolynomial[
    numerator, deltaSeries, variable, seriesOrder, modulus, reduce
  ];
  denominatorSeries = seriesComposePolynomial[
    denominator, deltaSeries, variable, seriesOrder, modulus, reduce
  ];
  numeratorValuation = FirstPosition[
    numeratorSeries, Except[0], {Infinity},
    {1}, Heads -> False
  ][[1]] - 1;
  denominatorValuation = FirstPosition[
    denominatorSeries, Except[0], {Infinity},
    {1}, Heads -> False
  ][[1]] - 1;
  shift = numeratorValuation - denominatorValuation;
  numeratorRegular = Drop[numeratorSeries, numeratorValuation];
  denominatorRegular = Drop[denominatorSeries, denominatorValuation];
  maximumRegularOrder = Max[0, maximumPower - shift];
  numeratorRegular = PadRight[
    numeratorRegular, maximumRegularOrder + 1
  ];
  denominatorRegular = PadRight[
    denominatorRegular, maximumRegularOrder + 1
  ];
  quotient = ConstantArray[0, maximumRegularOrder + 1];
  inverseConstant = inverse[First[denominatorRegular]];
  Do[
    quotient[[power + 1]] = reduce[
      inverseConstant (
        numeratorRegular[[power + 1]] -
        Sum[
          denominatorRegular[[index + 1]]
            quotient[[power - index + 1]],
          {index, 1, power}
        ]
      )
    ],
    {power, 0, maximumRegularOrder}
  ];
  Association@Table[
    power -> If[
      0 <= power - shift <= maximumRegularOrder,
      quotient[[power - shift + 1]],
      0
    ],
    {power, minimumPower, maximumPower}
  ]
];

buildTaylorSystem[
  coefficients_Association,
  expandedDimension_Integer,
  poleOrder_Integer,
  equationOrder_Integer,
  modulus_Integer
] := Module[
  {identity, zero, maximumJet, coefficient, blockRows},
  identity = SparseArray[IdentityMatrix[expandedDimension]];
  zero = SparseArray[{}, {expandedDimension, expandedDimension}];
  maximumJet = equationOrder + Max[poleOrder, 1];
  coefficient[index_] := Lookup[coefficients, index, zero];
  blockRows = Table[
    ArrayFlatten[{
      Table[
        SparseArray@Mod[
          coefficient[power - jet] -
            If[
              power >= 0 && jet === power + 1,
              (power + 1) identity,
              zero
            ],
          modulus
        ],
        {jet, 0, maximumJet}
      ]
    }],
    {power, -poleOrder, equationOrder}
  ];
  <|
    "MaximumJetOrder" -> maximumJet,
    "Matrix" -> SparseArray[Join @@ blockRows]
  |>
];

projectedJetBasis[
  systemMatrix_,
  retainedColumns_Integer,
  modulus_Integer
] := Module[{null, projected},
  null = NullSpace[systemMatrix, Modulus -> modulus];
  If[null === {}, Return[{}]];
  projected = null[[All, 1 ;; retainedColumns]];
  rowBasis[projected, modulus]
];

projectedTaylorBasisFuchsian[
  coefficients_Association,
  expandedDimension_Integer,
  retainedJetCount_Integer,
  equationOrder_Integer,
  modulus_Integer
] := Module[
  {
    identity, zero, residue, parameterization, freeDimension,
    oldCoefficient, newCoefficient, equationMatrix, nullColumns,
    embedding, projected
  },
  identity = IdentityMatrix[expandedDimension, SparseArray];
  zero = SparseArray[{}, {expandedDimension, expandedDimension}];
  residue = Lookup[coefficients, -1, zero];
  parameterization = Transpose[
    NullSpace[residue, Modulus -> modulus]
  ];
  If[
    parameterization === {} ||
      Dimensions[parameterization][[2]] === 0,
    Return[{}]
  ];
  Do[
    freeDimension = Dimensions[parameterization][[2]];
    oldCoefficient = ArrayFlatten[{
      Table[
        Lookup[coefficients, power - jet, zero],
        {jet, 0, power}
      ]
    }];
    newCoefficient = SparseArray@Mod[
      residue - (power + 1) identity,
      modulus
    ];
    equationMatrix = ArrayFlatten[{{
      SparseArray@Mod[
        oldCoefficient . parameterization,
        modulus
      ],
      newCoefficient
    }}];
    nullColumns = Transpose[
      NullSpace[equationMatrix, Modulus -> modulus]
    ];
    If[
      nullColumns === {} ||
        Dimensions[nullColumns][[2]] === 0,
      Return[{}]
    ];
    embedding = ArrayFlatten[{
      {
        parameterization,
        SparseArray[
          {},
          {
            Dimensions[parameterization][[1]],
            expandedDimension
          }
        ]
      },
      {
        SparseArray[
          {},
          {expandedDimension, freeDimension}
        ],
        identity
      }
    }];
    parameterization = SparseArray@Mod[
      embedding . nullColumns,
      modulus
    ],
    {power, 0, equationOrder}
  ];
  projected = Transpose[
    parameterization[[
      1 ;; retainedJetCount expandedDimension,
      All
    ]]
  ];
  rowBasis[Normal[projected], modulus]
];

multiplicationMatrixFor[
  element_,
  generator_,
  degree_Integer,
  reduce_
] := Transpose@Table[
  PadRight[
    CoefficientList[
      reduce[element generator^column],
      generator
    ],
    degree
  ],
  {column, 0, degree - 1}
];

expandKMatrix[
  matrix_,
  generator_,
  degree_Integer,
  reduce_
] := SparseArray@ArrayFlatten[
  Map[
    multiplicationMatrixFor[#, generator, degree, reduce] &,
    matrix,
    {2}
  ]
];

localFieldGeneratorBasis[
  rows_List,
  dimension_Integer,
  generator_,
  degree_Integer,
  modulus_Integer,
  reduce_
] := Module[
  {
    generatorBlock, selected = {}, span = {}, candidate,
    orbit, enlarged
  },
  If[rows === {}, Return[{}]];
  generatorBlock = KroneckerProduct[
    IdentityMatrix[dimension],
    multiplicationMatrixFor[
      generator, generator, degree, reduce
    ]
  ];
  Do[
    enlarged = rowBasis[Append[span, candidate], modulus];
    If[Length[enlarged] > Length[span],
      AppendTo[selected, candidate];
      orbit = NestList[
        Mod[# . Transpose[generatorBlock], modulus] &,
        candidate,
        degree - 1
      ];
      span = rowBasis[Join[span, orbit], modulus]
    ],
    {candidate, rows}
  ];
  selected
];

localFieldRowSpaceFp[
  matrix_,
  generator_,
  degree_Integer,
  modulus_Integer,
  reduce_
] := Module[{dimension, generatorBlock, coordinateRows, orbits},
  dimension = Length[First[matrix]];
  generatorBlock = KroneckerProduct[
    IdentityMatrix[dimension],
    multiplicationMatrixFor[
      generator, generator, degree, reduce
    ]
  ];
  coordinateRows = Flatten[
    PadRight[
      CoefficientList[reduce[#], generator],
      degree
    ] & /@ #,
    1
  ] & /@ matrix;
  orbits = Flatten[
    NestList[
      Mod[# . Transpose[generatorBlock], modulus] &,
      #,
      degree - 1
    ] & /@ coordinateRows,
    1
  ];
  rowBasis[orbits, modulus]
];

buildNumeratorTransform[
  deltaSeries_List,
  dimension_Integer,
  degree_Integer,
  dI_Integer,
  generator_,
  modulus_Integer,
  reduce_
] := Module[
  {
    fieldLayerBlocks, basisSeries, composed, localPower, rawPower,
    zeroFieldBlock, identityMasters, masterBlock
  },
  zeroFieldBlock = ConstantArray[0, {degree, degree}];
  identityMasters = IdentityMatrix[dimension];
  fieldLayerBlocks = Table[
    localPower = dI - localBlock + 1;
    rawPower = dI - rawBlock + 1;
    If[rawPower < localPower,
      zeroFieldBlock,
      Transpose@Table[
        basisSeries = ConstantArray[0, Length[deltaSeries]];
        basisSeries = seriesComposePolynomial[
          generator^basisPower,
          deltaSeries,
          generator,
          Length[deltaSeries] - 1,
          modulus,
          reduce
        ];
        PadRight[
          CoefficientList[
            reduce[basisSeries[[rawPower - localPower + 1]]],
            generator
          ],
          degree
        ],
        {basisPower, 0, degree - 1}
      ]
    ],
    {localBlock, 1, dI},
    {rawBlock, 1, dI}
  ];
  ArrayFlatten@Table[
    masterBlock = fieldLayerBlocks[[localBlock, rawBlock]];
    KroneckerProduct[identityMasters, masterBlock],
    {localBlock, 1, dI},
    {rawBlock, 1, dI}
  ]
];

buildConstraintMatrix[
  jetBasisRows_List,
  dimension_Integer,
  degree_Integer,
  dI_Integer,
  modulus_Integer,
  generator_,
  reduce_
] := Module[
  {
    expandedDimension = dimension degree, zeroMasterBlock,
    rows, constraintBlocks, solution, jetValue, outputBlocks,
    localPower, polePower,
    jetOrder, masterBlocks
  },
  If[jetBasisRows === {}, Return[
    ConstantArray[0, {0, dI expandedDimension}]
  ]];
  zeroMasterBlock =
    ConstantArray[0, {degree, expandedDimension}];
  constraintBlocks = Flatten[Table[
    solution = jetBasisRows[[solutionIndex]];
    outputBlocks = Table[
      localPower = dI - localBlock + 1;
      If[localPower < polePower,
        zeroMasterBlock,
        jetOrder = localPower - polePower;
        masterBlocks = Table[
          jetValue = Sum[
            solution[[
              jetOrder expandedDimension +
                (masterIndex - 1) degree + coordinate
            ]] generator^(coordinate - 1),
            {coordinate, 1, degree}
          ];
          multiplicationMatrixFor[
            reduce[jetValue], generator, degree, reduce
          ],
          {masterIndex, 1, dimension}
        ];
        ArrayFlatten[{masterBlocks}]
      ],
      {localBlock, 1, dI}
    ];
    ArrayFlatten[{outputBlocks}],
    {solutionIndex, 1, Length[jetBasisRows]},
    {polePower, 1, dI}
  ], 1];
  rows = Join @@ constraintBlocks;
  Mod[rows, modulus]
];

polynomialCoordinates[
  polynomial_,
  variable_,
  generator_,
  factor_,
  degree_Integer,
  modulus_Integer
] := PadRight[
  CoefficientList[
    PolynomialRemainder[
      polynomialOverField[polynomial, variable, modulus] /.
        variable -> generator,
      factor /. variable -> generator,
      generator,
      Modulus -> modulus
    ],
    generator
  ],
  degree
];

ComputeIrreduciblePoleConstraintStructures[
  matrix_,
  variable_Symbol,
  inputFactor_,
  maximumDI_Integer?Positive,
  OptionsPattern[]
] := Module[
  {
    modulus, equationStart, equationStep, equationLimit, stableRepeats,
    residueOnly,
    factor, degree, generator, factorGenerator, reduce, inverse,
    dimension, expandedDimension, entryData, poleOrder, seriesOrder,
    deltaSeries, deltaCheck, deltaDerivative, matrixSeries,
    coefficientKMatrices, coefficientFpMatrices, equationOrder,
    system, retainedColumns, jetBasis, previousBasis = Missing[],
    stableCount = 0, stabilizedQ = False, history = {}, structures,
    dI, constraintLocal, numeratorTransform, constraintRaw,
    constraintRank, allowedBasis, ambientFp, allowedFp, allowedK,
    numeratorTransformRank, residueK, residueFp, residueJetBasis,
    residueJetGenerators,
    residueConstraintLocal, residueNumeratorTransform,
    residueConstraintRaw, residueConstraintRank, residueAllowedBasis,
    residueAmbientFp, residueAllowedFp, residueAllowedK,
    residueTransformRank, resonantOrders
  },
  modulus = OptionValue["Modulus"];
  If[!IntegerQ[modulus] || modulus <= 1,
    Print["A positive prime Modulus is required."];
    Return[$Failed]
  ];
  equationStep = OptionValue["EquationOrderStep"];
  equationLimit = OptionValue["EquationOrderLimit"];
  stableRepeats = OptionValue["StableRepeats"];
  residueOnly = TrueQ[OptionValue["FuchsianResidueOnly"]];
  factor = monicPolynomial[inputFactor, variable, modulus];
  degree = Exponent[factor, variable];
  generator = Unique["alpha"];
  factorGenerator = factor /. variable -> generator;
  reduce[expression_] := PolynomialRemainder[
    PolynomialMod[Expand[expression], modulus],
    factorGenerator,
    generator,
    Modulus -> modulus
  ];
  inverse[element_] := Module[{gcd},
    gcd = PolynomialExtendedGCD[
      reduce[element],
      factorGenerator,
      generator,
      Modulus -> modulus
    ];
    If[First[gcd] =!= 1,
      Print["Noninvertible local-field element: ", InputForm[element]];
      Return[$Failed]
    ];
    reduce[gcd[[2, 1]]]
  ];
  dimension = Length[matrix];
  expandedDimension = dimension degree;
  entryData = Map[
    Function[entry,
      If[entry === 0,
        <|"NumeratorValuation" -> Infinity,
          "DenominatorValuation" -> 0|>,
        With[{together = Together[entry]},
          <|
            "NumeratorValuation" -> factorValuation[
              Numerator[together], factor, variable, modulus
            ],
            "DenominatorValuation" -> factorValuation[
              Denominator[together], factor, variable, modulus
            ]
          |>
        ]
      ]
    ],
    matrix,
    {2}
  ];
  poleOrder = Max[
    0,
    Max[
      Flatten[
        Map[
          If[
            #["NumeratorValuation"] === Infinity,
            0,
            #["DenominatorValuation"] -
              #["NumeratorValuation"]
          ] &,
          entryData,
          {2}
        ]
      ]
    ]
  ];
  If[residueOnly && poleOrder === 1 && maximumDI === 1,
    residueK = MapThread[
      Function[{entry, valuation},
        If[
          entry === 0 ||
            valuation["NumeratorValuation"] === Infinity ||
            valuation["DenominatorValuation"] -
              valuation["NumeratorValuation"] =!= 1,
          0,
          With[
            {
              together = Together[entry],
              numeratorPower = valuation["NumeratorValuation"],
              denominatorPower = valuation["DenominatorValuation"]
            },
            reduce[
              removeFactorPower[
                Numerator[together],
                factor,
                numeratorPower,
                variable,
                modulus
              ] /. variable -> generator
              inverse[
                reduce[
                  removeFactorPower[
                    Denominator[together],
                    factor,
                    denominatorPower,
                    variable,
                    modulus
                  ] D[factor, variable] /. variable -> generator
                ]
              ]
            ]
          ]
        ]
      ],
      {matrix, entryData},
      2
    ];
    residueFp = expandKMatrix[
      residueK, generator, degree, reduce
    ];
    residueAllowedBasis = localFieldRowSpaceFp[
      residueK,
      generator,
      degree,
      modulus,
      reduce
    ];
    residueAllowedFp = Length[residueAllowedBasis];
    residueAllowedK = residueAllowedFp/degree;
    residueConstraintLocal = rowBasis[
      NullSpace[
        residueAllowedBasis,
        Modulus -> modulus
      ],
      modulus
    ];
    residueNumeratorTransform = IdentityMatrix[
      expandedDimension
    ];
    residueConstraintRaw = residueConstraintLocal;
    residueConstraintRank = Length[residueConstraintRaw];
    residueAmbientFp = expandedDimension;
    residueTransformRank = expandedDimension;
    resonantOrders = Select[
      Range[equationLimit],
      MatrixRank[
        Mod[
          residueFp -
            # IdentityMatrix[expandedDimension],
          modulus
        ],
        Modulus -> modulus
      ] < expandedDimension &
    ];
    Return[<|
      "Variable" -> variable,
      "Factor" -> factor,
      "FactorCoefficients" ->
        CoefficientList[factor, variable],
      "Modulus" -> modulus,
      "Dimension" -> dimension,
      "LocalFieldDegree" -> degree,
      "ExpandedDimensionFp" -> expandedDimension,
      "DEPoleOrder" -> poleOrder,
      "MaximumDI" -> maximumDI,
      "StabilizedQ" -> resonantOrders === {},
      "ResidueOnlyQ" -> True,
      "PositiveIntegerResonanceScanLimit" -> equationLimit,
      "PositiveIntegerResonantOrdersFp" -> resonantOrders,
      "StabilityHistory" -> {
        <|
          "EquationOrder" -> 0,
          "MaximumJetOrder" -> 0,
          "RetainedJetDimensionFp" ->
            residueAmbientFp - residueAllowedFp,
          "RetainedJetDimensionK" ->
            dimension - residueAllowedK
        |>
      },
      "DeltaSeriesCoefficients" -> {generator},
      "TaylorJetBasisRowsFp" -> Missing[
        "NotStoredInResidueRowSpaceMode"
      ],
      "Structures" -> <|
        1 -> <|
          "dI" -> 1,
          "Variable" -> variable,
          "Factor" -> factor,
          "Modulus" -> modulus,
          "Dimension" -> dimension,
          "BaseField" -> "F_p",
          "LocalFieldDegree" -> degree,
          "AmbientDimensionFp" -> residueAmbientFp,
          "AmbientDimensionK" -> dimension,
          "ConstraintRankFp" -> residueConstraintRank,
          "ConstraintRankK" -> residueConstraintRank/degree,
          "AllowedDimensionFp" -> residueAllowedFp,
          "AllowedDimensionK" -> residueAllowedK,
          "LocalNumeratorConstraintMatrixFp" ->
            residueConstraintLocal,
          "ProperNumeratorTransformationMatrixFp" ->
            residueNumeratorTransform,
          "ProperNumeratorTransformationRankFp" ->
            residueTransformRank,
          "ProperNumeratorTransformationInvertibleQ" -> True,
          "ProperNumeratorConstraintMatrixFp" ->
            residueConstraintRaw,
          "AllowedProperNumeratorBasisFp" ->
            residueAllowedBasis,
          "DimensionDivisibilityCheckQ" ->
            Divisible[residueConstraintRank, degree] &&
              Divisible[residueAllowedFp, degree]
        |>
      |>
    |>]
  ];
  equationStart = Replace[
    OptionValue["EquationOrderStart"],
    Automatic :>
      Max[0, maximumDI - 1 - Max[poleOrder, 1]]
  ];
  seriesOrder =
    equationLimit + poleOrder + maximumDI + 4;
  deltaSeries = buildDeltaSeries[
    factor, variable, generator, seriesOrder, modulus,
    reduce, inverse
  ];
  deltaCheck = seriesComposePolynomial[
    factor, deltaSeries, variable, seriesOrder, modulus, reduce
  ];
  If[
    deltaCheck =!=
      Join[If[seriesOrder >= 1, {0, 1}, {0}],
        ConstantArray[0, Max[0, seriesOrder - 1]]],
    Print["Failed to construct d(delta(t))=t."];
    Return[$Failed]
  ];
  deltaDerivative = Table[
    reduce[(power + 1) deltaSeries[[power + 2]]],
    {power, 0, seriesOrder - 1}
  ];
  deltaDerivative = PadRight[deltaDerivative, seriesOrder + 1];
  matrixSeries = Map[
    rationalLocalSeries[
      #, variable, deltaSeries, -poleOrder, equationLimit,
      modulus, reduce, inverse
    ] &,
    matrix,
    {2}
  ];
  coefficientKMatrices = Association@Table[
    power -> Map[Lookup[#, power, 0] &, matrixSeries, {2}],
    {power, -poleOrder, equationLimit}
  ];
  (* Convert dM/d(variable)=A M to dM/dt=(d variable/dt) A M. *)
  coefficientKMatrices = Association@Table[
    power -> Map[
      reduce[#] &,
      Sum[
        Lookup[
          coefficientKMatrices,
          power - derivativePower,
          ConstantArray[0, {dimension, dimension}]
        ] deltaDerivative[[derivativePower + 1]],
        {derivativePower, 0, Min[seriesOrder, power + poleOrder]}
      ],
      {2}
    ],
    {power, -poleOrder, equationLimit}
  ];
  coefficientFpMatrices = Association@KeyValueMap[
    #1 -> expandKMatrix[#2, generator, degree, reduce] &,
    coefficientKMatrices
  ];
  retainedColumns = maximumDI expandedDimension;
  jetBasis = {};
  Do[
    If[
      poleOrder === 1,
      jetBasis = projectedTaylorBasisFuchsian[
        coefficientFpMatrices,
        expandedDimension,
        maximumDI,
        equationOrder,
        modulus
      ];
      system = <|
        "MaximumJetOrder" -> equationOrder + 1
      |>,
      system = buildTaylorSystem[
        coefficientFpMatrices,
        expandedDimension,
        poleOrder,
        equationOrder,
        modulus
      ];
      jetBasis = projectedJetBasis[
        system["Matrix"], retainedColumns, modulus
      ]
    ];
    AppendTo[
      history,
      <|
        "EquationOrder" -> equationOrder,
        "MaximumJetOrder" -> system["MaximumJetOrder"],
        "RetainedJetDimensionFp" -> Length[jetBasis],
        "RetainedJetDimensionK" -> Length[jetBasis]/degree
      |>
    ];
    If[
      !MissingQ[previousBasis] && jetBasis === previousBasis,
      stableCount++,
      stableCount = 0
    ];
    If[stableCount >= stableRepeats,
      stabilizedQ = True;
      Break[]
    ];
    previousBasis = jetBasis,
    {equationOrder, equationStart, equationLimit, equationStep}
  ];
  structures = Association@Table[
    constraintLocal = buildConstraintMatrix[
      jetBasis, dimension, degree, dI, modulus, generator, reduce
    ];
    numeratorTransform = buildNumeratorTransform[
      deltaSeries, dimension, degree, dI, generator, modulus, reduce
    ];
    constraintRaw = Mod[
      constraintLocal . numeratorTransform,
      modulus
    ];
    constraintRaw = rowBasis[constraintRaw, modulus];
    constraintRank = Length[constraintRaw];
    ambientFp = dI expandedDimension;
    numeratorTransformRank = MatrixRank[
      numeratorTransform,
      Modulus -> modulus
    ];
    allowedFp = ambientFp - constraintRank;
    allowedK = allowedFp/degree;
    allowedBasis = NullSpace[constraintRaw, Modulus -> modulus];
    dI -> <|
      "dI" -> dI,
      "Variable" -> variable,
      "Factor" -> factor,
      "Modulus" -> modulus,
      "Dimension" -> dimension,
      "BaseField" -> "F_p",
      "LocalFieldDegree" -> degree,
      "AmbientDimensionFp" -> ambientFp,
      "AmbientDimensionK" -> dI dimension,
      "ConstraintRankFp" -> constraintRank,
      "ConstraintRankK" -> constraintRank/degree,
      "AllowedDimensionFp" -> allowedFp,
      "AllowedDimensionK" -> allowedK,
      "LocalNumeratorConstraintMatrixFp" -> constraintLocal,
      "ProperNumeratorTransformationMatrixFp" -> numeratorTransform,
      "ProperNumeratorTransformationRankFp" ->
        numeratorTransformRank,
      "ProperNumeratorTransformationInvertibleQ" ->
        numeratorTransformRank === ambientFp,
      "ProperNumeratorConstraintMatrixFp" -> constraintRaw,
      "AllowedProperNumeratorBasisFp" -> allowedBasis,
      "DimensionDivisibilityCheckQ" ->
        Divisible[constraintRank, degree] &&
          Divisible[allowedFp, degree]
    |>,
    {dI, 1, maximumDI}
  ];
  <|
    "Variable" -> variable,
    "Factor" -> factor,
    "FactorCoefficients" -> CoefficientList[factor, variable],
    "Modulus" -> modulus,
    "Dimension" -> dimension,
    "LocalFieldDegree" -> degree,
    "ExpandedDimensionFp" -> expandedDimension,
    "DEPoleOrder" -> poleOrder,
    "MaximumDI" -> maximumDI,
    "StabilizedQ" -> stabilizedQ,
    "StabilityHistory" -> history,
    "DeltaSeriesCoefficients" -> deltaSeries,
    "TaylorJetBasisRowsFp" -> jetBasis,
    "Structures" -> structures
  |>
];

ProperNumeratorTupleContainedQ[
  tuple_List,
  structure_Association
] := Module[
  {
    variable, factor, modulus, degree, generator, dimension, dI,
    vector, residual
  },
  variable = Lookup[structure, "Variable", Missing["Variable"]];
  factor = Lookup[structure, "Factor", Missing["Factor"]];
  modulus = structure["Modulus"];
  degree = structure["LocalFieldDegree"];
  dimension = structure["Dimension"];
  dI = structure["dI"];
  If[MissingQ[variable] || MissingQ[factor],
    Print["The structure does not contain Variable/Factor metadata."];
    Return[False]
  ];
  generator = Unique["alpha"];
  If[Dimensions[tuple] =!= {dI, dimension},
    Print[
      "Expected tuple dimensions ", {dI, dimension},
      ", got ", Dimensions[tuple]
    ];
    Return[False]
  ];
  vector = Flatten[
    polynomialCoordinates[
      #, variable, generator, factor, degree, modulus
    ] & /@ Flatten[tuple]
  ];
  residual = Mod[
    structure["ProperNumeratorConstraintMatrixFp"] . vector,
    modulus
  ];
  And @@ (# === 0 & /@ residual)
];

End[];
EndPackage[];
