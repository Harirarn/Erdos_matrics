(* The MaxTrace function computes the maximal trace of a matrix over a
given list of permutations *)
MaxTrace[ A_ , permlist_] := Max[Tr[A . #] & /@ permlist ];
                                    
(* The splTrace function acts as a heuristic filter to efficiently eliminate a substantial portion of non - Erdos matrices before performing
the computationally expensive defining condition check . It computes a specialized trace by iteratively selecting the largest entries in the
matrix while ensuring they can be placed as diagonal elements under some permutation. *)
splTrace[A_]:= Module[{B=A, u=Length[A], sum=0, pos},
    Do[
        pos = FirstPosition[B, Max[B]]; (* Find position of maximum entry *)
        sum = sum + Max[B]; (* Add the max entry to sum *)
        B = Drop[B , {pos [[1]]} ,{pos [[2]]}]; (* Remove the row and column of the selected entry *),
        u
    ];
    sum];

 (* The ClassifyPermutationEquivalence module classifies a given list of
matrices into equivalence classes under permutation equivalence using
the list of all permutation matrices *)
 ClassifyPermutationEquivalence [ matricesList_ , permlist_ ]:= Module [{
permsSquare , equivClasses , equivClassfun} ,
13 (* Compute the Cartesian product of the set of permutation matrices *)
14 permsSquare= Tuples [ permlist ,2];
15 (* Define a function to compute all matrices equivalent to a given
matrix under the action of permutation matrices *)
16 equivClassfun[ permMatrix_ ]:= Table [ permsSquare [[j ,1]]. permMatrix.
permsSquare [[j ,2]] ,{j ,1 , Length [ permsSquare ]}];
17 equivClasses= Reap[ NestWhile[ Complement[# , Sow[ equivClassfun[ First [#]]\[
Intersection]#]]& , matricesList ,#!={}&]][[2 ,1]]; (* NestWhile
iteratively extracts and removes equivalence classes , with Sow logging
them and Reap collecting at the end *)
18 equivClasses (* Returns the equivalence classes of matrices *)]
19
20 (* Supply values for n and k . These can be modified as needed . *)
21 n =4; k =6;
22 KUSHWAHA AND TRIPATHI
22
23 (* Step 1: Generate the set of permutation matrices of size nxn *)
24 perms = Permutations[ IdentityMatrix[n ]]; (* All nxn permutation matrices *)
25 permSet= ArrayReshape[ Rest[ perms ] ,{ n ! -1 ,n ^2}]; (* Excluding identity matrix
*)
26
27 (* Step 2: Generate subsets of permutation matrices containing the identity
matrix *)
28 permSubsets= Subsets[ permSet ,{k -1}]; (* Choose (k -1) permutations *)
29 numSubsets= Length [ permSubsets ]; (* Number of subsets: binomial(n ! -1 ,k -1) *)
30
31 Print[" Number of subsets of ", IntegerName[ k] ," ",n ," x " ,n ," permutation
matrices including identity: " , numSubsets ];
32
33 (* Include the identity matrix in each subset *)
34 subsetList= Transpose[ Join [{ Table[ Flatten[ IdentityMatrix[ n ]] , numSubsets]} ,
Transpose[ permSubsets ]]];
35
36 (* Step 3: Find linearly independent subsets of permutation matrices *)
37 ranks = Table [ MatrixRank[ Transpose[ subsetList[[ i ]]]] ,{i ,1 , numSubsets }];
38
39 Print[" Rank distribution of the spans of these subsets: ",Tally [ ranks ]];
40
41 numIndepSubsets= Count [ ranks ,k ]; (* Number of linearly independent subsets*)
42
43 Print[" The number of linearly independent subsets , |" ,
44 Subscript[\[ ScriptCapitalI] , k ] , "| = ", numIndepSubsets ]
45
46 (* Extract only the linearly independent subsets*)
47 indepSubsets= Select [ subsetList , MatrixRank[ Transpose[#]]== k &];
48
49 (* Step 4: Compute Gram matrices for all independent subsets*)
50 GramMatrix[ u_ ]:= Outer [Dot , indepSubsets [[ u ]] , indepSubsets [[ u ]] ,1];
51 gramList= GramMatrix/ @Range [ numIndepSubsets ];
52
53 (* Solve Subscript[M , X] y =1 to obtain coefficients y *)
54 solYList= LinearSolve[# , ConstantArray [1 , k ]]&/ @gramList;
55
56 (* Step 5: Select subsets X for which y satisfies y_i >0 *)
57 validYIndex= Select [ Range [1 , numIndepSubsets ] , Positive@Min@solYList [[#]]&];
58 numValidY= Length [ validYIndex ]; (* Count of valid subsets*)
59
60 Print[" The number of subsets for which y_i > 0 , |" ,Subscript[\[
ScriptCapitalG] , k] ,"| = ", numValidY];
61
62 (* Normalize coefficients to obtain stochastic weights x *)
63 solYValid= solYList [[ validYIndex ]];
64 weights =(#/ Total [#]) &/ @solYValid;
65
66 (* Extract corresponding subsets of permutation matrices*)
MARCUS–REE INEQUALITY AND ERDOS MATRICES 23 ˝
67 validSubsets= indepSubsets [[ validYIndex ]];
68
69 (* Compute the bistochastic matrix Subscript[A , X] for each subset *)
70 BistochasticMatrix [ u_ ]:= Total[ weights [[ u ]] validSubsets [[ u ]]]// N; (*N
converts the rational entries of matrices to decimal form , enabling
faster numerical computations in Mathematica. *)
71 candidates= Union[ Table[ ArrayReshape[ BistochasticMatrix [ i ] ,{n ,n }] ,{i ,1 ,
numValidY }]];
72
73 (* Step 6: Apply the heuristic splTrace filter *)
74 filteredCandidates = Select [ candidates , splTrace[#] <= Norm[# ," Frobenius" ]^2&];
75
76 Print[" The number of filtered candidate matrices after splTrace: " ,Length [
filteredCandidates ]];
77
78 (* Step 7: Verify the defining condition to find actual Erdos matrices*)
79 erdosMatrices= Select [ filteredCandidates , MaxTrace[# , perms ]== Norm [# ,"
Frobenius" ]^2&];
80
81 numErdos= Length [ erdosMatrices ];
82
83 Print[" The number of Erdos matrices , |", Subscript[\[ ScriptCapitalE ] , k] ," |
= ", numErdos];
84
85 (* Step 8: Identify non - equivalent Erdos matrices *)
86 eqClassesErdos= ClassifyPermutationEquivalence [ erdosMatrices , perms ];
87
88 (* Count of non - equivalent Erdos matrices *)
89 noneqcount= Length [ eqClassesErdos ];
90 Print[" The number of non - equivalent Erdos matrices , | ",Subscript[\[
ScriptCapitalH] , k] ,"| = ", noneqcount ];
91
92 (* Display count of equivalent matrices in each class *)
93 Print[" Sizes of equivalence classes: ",Table [ Length [ eqClassesErdos [[ j ]]] ,{
j ,1 , noneqcount }]];
94
95 (* List of representative Erdos matrices*)
96 noneqErdosList= Table [ eqClassesErdos [[j ,1]] ,{j ,1 , noneqcount }]// Rationalize;
97 (* Rationalize restores the decimal numbers to their fractional
representation *)
98
99 Print[" List of all " ,n ," x ",n ," Erdos matrices ( up to equivalence)
generated from subsets of " ,k ," permutation matrices: ", MatrixForm/
@noneqErdosList ];
100
101 (* Identifying A with A^T , all non equivalent nxn Erdos matrices generated
from subsets of size k are obtained by the modifying the equivClassfun
function of the ClassifyPermutationEquivalence module in the following
manner *)
24 KUSHWAHA AND TRIPATHI
102 IdentifyUptoTranspose [ matricesList_ , permlist_ ]:= Module [{ permsSquare ,
equivClasses , equivClassfun} ,
103 permsSquare= Tuples [ permlist ,2];
104 equivClassfun[ permMatrix_ ]:= Union [ Table [ permsSquare [[j ,1]]. permMatrix.
permsSquare [[j ,2]] ,{j ,1 , Length [ permsSquare ]}] , Table[ permsSquare [[j ,1]].
Transpose[ permMatrix]. permsSquare [[j ,2]] ,{j ,1 , Length [ permsSquare ]}]];
105 equivClasses= Reap[ NestWhile[ Complement[# , Sow[ equivClassfun[ First [#]]\[
Intersection]#]]& , matricesList ,#!={}&]][[2 ,1]]]
106
107 Print[" Non equivalent Erdos matrices identfied upto transpose: ",
MatrixForm/ @First / @IdentifyUptoTranspose [ noneqErdosList , perms ]];
108 Print[" Number of classes: " ,Length [ IdentifyUptoTranspose [ noneqErdosList ,
perms ]]]
