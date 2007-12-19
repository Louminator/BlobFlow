/* This is a lookup table for asymmetric splitting parameters. */

/* alpha, rA, rB, gammaA, gammaB, aA2. */
/* Central element: gamma <- 1 - gammaA - gammaB */
/* Pair A: x <- x+-rA*s2, a2 <- a2*aA2, gamma <- gammaA/2 */
/* Pair B: x <- x+-rB*s2, gamma <- gammaB/2 */

static double split15asym_parms[40][6] =
  {
    { 6.0000000e-01,  1.7227257e+00,  1.6438763e+00,  3.3485127e-01,  3.2531355e-01,  5.2969246e-01},
    { 6.1000000e-01,  1.7320581e+00,  1.6396988e+00,  3.3092771e-01,  3.2366041e-01,  5.3535647e-01},
    { 6.2000000e-01,  1.7408716e+00,  1.6346415e+00,  3.2689258e-01,  3.2203268e-01,  5.4114224e-01},
    { 6.3000000e-01,  1.7491286e+00,  1.6286914e+00,  3.2274955e-01,  3.2043145e-01,  5.4707036e-01},
    { 6.4000000e-01,  1.7568122e+00,  1.6218622e+00,  3.1849895e-01,  3.1885111e-01,  5.5314393e-01},
    { 6.5000000e-01,  1.7639152e+00,  1.6140862e+00,  3.1414144e-01,  3.1730954e-01,  5.5937502e-01},
    { 6.6000000e-01,  1.7703722e+00,  1.6053944e+00,  3.0968598e-01,  3.1579449e-01,  5.6578254e-01},
    { 6.7000000e-01,  1.7761702e+00,  1.5957527e+00,  3.0512961e-01,  3.1431152e-01,  5.7237441e-01},
    { 6.8000000e-01,  1.7812652e+00,  1.5851400e+00,  3.0047809e-01,  3.1286593e-01,  5.7916537e-01},
    { 6.9000000e-01,  1.7856110e+00,  1.5735459e+00,  2.9573519e-01,  3.1145350e-01,  5.8617319e-01},
    { 7.0000000e-01,  1.7891559e+00,  1.5609177e+00,  2.9090388e-01,  3.1008491e-01,  5.9340349e-01},
    { 7.1000000e-01,  1.7918437e+00,  1.5472693e+00,  2.8598858e-01,  3.0875384e-01,  6.0088718e-01},
    { 7.2000000e-01,  1.7936221e+00,  1.5325575e+00,  2.8098959e-01,  3.0746980e-01,  6.0862669e-01},
    { 7.3000000e-01,  1.7944064e+00,  1.5167532e+00,  2.7591694e-01,  3.0622845e-01,  6.1665253e-01},
    { 7.4000000e-01,  1.7941240e+00,  1.4998257e+00,  2.7077351e-01,  3.0503434e-01,  6.2497285e-01},
    { 7.5000000e-01,  1.7926873e+00,  1.4817431e+00,  2.6556197e-01,  3.0388918e-01,  6.3360606e-01},
    { 7.6000000e-01,  1.7899914e+00,  1.4624738e+00,  2.6029027e-01,  3.0279375e-01,  6.4257871e-01},
    { 7.7000000e-01,  1.7859216e+00,  1.4419680e+00,  2.5496427e-01,  3.0175129e-01,  6.5191008e-01},
    { 7.8000000e-01,  1.7803589e+00,  1.4201949e+00,  2.4958825e-01,  3.0075950e-01,  6.6161891e-01},
    { 7.9000000e-01,  1.7731553e+00,  1.3970852e+00,  2.4417030e-01,  2.9982468e-01,  6.7172377e-01},
    { 8.0000000e-01,  1.7641394e+00,  1.3726041e+00,  2.3871837e-01,  2.9894255e-01,  6.8226181e-01},
    { 8.1000000e-01,  1.7531331e+00,  1.3466863e+00,  2.3323818e-01,  2.9811299e-01,  6.9323552e-01},
    { 8.2000000e-01,  1.7399258e+00,  1.3192475e+00,  2.2773907e-01,  2.9733970e-01,  7.0467847e-01},
    { 8.3000000e-01,  1.7242797e+00,  1.2902155e+00,  2.2222941e-01,  2.9662007e-01,  7.1661339e-01},
    { 8.4000000e-01,  1.7059098e+00,  1.2595026e+00,  2.1672050e-01,  2.9594838e-01,  7.2905797e-01},
    { 8.5000000e-01,  1.6845167e+00,  1.2269865e+00,  2.1121860e-01,  2.9532869e-01,  7.4202808e-01},
    { 8.6000000e-01,  1.6597291e+00,  1.1925477e+00,  2.0573746e-01,  2.9475302e-01,  7.5554863e-01},
    { 8.7000000e-01,  1.6311332e+00,  1.1560277e+00,  2.0028713e-01,  2.9421756e-01,  7.6963206e-01},
    { 8.8000000e-01,  1.5982417e+00,  1.1172344e+00,  1.9487926e-01,  2.9371836e-01,  7.8428784e-01},
    { 8.9000000e-01,  1.5604719e+00,  1.0759335e+00,  1.8952734e-01,  2.9325095e-01,  7.9952792e-01},
    { 9.0000000e-01,  1.5171556e+00,  1.0318504e+00,  1.8424179e-01,  2.9279791e-01,  8.1534462e-01},
    { 9.1000000e-01,  1.4674570e+00,  9.8460029e-01,  1.7903636e-01,  2.9235854e-01,  8.3173770e-01},
    { 9.2000000e-01,  1.4103626e+00,  9.3368211e-01,  1.7392285e-01,  2.9193211e-01,  8.4868816e-01},
    { 9.3000000e-01,  1.3445786e+00,  8.7848825e-01,  1.6891266e-01,  2.9149122e-01,  8.6618113e-01},
    { 9.4000000e-01,  1.2683937e+00,  8.1811334e-01,  1.6401790e-01,  2.9103200e-01,  8.8417910e-01},
    { 9.5000000e-01,  1.1794566e+00,  7.5125214e-01,  1.5924522e-01,  2.9055926e-01,  9.0263996e-01},
    { 9.6000000e-01,  1.0742484e+00,  6.7596591e-01,  1.5460404e-01,  2.9004616e-01,  9.2151675e-01},
    { 9.7000000e-01,  9.4701082e-01,  5.8895345e-01,  1.5010242e-01,  2.8948941e-01,  9.4075416e-01},
    { 9.8000000e-01,  7.8679376e-01,  4.8383188e-01,  1.4574295e-01,  2.8888364e-01,  9.6029155e-01},
    { 9.9000000e-01,  5.6587638e-01,  3.4424496e-01,  1.4152764e-01,  2.8822869e-01,  9.8006180e-01}
  };