// -----------------------------------------------------------------------
/// \brief Calcule un terme de contrainte a partir d'une homographie.
///
/// \param H: matrice 3*3 définissant l'homographie.
/// \param i: premiere colonne.
/// \param j: deuxieme colonne.
/// \return vecteur definissant le terme de contrainte.
// -----------------------------------------------------------------------
function v = ZhangConstraintTerm(H, i, j)
  // A modifier! 
  v = [H(1,i) * H(1,j); H(1,i) * H(2,j) + H(2,i) * H(1,j); H(2,i) * H(2,j); H(3,i) * H(1,j) + H(1,i) * H(3,j); H(3,i) * H(2,j) + H(2,i) * H(3,j); H(3,i) * H(3,j)]'
endfunction

// -----------------------------------------------------------------------
/// \brief Calcule deux equations de contrainte a partir d'une homographie
///
/// \param H: matrice 3*3 définissant l'homographie.
/// \return matrice 2*6 definissant les deux contraintes.
// -----------------------------------------------------------------------
function v = ZhangConstraints(H)
  v = [ZhangConstraintTerm(H, 1, 2); ...
    ZhangConstraintTerm(H, 1, 1) - ZhangConstraintTerm(H, 2, 2)];
endfunction

// -----------------------------------------------------------------------
/// \brief Calcule la matrice des parametres intrinseques.
///
/// \param b: vecteur resultant de l'optimisation de Zhang.
/// \return matrice 3*3 des parametres intrinseques.
// -----------------------------------------------------------------------
function A = IntrinsicMatrix(b)
  v0 = (b(2)*b(4) - b(1)*b(5))/(b(1)*b(3) - b(2) * b(2))
  L = b(6) - [b(4)*b(4) + v0 * (b(2)*b(4) - b(1) * b(5))] / b(1)
  a = sqrt(L/b(1))
  B = sqrt(L * b(1) / (b(1) * b(3) - b(2)*b(2)))
  G = - b(2) * a * a * B / L
  u0 = G * v0 / B - b(4) * a * a / L
  A = [a,G,u0;
       0,B,v0;
       0,0,1];
endfunction

// -----------------------------------------------------------------------
/// \brief Calcule la matrice des parametres extrinseques.
///
/// \param iA: inverse de la matrice intrinseque.
/// \param H: matrice 3*3 definissant l'homographie.
/// \return matrice 3*4 des parametres extrinseques.
// -----------------------------------------------------------------------
function E = ExtrinsicMatrix(iA, H)
  // A modifier!
  E = rand(3, 4);
endfunction

