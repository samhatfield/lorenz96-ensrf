*> \brief \b DGETRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRF2 computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>            
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = min(m,n)
*>        [  A21 | A22  ]       n2 = n-n1
*>            
*>                                       [ A11 ]
*> The subroutine calls itself to factor [ --- ],
*>                                       [ A12 ]
*>                 [ A12 ]
*> do the swaps on [ --- ], solve A12, update A22,
*>                 [ A22 ]
*>
*> then calls itself to factor A22 and do the swaps on A21.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE dgetrf2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, TEMP
      INTEGER            I, IINFO, N1, N2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IDAMAX
      EXTERNAL           dlamch, idamax
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemm, dscal, dlaswp, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGETRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN

      IF ( m.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         ipiv( 1 ) = 1
         IF ( a(1,1).EQ.zero )
     $      info = 1
*
      ELSE IF( n.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         sfmin = dlamch('S')
*
*        Find pivot and test for singularity
*
         i = idamax( m, a( 1, 1 ), 1 )
         ipiv( 1 ) = i
         IF( a( i, 1 ).NE.zero ) THEN
*
*           Apply the interchange
*
            IF( i.NE.1 ) THEN
               temp = a( 1, 1 )
               a( 1, 1 ) = a( i, 1 )
               a( i, 1 ) = temp
            END IF
*
*           Compute elements 2:M of the column
*
            IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
               CALL dscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
            ELSE
               DO 10 i = 1, m-1
                  a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            info = 1
         END IF
*
      ELSE
*
*        Use recursive code
*
         n1 = min( m, n ) / 2
         n2 = n-n1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL dgetrf2( m, n1, a, lda, ipiv, iinfo )

         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL dlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
*
*        Solve A12
*
         CALL dtrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda, 
     $               a( 1, n1+1 ), lda )
*
*        Update A22
*
         CALL dgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda, 
     $               a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
*
*        Factor A22
*
         CALL dgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),
     $                 iinfo )
*
*        Adjust INFO and the pivot indices
*
         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo + n1
         DO 20 i = n1+1, min( m, n )
            ipiv( i ) = ipiv( i ) + n1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL dlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
*
      END IF
      RETURN
*
*     End of DGETRF2
*
      END
