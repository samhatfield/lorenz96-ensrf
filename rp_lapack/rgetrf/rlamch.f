*> \brief \b RLAMCH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION RLAMCH( CMACH )
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> RLAMCH determines double precision machine parameters.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CMACH
*> \verbatim
*>          Specifies the value to be returned by RLAMCH:
*>          = 'E' or 'e',   RLAMCH := eps
*>          = 'S' or 's ,   RLAMCH := sfmin
*>          = 'B' or 'b',   RLAMCH := base
*>          = 'P' or 'p',   RLAMCH := eps*base
*>          = 'N' or 'n',   RLAMCH := t
*>          = 'R' or 'r',   RLAMCH := rnd
*>          = 'M' or 'm',   RLAMCH := emin
*>          = 'U' or 'u',   RLAMCH := rmin
*>          = 'L' or 'l',   RLAMCH := emax
*>          = 'O' or 'o',   RLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
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
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      FUNCTION RLAMCH( CMACH )

      USE RP_EMULATOR
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015

*     .. Return Type
      TYPE(RPE_VAR) RLMAMCH
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      TYPE(RPE_VAR)   ONE, ZERO
*     ..
*     .. Local Scalars ..
      TYPE(RPE_VAR)   RND, EPS, SFMIN, SMALL, RMACH, LARGE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intialise parameters
      ONE = 1.0d+0
      ZERO = 0.0d+0
      LARGE = 2.0d+0**1023.0d+0
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      RND = ONE
*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
*
*	  The large IF statement that used to be here has been removed, as
*     this function will only ever be called with the 'S' option
      SFMIN = TINY(ZERO)
	  ! Currently, the emulator returns infinity when calling HUGE(rpe_instance)
	  ! Until a cleverer solution is figured out, just use a
      ! conservative estimate for HUGE(rpe_var) (i.e. the max value for
      ! no significand bits)
      !SMALL = ONE / HUGE(ZERO)
	  SMALL = ONE / LARGE
      IF( SMALL.GE.SFMIN ) THEN
*
*        Use SMALL plus a bit, to avoid the possibility of rounding
*        causing overflow when computing  1/sfmin.`
*
         SFMIN = SMALL*( ONE+EPS )
      END IF
      RMACH = SFMIN

      RLAMCH = RMACH
      RETURN
*
*     End of RLAMCH
*
      END
************************************************************************
*> \brief \b DLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date November 2015
*> \ingroup auxOTHERauxiliary
*>
*> \param[in] A
*> \verbatim
*>          A is a DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is a DOUBLE PRECISION
*>          The values A and B.
*> \endverbatim
*>
      FUNCTION RLAMC3( A, B )

      USE RP_EMULATOR
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      TYPE(RPE_VAR)   A, B

      TYPE(RPE_VAR) RLAMC3
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************