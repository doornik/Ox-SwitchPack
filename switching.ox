#include <oxstd.oxh>
#include <oxprob.oxh>
#include <oxfloat.oxh>
#include "switching.oxh"
#import <maximize>
#import "maxscalar"

///////////////////////////////////////////////////////////////////////
//

/**    Constructor of the Switching class for function maximization using the alternating variables approach.
*/
Switching::Switching()
{
	m_mMaxTrace = <>;
	m_iResult = MAX_NOCONV;
	m_iLineSearchMode = LS_1STEP;
	m_cLineSearchWarmUp = 0;
	m_dLineSearchCrit = 1e-4;
	m_dLineSearchTol = 0.5;
	m_cExtraUpdate = 0;
	m_cIterCount = m_cUpdateCount = m_cLogLikCount = 0;
}

/**    Returns the version number as a string.
*/
Switching::GetVersion()
{
	return "2.0";
}
/**    Returns the maximization trace.
*/
Switching::GetTrace()
{
	return m_mMaxTrace;
}
/**    Returns the maximization statistics.
@returns row vector with maximization result code, no iterations, no of updates, no of loglik calls
*/
Switching::GetIterInfo()
{
	return m_iResult ~ m_cIterCount ~ m_cUpdateCount ~ m_cLogLikCount;
}
/**    Returns the name of the specified line search method.
*/
Switching::GetLineSearchName(const iLineSearchMode)
{
	switch_single (iLineSearchMode)
	{
		case LS_NONE:   	return "None";
		case LS_1STEP:  	return "L1Step";
		case LS_1STEP1:  	return "L1Step*";
		case LS_2STEP:  	return "L2Step";
		case LS_2STEPSQ3:	return "L1Step+SqS3";
		case LS_1STD:   	return "L1Std";
		case LS_SQS3:		return "LSqS3";
		case LS_SQS3G:		return "LgSqS3";
		case LS_PARABOLIC:	return "LParabolic";
		case LS_PARABOLIC2:	return "LParabolic2";
		case LS_POWELL:		return "LPowell";
		case LS_BRENT:		return "LBrent";
		case LS_BRENTSTD:	return "LBrentStd";
		case LS_QSTEP:		return "LQStep";
		case LS_USER1:		return "LUser1";
		case LS_USER2:		return "LUser2";
		case LS_USER3:		return "LUser3";
	}
	return "??";
}
/**    Returns the line search criterion.
*/
Switching::GetLineSearchCrit()
{
	return m_dLineSearchCrit;
}
/**    Returns the line search tolerance.
*/
Switching::GetLineSearchTol()
{
	return m_dLineSearchTol;
}
/**    Sets the line search criterion: the line search is entered if relative progress is
below this value (0: always enter, but see `Switching::SetLineSearchWarmUp').
*/
Switching::SetLineSearchCrit(const dLineSearchCrit)
{
	m_dLineSearchCrit = dLineSearchCrit;
}
/**    Sets the line search tolerance: only used for Powell and Brent.
*/
Switching::SetLineSearchTol(const dLineSearchTol)
{
	m_dLineSearchTol = dLineSearchTol;
}
/**    Sets the number of iteration of warm up: the line search is skipped during the warm up.
*/
Switching::SetLineSearchWarmUp(const cLineSearchWarmUp)
{
	m_cLineSearchWarmUp = cLineSearchWarmUp;
}
/**    Returns the current line search mode.
*/
Switching::GetLineSearchMode()
{
	return m_iLineSearchMode;
}
/**    Sets the line search mode: one of the Switching::LS_... constants.
*/
Switching::SetLineSearchMode(const iLineSearchMode)
{
	m_iLineSearchMode = iLineSearchMode;
}
/**    Can be used to insert an additional update.
*/
Switching::SetExtraUpdate(const cExtraUpdate)
{
	m_cExtraUpdate = cExtraUpdate;
}
/**    Can be used to switch on the intermediate trace of maximization progress (on a log10 scale,
otherwise it is only done at the beginning and end).
*/
static Switching::SetTraceProgress(const bTrace)
{
	sm_bTraceProgress = bTrace;
}

/**    Checks for convergence, see eqns (10) and (11).
@param dFunc: double, the new value of the objective function
@param dFuncPrev: double, the previous value of the objective function
@param mPi: matrix, the new value of coefficients
@param mPiPrev: matrix, the previous value of coefficients
@param dCrit: double, convergence tolerance \epsilon_1
@param iMaxIt: int, maximum number of iterations
@param iInfo: int, <=0: don't print anything, 1: print every iteration, >1: print iteration on a log10 scale
@param iItno: int, iteration number
@param mIterInfo: matrix, previous iteration info as maintained by this function
@param dStep: double, step length used in new coefficients
@param bForceTerminate: TRUE: return termination code
@returns array with 5 values: convergence criterion of objective function, convergence criterion of
coefficients, TRUE if finished (FALSE otherwise), return code
(MAX_ITERATING, MAX_CONV, MAX_WEAK_CONV, MAX_MAXIT, MAX_NOCONV, MAX_FUNC_FAIL), updated iteration info};
@comment The coefficients do not have to be the actual parameters over which the objective is
maximized. It can be a derived coefficient instead, e.g. one that is identified (identification is preferred,
otherwise the coefficient can change but the objective stay the same).
*/
static Switching::CheckConvergence(const dFunc, const dFuncPrev, const mPi, const mPiPrev,
	const dCrit, const iMaxIt, const iInfo, const iItno, mIterInfo, const dStep, const bForceTerminate)
{
	decl dconvcrit = (dFunc - dFuncPrev) / (1 + fabs(dFuncPrev));
	decl dpicrit = norm( (mPi - mPiPrev) ./ (1 + fabs(mPiPrev)) );
	decl bconverged = iItno > 0 && fabs(dconvcrit) <= dCrit && dpicrit <= sqrt(dCrit);
	decl bfinished = iItno >= iMaxIt || bconverged || bForceTerminate || ismissing(dFunc);
	decl retval = MAX_ITERATING, iprint = iInfo;

	if (mIterInfo == <>)
		mIterInfo |= iItno - 1 ~ min(dFunc, dFuncPrev) ~ MAX_INITIALIZING ~ dCrit;

	if (bfinished)
	{
		if (ismissing(dFunc))
			retval = MAX_FUNC_FAIL;
		else
		{
			retval = bconverged ? MAX_CONV : iItno >= iMaxIt ? MAX_MAXIT
				: fabs(dconvcrit) <= 100 * dCrit ? MAX_WEAK_CONV : MAX_NOCONV;
			if (dStep == 0)				  // no improvement in linesearch, in switching that could also mean failure
				retval = MAX_WEAK_CONV;
		}
		mIterInfo |= iItno ~ dFunc ~ retval ~ dCrit;
	}
	else
	{	// use MAX_ITERATING as a `still iterating' code, will not be used upon termination
		if (iItno == 1)
			mIterInfo |= iItno ~ dFunc ~ retval ~ dCrit;
		else if (sm_bTraceProgress && iItno > 1)
		{
			if (iItno < 10 || (iItno < 100 && imod(iItno, 10) == 0) || (iItno < 1000 && imod(iItno, 100) == 0)
				|| (iItno < 10000 && imod(iItno, 1000) == 0) || (iItno < 100000 && imod(iItno, 10000) == 0) )
			mIterInfo |= iItno ~ dFunc ~ retval ~ dCrit;
		}
	}

	if (iprint > 0)
	{
		if (iItno >= 10000)
			iprint = max(10000, iprint);
		else if (iItno >= 1000)
			iprint = max(1000, iprint);
		else if (iItno >= 100)
			iprint = max(100, iprint);
		else if (iItno >= 10)
			iprint = max(10, iprint);
		if (iItno > 0 && (iItno == 1 || imod(iItno, iprint) == 0 || bfinished || iInfo == 1) )
		{
			println("it", "%-5d", iItno, " f=", "%#16.10g", dFunc, "  change=", "%#12.4g", dconvcrit,
				" pi-change", "%#12.4g", dpicrit, dStep != 1 ? sprint("  step=", dStep) : "");
			if (bfinished)
				println(retval == MAX_WEAK_CONV ?
					"Weak convergence (downward step or stalled in line search)" : ::MaxConvergenceMsg(retval));
		}
	}

	return {dconvcrit, dpicrit, bfinished, retval, mIterInfo};
}

/**    Implements the grid search part of the 1Step line search.
*/
Switching::LineSearchGrid(FnLogLikAt, const vP, const vP0, const dF, const vGrid, const bTryAll)
{
	decl vp = vP, vp_0 = vP0, f0 = dF, vp_step, vp_d, fstep, k, step_k, stage, dstep = 1;
	decl cloglik = 0;

//if (bTryAll)
//	println("\ntrying all");
	vp_d = vp - vp_0;
	foreach (step_k in vGrid[k])
	{
		vp_step = vp_0 + step_k * vp_d;
		fstep = FnLogLikAt(vp_step);
//if (bTryAll)
//	println("step=", "%10.3g", step_k, "%#25.15g", f0, "%#25.15g", fstep - dF); 
		++cloglik;
		if (fstep > f0)			// yes: better function value
		{	
			vp = vp_step;
			dstep = step_k;
			f0 = fstep;
		}
		else if (step_k > 1)	// assume no (further) improvement possible
		{
			if (!bTryAll) break;
		}
	}
	return {vp, f0, dstep, 0, cloglik};
}
/**    Implements the 1Step line search.
*/
Switching::LineSearch(FnLogLikAt, const vP, const vPprev, const vP0, const dF, const dFprev)
{
	if (dF - dFprev > 1e-12)
		return LineSearchGrid(FnLogLikAt, vP, vP0, dF, sm_vStep);

	decl retval = LineSearchGrid(FnLogLikAt, vP, vP0, dFprev, sm_vStepDown, TRUE);
	if (retval[1] <= dFprev)
		return {vPprev, dFprev, 0, 0, retval[4]};
	return retval;
}

/**    Implements the QStep line search.
*/
Switching::LineSearchQstep(FnLogLikAt, const vP, const vPprev, const dF, const dFprev, const dA, const dB)
{
	decl dtolf = DBL_EPSILON * 1e4 * (fabs(dF) + fabs(dFprev)) / 2, dtolx = 0.3, cfunc = 1;
	decl vp_d = vP - vPprev, df10 = dF - dFprev;
	
	decl f0 = dFprev, f1 = dF, f2 = FnLogLikAt(vPprev + 2 * vp_d);
	decl xq, fq, q = -f0 + 2 * f1 - f2;
	if (df10 > dtolf && f2 - f1 > df10 + dtolf)	// upward acceleration: try upperbound
		xq = dB;
	else if (q <= dtolf)		// almost flat or going to minimum
		xq = df10 > -dtolf ? dB / 2 : dA / 2;
	else						// try quadratic approximation
		xq = double(setbounds(0.5 * (-3 * f0 + 4 * f1 - f2) / q, dA, 0.5 * (dB + 2)));

	decl x = <0;1;2> | xq, f = f0 | f1 | f2 | -.Inf, imax = maxcindex(f), xmax = x[imax];

	// if xq is a good prediction: skip function evaluation, else use it
	if (fabs(xq - xmax) > dtolx)
	{
		f[3] = fq = FnLogLikAt(vPprev + xq * vp_d);		
		++cfunc;
		if (fq > f[imax])
			imax = 3, xmax  = xq;
	}

	return {vPprev + xmax * vp_d, f[imax], xmax, 0, cfunc};
}

/**    Implements variant 3 of the SQUAREM line search.
*/
Switching::Squarem3(FnLogLikAt, FnUpdatePar, const vP1, const vP0, const dF0, const bGlobal, const vLo, const vHi)
{
	decl vp2 = FnUpdatePar(vP1);
	decl r = vP1 - vP0, v = (vp2 - vP1) - r;
	decl alpha = -norm(r) / norm(v);
	decl vp1 = vP0 - 2 * alpha * r + alpha^2 * v;
	// problem: the loglik becomes undefined if a parameter moves out of bounds
	// solution: shrink alpha while the parameters are outside their bounds,
	// or use global adjustment below
	if (sizerc(vLo))   	// also useful when using global version
	{
		while (alpha != -1 && any(vp1 .<= vLo .|| vp1 .>= vHi))
		{
			alpha = (alpha - 1) / 2;
			vp1 = vP0 - 2 * alpha * r + alpha^2 * v;
		}
	}
	vp2 = FnUpdatePar(vp1);
	decl df = FnLogLikAt(vp2), cupdate = 2;
	if (bGlobal)
	{
		if (alpha > -1)
			alpha = -1;
		else if (ismissing(df) || df < dF0)
		{
			while (df <= dF0 && alpha < -1)
			{
				alpha = (alpha - 1) / 2;
				vp2 = FnUpdatePar(vP0 - 2 * alpha * r + alpha^2 * v);
				++cupdate;
				df = FnLogLikAt(vp2);
			}
		}
	}
	return {vp2, df, alpha, cupdate, cupdate - 1};
}

/**    Implements the parabolic line search.
*/
Switching::Parabolic(FnLogLikAt, const vP2, const vP1, const vP0)
{
	decl i, t = 1, t_best = 1, vp_new, df_new, vp_best, df_best;
	
	vp_new = 0.01 * vP0 - 0.22 * vP1 + 1.21 * vP2;
	df_new = FnLogLikAt(vp_new);
	vp_best = vP2;
	df_best = FnLogLikAt(vP2);

	for (i = 1; df_new > df_best; ++i)
	{
		vp_best = vp_new;  df_best = df_new;  t_best = t;
		t = 1 + 1.5^i * 0.1;
		vp_new = (1 - t)^2 * vP0 + 2 * t * (1 - t) * vP1 + t^2 * vP2;
		df_new = FnLogLikAt(vp_new);
	}
	return {vp_best, df_best, t_best, 0, 2 + i};
}

/**    Implements the Powell line search.
*/
Switching::Powell(FnLogLikAt, const vP, const vPprev, const dF, const dFprev)
{
	decl vp_d = vP - vPprev, dstep, fstep, cfunc;

	decl func = [=](const dStep)
	{
		return FnLogLikAt(vPprev + dStep * vp_d);	
	};
	
	[dstep, fstep, cfunc] = MaxScalarPowell(func, 0.0, 1.0, dFprev, dF, -1, 8, m_dLineSearchTol);

	return {vPprev + dstep * vp_d, fstep, dstep, 0, cfunc};
}

/**    Implements the Brent line search.
*/
Switching::Brent(FnLogLikAt, const vP, const vPprev, const dF, const dFprev)
{
	decl vp_d = vP - vPprev, dstep, fstep, cfunc;

	decl func = [=](const dStep)
	{
		return FnLogLikAt(vPprev + dStep * vp_d);	
	};
	
	[dstep, fstep, cfunc] = MaxScalarBrent(func, 0.0, 1.0, dFprev, dF, -1, 8, m_dLineSearchTol);

	return {vPprev + dstep * vp_d, fstep, dstep, 0, cfunc};
}

/**    Implements the switching algorithm, see Table 4.
@param FnLogLikAt: function that evaluates the objective function, FnLogLikAt(vP) returns double
@param FnUpdatePar: function that evaluates the update, FnUpdatePar(vP)	returns updated vP
@param FnCheckPar: function that returns the coefficients to check convergence, FnCheckPar(vP) returns vP or a derived value
@param vP: column vector, the initial values
@param vLo: <> or column vector with lower bounds of coefficients (only used with LS_2STEPSQ3,LS_SQS3,LS_SQS3G)
@param vHi: <> or column vector with upper bounds of coefficients (only used with LS_2STEPSQ3,LS_SQS3,LS_SQS3G)
@param dCrit: double, convergence tolerance \epsilon_1
@param iMaxIt: int, maximum number of iterations
@param iInfo: int, <=0: don't print anything, >1: print iteration on a log10 scale
@param FnUserLineSearch: user-defined line search function (only used with LS_USERi)
@returns array with 3 values: convergenced coefficients, convergenced function value, return code
(MAX_ITERATING, MAX_CONV, MAX_WEAK_CONV, MAX_MAXIT, MAX_NOCONV, MAX_FUNC_FAIL).
*/
Switching::Estimate(FnLogLikAt, FnUpdatePar, FnCheckPar, vP, const vLo, const vHi, const dCrit, const iMaxIt, const iInfo,
	FnUserLineSearch)
{
	decl retval, bfinished = FALSE, itno, dfunc, dfunc_candidate, dfunc_prev, dfunc_0, vp, vp_candidate, vp_0, vp_prev,
		vp_0_0, dconvcrit, dstep = 1, iter_info = <>, dpicrit, bforceterminate = FALSE;
	decl cupdate_it, cloglik_it, cupdate_ls, cloglik_ls, cupdate = 0, cloglik = 0;
	decl bdoeval = !(m_iLineSearchMode == LS_PARABOLIC || m_iLineSearchMode == LS_PARABOLIC2);
	decl cextra = m_iLineSearchMode == LS_1STEP1 ? 1 : m_cExtraUpdate;

	// ideally, if the initial values don't satisfy the model, we set dfunc to -.Inf
	vp = vp_candidate = vp_0 = vec(vP);
	dfunc = dfunc_candidate = dfunc_0 = FnLogLikAt(vp);
	cloglik = 1;
	
	if (iInfo > 0)
	{
		println("\nSWITCHING ALGORITHM v", GetVersion(), ", eps=", dCrit, " line search=", GetLineSearchName(m_iLineSearchMode));
		if (m_cExtraUpdate)
			println("With extra update");
		println("start   f=", "%#16.10g", dfunc);
	}
	
	for (itno = 0; !bfinished; )
	{
		vp_0_0 = vp_0;
		dfunc_0 = dfunc_candidate;  vp_0 = vp_candidate;
		dfunc_prev = dfunc;  vp_prev = vp;

		// EM step
		cloglik_it = cupdate_ls = cloglik_ls = 0;
		vp = FnUpdatePar(vp);  cupdate_it = 1;
		if (cextra)
			vp = FnUpdatePar(vp), ++cupdate_it;
		vp_candidate = vp;
		if (bdoeval)
			dfunc = dfunc_candidate = FnLogLikAt(vp), ++cloglik_it;

		// optional linesearch: must enter when dstepcrit = (dfunc - dfunc_prev) / (1 + fabs(dfunc_prev))<=0
		// because we would go down otherwise
		dstep = 1;
		if (itno >= m_cLineSearchWarmUp && (m_dLineSearchCrit <= 0 ||
			(dfunc - dfunc_prev) <= m_dLineSearchCrit * (1 + fabs(dfunc_prev))))
		{
			switch_single (m_iLineSearchMode)
			{
				case LS_1STEP:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearch(FnLogLikAt, vp, vp_prev, vp_0, dfunc, dfunc_prev);
				case LS_1STEP1:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearch(FnLogLikAt, vp, vp_prev, vp_0, dfunc, dfunc_prev);
				case LS_2STEP:
				{
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearch(FnLogLikAt, vp, vp_prev, vp_0, dfunc, dfunc_prev);
					decl cupdate_2, cloglik_2;
					[vp, dfunc, dstep, cupdate_2, cloglik_2] = LineSearch(FnLogLikAt, vp, vp_prev, vp_prev, dfunc, dfunc_prev);
					cupdate_ls += cupdate_2, cloglik_ls += cloglik_2;
				}
				case LS_2STEPSQ3:
				{
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearch(FnLogLikAt, vp, vp_prev, vp_0, dfunc, dfunc_prev);
					decl cupdate_2, cloglik_2;
					[vp, dfunc, dstep, cupdate_2, cloglik_2] = Squarem3(FnLogLikAt, FnUpdatePar, vp, vp_prev, dfunc_prev, FALSE, vLo, vHi);
					cupdate_ls += cupdate_2, cloglik_ls += cloglik_2;
				}
				case LS_1STD:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearch(FnLogLikAt, vp, vp_prev, vp_prev, dfunc, dfunc_prev);
				case LS_SQS3:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Squarem3(FnLogLikAt, FnUpdatePar, vp, vp_prev, dfunc_prev, FALSE, vLo, vHi);
				case LS_SQS3G:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Squarem3(FnLogLikAt, FnUpdatePar, vp, vp_prev, dfunc_prev, TRUE, vLo, vHi);
				case LS_PARABOLIC:
					// extra update
					vp_candidate = vp = FnUpdatePar(vp);	 ++cupdate_it;
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Parabolic(FnLogLikAt, vp, vp_0, vp_0_0);
				case LS_PARABOLIC2:
					// extra update
					vp = FnUpdatePar(vp);	 ++cupdate_it;
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Parabolic(FnLogLikAt, vp, vp_prev, vp_0_0);
				case LS_POWELL:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Powell(FnLogLikAt, vp, vp_0, dfunc, dfunc_0);
				case LS_BRENT:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Brent(FnLogLikAt, vp, vp_0, dfunc, dfunc_0);
				case LS_BRENTSTD:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = Brent(FnLogLikAt, vp, vp_prev, dfunc, dfunc_prev);
				case LS_QSTEP:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = LineSearchQstep(FnLogLikAt, vp, vp_0, dfunc, dfunc_0, -1, 8);
				case LS_USER1:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = FnUserLineSearch(this, vp, vp_prev, vp_0, dfunc, dfunc_prev, dfunc_0);
				case LS_USER2:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = FnUserLineSearch(this, vp, vp_prev, vp_0, dfunc, dfunc_prev, dfunc_0);
				case LS_USER3:
					[vp, dfunc, dstep, cupdate_ls, cloglik_ls] = FnUserLineSearch(this, vp, vp_prev, vp_0, dfunc, dfunc_prev, dfunc_0);
			}
		}
		++itno;
		cupdate += cupdate_it + cupdate_ls;
		cloglik += cloglik_it + cloglik_ls;

		if (dfunc < dfunc_prev && itno > 0)
		{	// assume initial values satisfy restrictions 
			bforceterminate = TRUE;
			vp = vp_prev;  dfunc = dfunc_prev;
			dstep = 0;
		}

		[dconvcrit, dpicrit, bfinished, retval, iter_info] = CheckConvergence(dfunc, dfunc_prev,
			FnCheckPar(vp), FnCheckPar(vp_prev), dCrit, iMaxIt, iInfo, itno, iter_info, dstep, bforceterminate);
	}
	
	// track performance
	m_mMaxTrace = iter_info;
    m_iResult = retval;
	m_cIterCount = itno;
	m_cUpdateCount = cupdate;
	m_cLogLikCount = cloglik;

	return {vp, dfunc, retval};
}
///////////////////////////////////////////////////////////////////////
