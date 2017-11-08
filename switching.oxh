#ifndef SWITCHING_OXH
#define SWITCHING_OXH

#include <maximize.oxh>

class Switching
{
	Switching();
	GetTrace();
	GetIterInfo();
	static GetLineSearchName(const iLineSearchMode);
	GetVersion();
	GetLineSearchMode();
	SetLineSearchMode(const iLineSearchMode);
	SetLineSearchWarmUp(const cLineSearchWarmUp);
	GetLineSearchCrit();
	SetLineSearchCrit(const iLineSearchCrit);
	GetLineSearchTol();
	SetLineSearchTol(const dLineSearchTol);
	SetExtraUpdate(const cExtraUpdate);
	static SetTraceProgress(const bTrace);
	
	static CheckConvergence(const dFunc, const dFuncPrev, const mPi, const mPiPrev, const dCrit,
		const iMaxIt, const iInfo, const iItno, mIterInfo, const dStep=1, const bForceTerminate=FALSE);
	static LineSearchGrid(FnLogLikAt, const vP, const vP0, const dF, const vGrid, const bTryAll=FALSE);
	LineSearch(FnLogLikAt, const vP, const vPprev, const vP0, const dF, const dFprev);
	LineSearchQstep(FnLogLikAt, const vP, const vPprev, const dF, const dFprev, const dA, const dB);
	Squarem3(FnLogLikAt, FnUpdatePar, const vP1, const vP0, const dF0, const bGlobal=TRUE, const vLo=<>, const vHi=<>);
	Parabolic(FnLogLikAt, const vP2, const vP1, const vP0);
	Powell(FnLogLikAt, const vP, const vPprev, const dF, const dFprev);
	Brent(FnLogLikAt, const vP, const vPprev, const dF, const dFprev);

	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	Estimate(FnLogLikAt, FnUpdatePar, FnCheckPar, vP, const vLo=<>, const vHi=<>, const dCrit=1e-6,
		const iMaxIt=10000, const iInfo=1, FnUserLineSearch=0);

	static decl sm_vStep = <1.2,2,4,8>;
	static decl sm_vStepDown = <-1e-1,-1e-2,-1e-4,0.01,0.1,0.4,0.99,1.01,1.1>;
	static decl sm_bTraceProgress = FALSE;
	
	decl m_iLineSearchMode, m_cLineSearchWarmUp, m_dLineSearchCrit, m_dLineSearchTol;
	decl m_mMaxTrace;
    decl m_iResult;
	decl m_cIterCount, m_cUpdateCount, m_cLogLikCount;
	decl m_cExtraUpdate;
	
	enum
	{
	    MAX_ITERATING = MAX_NOCONV + 10, MAX_INITIALIZING = -1
	};
public:	
	enum
	{
	    LS_NONE = 0, LS_1STEP, LS_1STEP1, LS_2STEP, LS_2STEPSQ3, LS_1STD, LS_SQS3, LS_SQS3G,
		LS_PARABOLIC, LS_PARABOLIC2, LS_POWELL, LS_BRENT, LS_BRENTSTD, LS_QSTEP, LS_USER1, LS_USER2, LS_USER3 
	};
}

#endif //SWITCHING_OXH