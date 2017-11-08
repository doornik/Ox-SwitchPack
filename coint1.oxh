#ifndef COINT1_H_INCLUDED
#define COINT1_H_INCLUDED

class Coint
{
	Coint(const mY, const cLag=2, const iDet=2);
	static RRR(const mY, const mX, const cR);
	MapToArgs(const vP);
	MapToPar(const mA, const mB);
	SetRank(const iR);
	Restrict(const amG, const amH);
	GenerateVAR1(const mEps, const mAlpha, const mBeta);
	UpdateAlpha(const mBeta, const mOmega);
	virtual UpdatePar(vP);
	UpdatePar_BetaSwitching(vP);
	OmegaAt(const mA, const mB);
	LogDetAt(const mA, const mB);
	virtual LogLikAt(const vP);
	LineSearchBetaType(const objSwitching, const bGrid, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0);
	LineSearch1Beta(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0);
	LineSearchQBeta(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0);
	LineSearchOptLsq(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0, const dA=-.Inf, const dB=.Inf);
	LineSearchOpt(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0);
	CheckPar(const vP);

	decl m_mY, m_mZ0, m_mZ1, m_cZ0, m_cZ1, m_iR;
	decl m_mAlpha_I1, m_mBeta_I1;
	decl m_cA, m_cB, m_amG, m_amH, m_mG, m_mH, m_mGbar;
};

#endif // COINT1_H_INCLUDED