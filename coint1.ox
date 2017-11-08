#include <oxstd.oxh>
#include <arma.oxh>

#include "coint1.oxh"


///////////////////////////////////////////////////////////////////////
//
/**	Solves the I(1) pencil |\lambda S_{11} - S_{10}S_{00}^{-1}S_{01}|.
*/
Coint::RRR(const mY, const mX, const cR)
{
	decl ct = sizer(mY), s01 = (mY'mX) / ct, s11 = (mX'mX) / ct, s00_inv, eval, evec,
		alpha, beta, dlogdet, scale;
	
	s00_inv = invertsym((mY'mY) / ct, &dlogdet);

	if (eigensymgen(s01' s00_inv * s01, s11, &eval, &evec) < 0)
	{	
		return {zeros(sizec(mY), cR), zeros(sizec(mX), cR), dlogdet};
	}
	beta = evec[][ : cR - 1];
	alpha = s01 * beta;
	scale = maxc(fabs(beta));
	beta ./= scale;  alpha .*= scale;
	dlogdet += sumr(log(1 - eval[ : cR - 1]));
	
	return {alpha, beta, -double(dlogdet)};
}

Coint::Coint(const mY, const cLag, const iDet)
{
	m_mY = mY;
	m_mZ0 = diff(mY, 1);	   // DY 
	m_mZ1 = lag(mY, 1);		   // Y_1
	
	decl z2 = <>;
	if (iDet == 1)			   // restricted constant
		m_mZ1 ~= 1;
	else if (iDet == 2)		   // unrestricted constant, restricted trend
	{
		m_mZ1 ~= range(1, sizer(mY))' - cLag;
		z2 ~= ones(sizer(mY), 1);
	}	
	if (cLag > 1)
	{	// add lagged differences
		z2 ~= lag(m_mZ0, range(1, cLag - 1));
	}
	// drop lost observations from Z0,Z1,Z2
	m_mZ0 = m_mZ0[cLag : ][];
	m_mZ1 = m_mZ1[cLag : ][];
	z2 = z2[cLag : ][];
	
	if (sizer(z2) > 1)
	{	// partial out lagged differences and unrestricted terms
		decl b;
		olsc(m_mZ0, z2, &b);  m_mZ0 -= z2 * b;
		olsc(m_mZ1, z2, &b);  m_mZ1 -= z2 * b;
	}
	m_cZ0 = sizec(m_mZ0);
	m_cZ1 = sizec(m_mZ1);

	// no restrictions yet
	m_mG = m_mH = 1;
	m_mAlpha_I1 = m_mBeta_I1 = <>;
	m_iR = 0;
}
Coint::MapToArgs(const vP)
{
	decl ma = shape(m_mG * vP[ : m_cA - 1], m_cZ0, m_iR), mb = shape(m_mH * vP[m_cA :], m_cZ1, m_iR);
	return { ma, mb } ;
}
Coint::MapToPar(const mA, const mB)
{
	decl theta, phi;
	if (m_mG == 1)
		theta = vec(mA);
	else
		olsc(vec(mA), m_mG, &theta); 

	if (m_mH == 1)
		phi = vec(mB);
	else
		olsc(vec(mB), m_mH, &phi); 

	return theta | phi;
}
Coint::SetRank(const iR)
{
	m_iR = iR;	   //  rank

	decl alpha, beta, dlogdet;
	[alpha, beta, dlogdet] = RRR(m_mZ0, m_mZ1, iR);

	m_mAlpha_I1 = alpha;
	m_mBeta_I1 = beta;
	m_cA = sizerc(alpha);
	m_cB = sizerc(beta);

	return MapToPar(alpha, beta);
}
Coint::Restrict(const amG, const amH)
{
	decl theta, phi, i;
	
	m_mGbar = <>;
	
	if (sizeof(amG) == m_iR)
	{
		m_amG = amG;
		for (i = 1, m_mG = m_amG[0]; i < m_iR; ++i)
			m_mG = diagcat(m_mG, m_amG[i]); 
		m_cA = sizec(m_mG);
	}
	else if (sizeof(amG) == 1)	// common restriction
	{
		m_amG = amG;
		decl mg = m_amG[0];
		for (i = 1, m_mG = mg; i < m_iR; ++i)
			m_mG = diagcat(m_mG, mg); 
		m_cA = sizec(m_mG);
		m_mGbar = mg * invertsym(mg'mg);
	}
	else
	{
		m_amG = {};
		m_mG = 1;
		m_cA = m_cZ0 * m_iR;
		if (sizeof(amG) != 0)
			oxrunerror("alpha restrictions");
	}
	if (sizeof(amH) == m_iR)
	{
		m_amH = amH;
		for (i = 1, m_mH = m_amH[0]; i < m_iR; ++i)
			m_mH = diagcat(m_mH, m_amH[i]); 
		m_cB = sizec(m_mH);
	}
	else
	{
		m_amH = {};
		m_mH = 1;
		m_cB = m_cZ1 * m_iR;
		if (sizeof(amH) != 0)
			oxrunerror("beta restrictions");
	}
	return MapToPar(m_mAlpha_I1, m_mBeta_I1);
}
Coint::GenerateVAR1(const mEps, const mAlpha, const mBeta)
{
	decl cy = sizec(m_mY), ct = sizer(mEps), mpi = unit(cy) + mAlpha * mBeta[ : cy - 1][]';
	decl mean = 0;
	// then add the restricted part to mean
	if (sizec(mBeta) && sizer(mBeta) > cy)
		mean += (range(1, ct)' - 2) * (mBeta[cy : ][] * mAlpha');

	// pre-sample y_t, rest is zero; indices of model
	decl y0 = m_mY[ : 1][] | zeros(ct - 1, cy);
	decl idx_y = range(0, cy - 1), idx_var = idx_y, idx_lag = ones(1, cy);

//	println(m_mY ~ modelforc(mean + mEps, y0, idx_y, idx_var, idx_lag, mpi, 1));
	return modelforc(mean + mEps, y0, idx_y, idx_var, idx_lag, mpi, 1);
}

Coint::CheckPar(const vP)
{
	decl ma, mb;
	[ma, mb] = MapToArgs(vP);
	return vec(ma * mb');
}
Coint::UpdateAlpha(const mBeta, const mOmega)
{
	decl alpha, p_inv, theta;

	// alpha|beta step
	if (sizeof(m_amG) == 0)	 // unrestricted alpha
	{
		olsc(m_mZ0, m_mZ1 * mBeta, &alpha);
		alpha = alpha';
		theta = vec(alpha);
	}
	else if (sizeof(m_amG) == 1)	// common alpha restriction
	{	
		decl omega_inv = invertsym(mOmega), mg = m_amG[0];
		decl mg_dbar = omega_inv * mg * invertsym(outer(mg', omega_inv));
		olsc(m_mZ0 * mg_dbar, m_mZ1 * mBeta, &alpha);
		alpha = alpha';
		theta = vec(alpha);
		alpha = mg * alpha;
	}
	else  // restricted alpha
	{
		p_inv = solvelu(choleski(mOmega), 0, 0, unit(m_cZ0));
		olsc(vec(p_inv * m_mZ0'), ((m_mZ1 * mBeta) ** p_inv) * m_mG, &theta);
		alpha = shape(m_mG * theta, m_cZ0, m_iR);
	}

	return {alpha, theta};
}
Coint::UpdatePar(vP)
{
	decl alpha, beta, e, omega, p_inv, phi, theta;
	[alpha, beta] = MapToArgs(vP);

	e = m_mZ0 - m_mZ1 * beta * alpha';
	omega = (e'e) / sizer(e);

	// beta|alpha step
	p_inv = solvelu(choleski(omega), 0, 0, unit(m_cZ0));
	olsc(vec(m_mZ0 * p_inv'), ((p_inv * alpha) ** m_mZ1) * m_mH, &phi);
	beta = shape(m_mH * phi, m_cZ1, m_iR);

	// alpha|beta step
	e = m_mZ0 - m_mZ1 * beta * alpha';
	[alpha, theta] = UpdateAlpha(beta, (e'e) / sizer(e));

	return theta | phi;
}
Coint::UpdatePar_BetaSwitching(vP)
{
	decl alpha, beta, e, omega, z0, z1, z1b, beta_not_i, alpha_i, phi_i, phi, i, b, theta, mg_dbar = 1;
	[alpha, beta] = MapToArgs(vP);

	e = m_mZ0 - m_mZ1 * beta * alpha';
	omega = (e'e) / sizer(e);
	if (sizeof(m_amG) == 1)
	{
		decl omega_inv = invertsym(omega), mg = m_amG[0];
		mg_dbar = omega_inv * mg * invertsym(outer(mg', omega_inv));
	}
	// beta step
	phi = <>;
	for (i = 0; i < m_iR; ++i)
	{
		beta_not_i = dropc(beta, i);
		z1b = m_mZ1 * beta_not_i;
		// partial other betas out
		olsc(m_mZ0 * mg_dbar, z1b, &b);  z0 = m_mZ0 * mg_dbar - z1b * b;
		olsc(m_mZ1, z1b, &b);  z1 = m_mZ1 - z1b * b;

		[alpha_i, phi_i] = RRR(z0, z1 * m_amH[i], 1);
		beta[][i] = m_amH[i] * phi_i;
		phi |= phi_i;
	}

	// alpha|beta step
	e = m_mZ0 - m_mZ1 * beta * alpha';
	[alpha, theta] = UpdateAlpha(beta, omega);

	return theta | phi;
}
Coint::OmegaAt(const mA, const mB)
{
	decl e = m_mZ0 - m_mZ1 * (mB * mA');
	return (e'e) / sizer(e);
}
Coint::LogDetAt(const mA, const mB)
{
	decl e, dlogdet, omega;
	e = m_mZ0 - m_mZ1 * (mB * mA');
	omega = (e'e) / sizer(e);
	invertsym(omega, &dlogdet);
	return -dlogdet; 
}
Coint::LogLikAt(const vP)
{
	decl ma, mb, mpi, e, dlogdet;
	[ma, mb] = MapToArgs(vP);
	return LogDetAt(ma, mb);
}
Coint::LineSearchBetaType(const objSwitching, const bGrid, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0)
{
	decl ma, mb, mpi, e, omega, dlogdet;
	[ma, mb] = MapToArgs(vP);
	e = m_mZ0 - m_mZ1 * (mb * ma');
	omega = (e'e) / sizer(e);

	// logdet for beta only, keep omega fixed (restricted alpha)
	decl fn_loglik = [=](const vP)
	{
		decl beta0 = MapToArgs(zeros(m_cA, 1) | vP)[1];
		return LogDetAt(UpdateAlpha(beta0, omega)[0], beta0); 
	};
	decl retval;
	if (bGrid)
		retval = objSwitching.LineSearch(fn_loglik, vP[m_cA : ], vPprev[m_cA : ], vP0[m_cA : ], dF, dFprev);
	else
		retval = objSwitching.LineSearchQstep(fn_loglik, vP[m_cA : ], vP0[m_cA : ], dF, dF0, -1, 8);
	// add alpha back in
	if (retval[2] == 1)
		retval[0] =	vP;
	else if (retval[2] == 0)
		retval[0] =	vPprev;
	else
		retval[0] =	UpdateAlpha(MapToArgs(zeros(m_cA, 1) | retval[0])[1], omega)[1] | retval[0];
	return retval;
}
Coint::LineSearch1Beta(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0)
{
	return LineSearchBetaType(objSwitching, TRUE, vP, vPprev, vP0, dF, dFprev, dF0);
}
Coint::LineSearchQBeta(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0)
{
	return LineSearchBetaType(objSwitching, FALSE, vP, vPprev, vP0, dF, dFprev, dF0);
}
Coint::LineSearchOptLsq(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0, const dA, const dB)
{
	decl fn_solve_linesearch = [=](const vP, const vPprev, const mOmegaInv)
	{
		decl alpha_0, beta_0, alpha_l, beta_l, dalpha, dbeta, roots; 
		[alpha_0, beta_0] = MapToArgs(vPprev);
		[dalpha, dbeta] = MapToArgs(vP - vPprev);
		decl e_0 = m_mZ0 - m_mZ1 * (beta_0 * alpha_0');
		decl z1da = (m_mZ1 * beta_0) * dalpha', z1db = (m_mZ1 * dbeta) * alpha_0';
		decl z1pi1 = z1da + z1db, z1pi2 = m_mZ1 * (dbeta * dalpha');
		decl d = trace(e_0'z1pi1 * mOmegaInv), c = trace((2 * e_0'z1pi2 - z1pi1'z1pi1) * mOmegaInv),
			b = -3 * trace(z1pi2'z1pi1 * mOmegaInv), a = -2 * trace(z1pi2'z1pi2 * mOmegaInv);
	
		polyroots(d ~ c ~ b ~ a, &roots);
		// only keep real roots
		roots = deleteifc(roots, fabs(roots[1][]) .> 1e-10);
		roots = 1 ./ roots[0][];
		// select the solution that is closest to one
		roots = sortbyc(roots' ~ fabs(roots' - 1), 1);
		return roots[][0];
	};

	decl ma, mb, e, omega_inv, dlogdet, dstep, vpnew, dfunc_step;
	[ma, mb] = MapToArgs(vP);
	e = m_mZ0 - m_mZ1 * (mb * ma');
	omega_inv = invertsym((e'e) / sizer(e));

	dstep = fn_solve_linesearch(vP, vP0, omega_inv / sizer(e))[0][0];
	dstep = double(setbounds(dstep, dA, dB));
	vpnew = vP0 + dstep * (vP - vP0);
	// only accept if better
	dfunc_step = LogLikAt(vpnew);
	if (dfunc_step < dF)
	{
		dstep = 1;
		vpnew = vP;
		dfunc_step = dF;
	}
	if (dfunc_step < dFprev)
	{
		dstep = 0;
		vpnew = vPprev;
		dfunc_step = dFprev;
	}
	return {vpnew, dfunc_step, dstep, 0, 1};
}
Coint::LineSearchOpt(const objSwitching, const vP, const vPprev, const vP0, const dF, const dFprev, const dF0)
{
	return Coint::LineSearchOptLsq(objSwitching, vP, vPprev, vP0, dF, dFprev, dF0);
}
///////////////////////////////////////////////////////////////////////
