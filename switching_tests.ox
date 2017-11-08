#include <oxstd.oxh>
#include <oxprob.oxh>
#import <switching>


///////////////////////////////////////////////////////////////////////
//
class PoissonMix
{
	PoissonMix(const mData);
	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	CheckPar(const vP);

	decl m_mData;
};

PoissonMix::PoissonMix(const mData)
{
	m_mData = mData;
}
PoissonMix::CheckPar(const vP)
{
	return vP;
}
PoissonMix::UpdatePar(vP)
{
	decl p = vP[0], mu1 = vP[1], mu2 = vP[2], pi1, pi2, eps = 1e-12;
	decl vi = range(0, sizerc(m_mData) - 1)';

	decl pfac = p .* exp(vi .* (log(mu1) - log(mu2)) - (mu1 - mu2));
	pi1 = pfac ./ (pfac + (1 - p));
	pi2 = p ./ (pfac + (1 - p));
	
	// update parameters
	decl npi1 = m_mData'pi1, npi2 = m_mData'pi2, vin = vi .* m_mData;
	p = npi1 / sumc(m_mData);
	p = setbounds(p, eps, 1 - eps);

	mu1 = (vin'pi1) / npi1;
	mu2 = (vin'pi2) / npi2;
	vP[0] = p; vP[1] = mu1; vP[2] = mu2 > eps ? mu2 : eps;
	return vP;
}
PoissonMix::LogLikAt(const vP)
{
	decl p = vP[0], mu1 = vP[1], mu2 = vP[2];
	decl vi = range(0, sizerc(m_mData) - 1)';

	decl fac1 = exp(-mu1 + vi .* log(mu1) - loggamma(vi + 1));
	decl fac2 = exp(-mu2 + vi .* log(mu2) - loggamma(vi + 1));
	decl vloglik = m_mData .* log(p .* fac1 + (1 - p) .* fac2); 
	return double(sumc(vloglik));
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
class MultivariateT
{
	MultivariateT(const mData, const dNu, const bUsePX);
	MapToArgs(const vP);
	MapToPar(const vMu, const mSigma);
	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	CheckPar(const vP);

	decl m_mData, m_dNu;
	decl m_bUsePX;
};

MultivariateT::MultivariateT(const mData, const dNu, const bUsePX)
{
	m_mData = mData;
	m_dNu = dNu;
	m_bUsePX = bUsePX;
}
MultivariateT::CheckPar(const vP)
{
	return vP;
}
MultivariateT::MapToArgs(const vP)
{
	decl cp = sizec(m_mData);
	decl mu = vP[ : cp - 1], chol = setupper(unvech(vP[cp : ]), 0);
	decl chol_inv = solvelu(chol, 0, 0, unit(cp));
	decl sigma_inv = chol_inv'chol_inv;
	return { mu', sigma_inv, chol} ;
}
MultivariateT::MapToPar(const vMu, const mSigma)
{
	decl chol = choleski(mSigma);
	if (chol == 0)
		return vec(vMu) | vech(nans(sizer(mSigma), sizec(mSigma)));
	return vec(vMu) | vech(chol);
}
MultivariateT::UpdatePar(vP)
{
	decl cp = sizec(m_mData), mu, sigma, sigma_inv;
	[mu, sigma_inv] = MapToArgs(vP);
	// E step
	decl w = (m_dNu + cp) ./ (m_dNu + outer(m_mData - mu, sigma_inv, 'd'))';
	decl sumw = sumc(w);
	// M step
	mu = sumc(w .* m_mData) / sumw;
	decl fac = m_bUsePX ? sumw : sizer(m_mData);
	sigma = outer( sqrt(w) .* (m_mData - mu), <>, 'o') / fac;
	
	return MapToPar(mu, sigma);
}
MultivariateT::LogLikAt(const vP)
{
	decl cp = sizec(m_mData), mu, chol, sigma_inv;
	[mu, sigma_inv, chol] = MapToArgs(vP);
	decl logdet = 2 * sumr(log(diagonal(chol))) + (m_dNu + cp) * meanr(log(m_dNu + outer(m_mData - mu, sigma_inv, 'd')));
	return double(-logdet);
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
class Parafac
{
	Parafac(const mData_I_JK, const cI, const cJ, const cK, const cF, const bConcentrate=TRUE);
	static Left(const mData_I_JK, const cI, const cJ, const cK);
	static Transpose(const mData_I_JK, const cI, const cJ, const cK);
	static Olsr(const mData, const mB, const mC);
	static ColKron(const mA, const mB);
	static Normalize(const mA, const mB, const mC);
	MapToArgs(const vP);
	MapToPar(const mA, const mB, const mC);
	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	CheckPar(const vP);

	decl m_mData_I_JK, m_mData_J_IK, m_mData_K_IJ;
	decl m_cI, m_cJ, m_cK, m_cF;
	decl m_bConcentrate;
};

Parafac::Parafac(const mData_I_JK, const cI, const cJ, const cK, const cF, const bConcentrate)
{
	m_mData_I_JK = mData_I_JK;
	m_mData_J_IK = Transpose(mData_I_JK, cI, cJ, cK);
	m_mData_K_IJ = Left(mData_I_JK, cI, cJ, cK);
	m_cI = cI;
	m_cJ = cJ;
	m_cK = cK;
	m_cF = cF;
	m_bConcentrate = bConcentrate;
}
Parafac::Left(const mData_I_JK, const cI, const cJ, const cK)
{
	// map I_JK to K_JI	or K_IJ
	//	 =0	|		|		|		  |		 =0	|		|		|		  |  		 =0	|		|		|		  |
	//	i=1	|	k=0	|	k=1	|	k=2	  |;vecr:k=1|	i=0	|	i=1	|	i=2	  |; or vec:k=1	|	j=0	|	j=1	|	j=2	  |  
	//	 =2	|		|		|		  |		 =2	|		|		|		  |  		 =2	|		|		|		  |
	//		 j=0,1,	 j=0,1,	 j=0,1,				 j=0,1,	 j=0,1,	 j=0,1,     			 i=0,1,	 i=0,1,	 i=0,1,    
	// Get k as first index: extract K matrices and store in rows
	decl data = new matrix[cK][cI * cJ];
	for (decl k = 0; k < cK; ++k)
		data[k][] = vec(mData_I_JK[][k * cJ : (k + 1) * cJ - 1])';
	return data;
}
Parafac::Transpose(const mData_I_JK, const cI, const cJ, const cK)
{
	// map I_JK to J_IK
	decl data = new matrix[cJ][cI * cK];
	for (decl k = 0; k < cK; ++k)
		data[][k * cI : (k + 1) * cI - 1] = mData_I_JK[][k * cJ : (k + 1) * cJ - 1]';
	return data;
}
Parafac::Olsr(const mData, const mB, const mC)
{
	return (mData * ColKron(mC, mB)) * invertgen((mC'mC) .* (mB'mB), 3)';
}
Parafac::ColKron(const mA, const mB)
{
	decl cf = sizec(mA), mc = new matrix[sizer(mA) * sizer(mB)][cf];
	for (decl i = 0; i < cf; ++i)
		mc[][i] = mA[][i] ** mB[][i];
	return mc;
}
Parafac::Normalize(const mA, const mB, const mC)
{
	decl norm_a = sqrt(sumsqrc(mA)), norm_b = sqrt(sumsqrc(mB));
	decl skipc = mC == <>, norm_c = skipc ? ones(1, sizec(mA)) : sqrt(sumsqrc(mC));
	decl scale = pow(norm_a .* norm_b .* norm_c, skipc ? 1/2 : 1/3);
	return {mA .* (scale ./ norm_a), mB .* (scale ./ norm_b), skipc ? <> : mC .* (scale ./ norm_c) };
}
Parafac::CheckPar(const vP)
{
	return vP;
}
Parafac::MapToArgs(const vP)
{
	decl ma = shape(vP, m_cI, m_cF), mb = shape(vP[m_cI * m_cF :], m_cJ, m_cF);
	if (m_bConcentrate)
		return { ma, mb, Olsr(m_mData_K_IJ, ma, mb)} ;
	return { ma, mb, shape(vP[(m_cI + m_cJ) * m_cF :], m_cK, m_cF)} ;
}
Parafac::MapToPar(const mA, const mB, const mC)
{
	if (m_bConcentrate)
		return vec(mA) | vec(mB);
	return vec(mA) | vec(mB) | vec(mC);
}
Parafac::UpdatePar(vP)
{
	decl ma, mb, mc;
	[ma, mb, mc] = MapToArgs(vP);

//	olsr(m_mData_I_JK, ColKron(mc, mb)', &ma);
//	olsr(m_mData_J_IK, ColKron(mc, ma)', &mb);
//	olsr(m_mData_K_IJ, ColKron(mb, ma)', &mc);

	ma = Olsr(m_mData_I_JK, mb, mc);
	mb = Olsr(m_mData_J_IK, ma, mc);
	if (m_bConcentrate)
		[ma, mb] = Normalize(ma, mb, <>);
	else
	{
		mc = Olsr(m_mData_K_IJ, ma, mb);
		[ma, mb, mc] = Normalize(ma, mb, mc);
	}
	return MapToPar(ma, mb, mc);
}
Parafac::LogLikAt(const vP)
{
	decl ma, mb, mc;
	[ma, mb, mc] = MapToArgs(vP);
	// same as -trace(eps'eps/T);
	return -norm(m_mData_I_JK - ma * ColKron(mc, mb)', 'F') / sqrt(sizerc(m_mData_I_JK));
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
class LowRank
{
	LowRank(const mY, const mX, const cQ, const iR, const cS);
	MapToArgs(const vP);
	MapToPar(const mA, const mB);
	virtual StartPar(iR=0);
	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	CheckPar(const vP);

	decl m_mY, m_mX, m_cQ, m_iR, m_cS;
};

LowRank::LowRank(const mY, const mX, const cQ, const iR, const cS)
{
	m_mY = mY;	   //  n x 1
	m_mX = mX;	   //  n x qs
	m_cQ = cQ;	   //  A[q][r]
	m_iR = iR;	   //  rank
	m_cS = cS;	   //  B[s][r]
}
LowRank::MapToArgs(const vP)
{
	decl ma = shape(vP, m_cQ, m_iR), mb = shape(vP[m_cQ * m_iR :], m_cS, m_iR);
	return { ma, mb } ;
}
LowRank::MapToPar(const mA, const mB)
{
	return vec(mA) | vec(mB);
}
LowRank::StartPar(iR)
{
	decl vc, mc, mu, mw, mv, ma, mb, cxs;

	cxs = max(sizec(m_mX) - sizer(m_mX), 0);
	olsc(m_mY | zeros(cxs, 1), m_mX | zeros(cxs, sizec(m_mX)), &vc);
	mc = shape(vc, m_cQ, m_cS);
	
	// map to lower rank
	decsvd(mc, &mu, &mw, &mv);
	if (iR == 0)
		iR = m_iR;
	if (iR < sizerc(mw)) mw[iR : ] = 0;
	ma = (mu .* mw)[][ : m_iR - 1];  mb = mv[][ : m_iR - 1];

	return MapToPar(ma, mb);
}
LowRank::CheckPar(const vP)
{
	decl ma, mb;
	[ma, mb] = MapToArgs(vP);
	return vec(ma * mb');
}
LowRank::UpdatePar(vP)
{
	decl ma, mb, va, vbt, b, e;
	[ma, mb] = MapToArgs(vP);

	// estimate new B | A
	olsc(m_mY, m_mX * (unit(m_cS) ** ma), &vbt);
	mb = reshape(vbt, m_cS, m_iR);

	// estimate new A | B
	olsc(m_mY, m_mX * (mb ** unit(m_cQ)), &va);
	ma = shape(va, m_cQ, m_iR);

	return MapToPar(ma, mb);
}
LowRank::LogLikAt(const vP)
{
	decl ma, mb, mc;
	[ma, mb] = MapToArgs(vP);
	mc = ma * mb';

	return -sqrt(double(sumsqrc(m_mY - m_mX * vec(mc))) / sizer(m_mY)); 
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
class MySwitching : Switching
{
	MySwitching(iLineSearchMode);
	static Report(mResult, time);
}
MySwitching::MySwitching(iLineSearchMode)
{
	Switching();
	SetLineSearchMode(iLineSearchMode);
	SetLineSearchCrit(0);
	SetLineSearchWarmUp(0);
}
MySwitching::Report(mResult, time)
{
	mResult[][0] -= max(mResult[][0]);
	println("%cf", {"%10.0f","%10.0f","%10.0f","%10.0f","%10.0f","%10.2f"}, "%c", {"iters","updates","logliks","failures","below","CPU(s)"},
		meanc(mResult[][2:]) ~ sumc(mResult[][1] .> MAX_WEAK_CONV) ~ .NaN ~ (timer() - time) / 100);
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
RunPoissonExp(const iLineSearchMode, const cM, bExtra=FALSE)
{
	decl data = <162;267;271;185;111;61;27;8;3;1>;
	decl model, switching, mresult = zeros(cM, 5), vp, dfunc, retval;
	decl time = timer();

	// run all experiments on the same seed
	ranseed(-1);
	
	parallel
	for (decl i = 0; i < cM; ++i)
	{
		model = new PoissonMix(data);
		switching = new MySwitching(iLineSearchMode);
		switching.SetLineSearchWarmUp(3);
		if (bExtra)
			switching.SetExtraUpdate(1);
		vp = 0.05 + ranu(1, 1) * 0.9 | ranu(2, 1) * 100;
		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt, model.UpdatePar, model.CheckPar,
			vp, constant(1e-12, 3, 1), <1;100;100> - 1e-12, 1e-12, 10000, -1);

		mresult[i][] = dfunc ~ switching.GetIterInfo();

		delete model;
		delete switching;
	}
	print("Poisson-mixture ", Switching::GetLineSearchName(iLineSearchMode), " M=", cM, bExtra ? " with extra update" : "");
	MySwitching::Report(mresult, time);

	return mresult;
}
RunCauchyExp(const iLineSearchMode, const cM, const bUsePX=TRUE, bExtra=FALSE)
{
	decl data, mresult = zeros(cM, 5), cp = 50, nu = 1, ct = 100, vp, dfunc, retval, model, switching;
	decl time = timer();

	// run all experiments on the same seed
	ranseed(-1);

	parallel
	for (decl i = 0; i < cM; ++i)
	{
		data = rant(ct, cp, nu);
		model = new MultivariateT(data, nu, bUsePX);
		switching = new MySwitching(iLineSearchMode);
		switching.SetLineSearchWarmUp(3);
		if (bExtra)
			switching.SetExtraUpdate(1);
		
		vp = meanc(data)' | vech(choleski(variance(data)));
		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt, model.UpdatePar, model.CheckPar,
			vp, <>, <>, 1e-12, 10000, -1);

		mresult[i][] = dfunc ~ switching.GetIterInfo();

		delete model;
		delete switching;
	}
	print("Multivariate-t ", Switching::GetLineSearchName(iLineSearchMode), " PX=", bUsePX, " M=", cM,
		" p=", cp, " T=", ct, " no params=", sizerc(vp), bExtra ? " with extra update" : "");
	MySwitching::Report(mresult, time);

	return mresult;
}																				
RunParafacExp(const iLineSearchMode, const cM, const cI, const cF, const dSigma=0.1, const bConcentrate=FALSE)
{
	decl data, mresult = zeros(cM, 5), vp, dfunc, retval, model, switching, ma, mb, mc, cj = cI, ck = cI, mu, mw;

	// run all experiments on the same seed
	ranseed(-1);

	ma = 3 + unit(cI, cF);  mb = 2 + unit(cj, cF);  mc = 1 + unit(ck, cF);
	[ma, mb, mc] = Parafac::Normalize(ma, mb, mc);
	data = ma * Parafac::ColKron(mc, mb)' + rann(cI, cj * ck) * dSigma;
	decl time = timer();

	parallel
	for (decl i = 0; i < cM; ++i)
	{
		model = new Parafac(data, cI, cj, ck, cF, bConcentrate);
		switching = new MySwitching(iLineSearchMode);

		vp = model.MapToPar(ma + (ranu(cI, cF) - 0.5)/10, mb + (ranu(cj, cF) - 0.5)/10, mc + (ranu(ck, cF) - 0.5)/10);
		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt, model.UpdatePar, model.CheckPar,
			vp, <>, <>, 1e-12, 100000, -1);

		mresult[i][] = dfunc ~ switching.GetIterInfo();
//		decl abc = model.MapToArgs(vp);
//		println(ma, mb, mc, Parafac::Normalize(abc[0], abc[1], abc[2]));

		delete model;
		delete switching;
	}
	print("Parafac ", Switching::GetLineSearchName(iLineSearchMode), " I=J=K=", cI, " F=", cF, " M=", cM,
		" concentrate=", bConcentrate, " sigma=", dSigma);
	MySwitching::Report(mresult, time);

	return mresult;
}

RunLowRankExp(const iLineSearchMode, const cM, const cN, const cQ, const iR, const cS)
{
	decl data, mresult = zeros(cM, 5), vp, dfunc, retval, model, switching, my, mx, mc,
		mc_start, ma_start, mb_start, mu, mv, mw;

	// run all experiments on the same seed
	ranseed(-1);

	mx = rann(cN, cQ * cS);
	mc = unit(cQ, iR) * unit(iR, cS);

	decl time = timer();

	parallel
	for (decl i = 0; i < cM; ++i)
	{
		my = mx * vec(mc) + 1 * rann(cN, 1);
		model = new LowRank(my, mx, cQ, iR, cS);
		switching = new MySwitching(iLineSearchMode);

		vp = model.StartPar();

		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt, model.UpdatePar, model.CheckPar,
			vp, <>, <>, 1e-12, 10000, -1);

		mresult[i][] = dfunc ~ switching.GetIterInfo();
//		decl abc = model.MapToArgs(vp);
//		decl mc_end = abc[0] * abc[1]';
//		println("C true=", "%6.0f", mc, "C max=", "%6.2f", mc_end, "%6.2f", mc_end_t);
//		decsvd(mc, &mu, &mw);  println("C sv:", mw);
//		decsvd(mc_end, &mu, &mw);  println("C max sv:", mw);

		delete model;
		delete switching;
	}
	print("LowRank ", Switching::GetLineSearchName(iLineSearchMode), " n=", cN, " q=s=", cQ,
		" rank=", iR, " M=", cM);
	MySwitching::Report(mresult, time);

	return mresult;
}
///////////////////////////////////////////////////////////////////////

main()
{
	format(1000);

	decl als0 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_1STEP1, Switching::LS_SQS3, Switching::LS_SQS3G, Switching::LS_PARABOLIC, Switching::LS_BRENT, Switching::LS_POWELL, Switching::LS_QSTEP};
	decl als1 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_1STEP1, Switching::LS_SQS3G, Switching::LS_PARABOLIC, Switching::LS_BRENT, Switching::LS_POWELL, Switching::LS_QSTEP};
	
	foreach (decl ls in als0)
		RunPoissonExp(ls,  5000);

	foreach (decl ls in als1)
		RunCauchyExp(ls,  1000, TRUE);

	foreach (decl ls in als1)
		RunParafacExp(ls, 100, 20, 5);

	foreach (decl ls in als1)
		RunLowRankExp(ls,  100, 200, 20, 5, 10);

exit(0);		

	// compare results for different methods
	decl als2 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_SQS3G, Switching::LS_PARABOLIC, Switching::LS_POWELL, Switching::LS_QSTEP};
	decl mlogdets = <>;
	foreach (decl ls in als2)
	{
		mlogdets ~= RunPoissonExp(ls,  5000)[][0];
//		mlogdets ~= RunCauchyExp(ls,  1000, TRUE)[][0];
//		mlogdets ~= RunParafacExp(ls, 100, 20, 5)[][0];
//		mlogdets ~= RunLowRankExp(ls,  100, 200, 20, 5, 10)[][0];
	}
	decl mgap = fabs(mlogdets - maxr(mlogdets));
	// print the largest deviation from the maximum over all methods
	println("Maximum of the likelihoods in deviation from the maximum of those achieved by each method:",
		maxc(mgap));
	println("Premature convergence, deviations:");
	println((range(0, sizer(mgap) - 1)' ~ fabs(mgap))[vecindex(maxr(mgap) .> 1e-5)][] );
	println("Premature convergence, likelihoods:");
	println("%cf", {"%8.0f","%22.10g"}, (range(0, sizer(mgap) - 1)' ~ fabs(mlogdets))[vecindex(maxr(mgap) .> 1e-5)][] );
}