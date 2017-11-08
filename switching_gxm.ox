#include <oxstd.oxh>
#include <oxfloat.oxh>
#import <switching>


///////////////////////////////////////////////////////////////////////
//
class GaussianMix
{
	GaussianMix(const mData, const cS);
	static Generate(const cN, const vProb, const avMu, const avSD);
	static MapToArgs(const vP, const cS, const cP);
	static MapToPar(const vProb, const avMu, const avSD);
	GetStartPar();
	GetLowerBound();
	GetUpperBound();
	GetLogGaussian(const vP);
	virtual UpdatePar(vP);
	virtual LogLikAt(const vP);
	CheckPar(const vP);

	decl m_mData;
	decl m_cS, m_cN, m_cP;
};

GaussianMix::GaussianMix(const mData, const cS)
{
	m_mData = mData;
	m_cS = cS;				// no of mixtures
	m_cN = sizer(mData);	// no of observations
	m_cP = sizec(mData);	// dimension of multivariate normal
}
/** Generates a random sample.
@param cN no of observations N
@param vProb[S] probabilities of each state
@param avMu[S] array[S] with means for each state (for p-dimensional normal: matrix[1][p] with means for each state)
@param avSd[S] array[S] with standard deviations for each state (for p-dimensional normal: matrix[1][p] with standard
deviations for each state, so assuming diagonal variance matrix)
@returns matrix[N][p]
*/
GaussianMix::Generate(const cN, const vProb, const avMu, const avSD)
{
	decl cs = sizerc(vProb), cp = sizerc(avMu[0]);
	// draw the states: this is a unit vector for each, e.g. 0100
	decl ms = countr(ranu(cN, 1), cumulate(vec(vProb[ : cs - 2])));
	decl mx = rann(cN, cp);
	// loop over the states, drawing from the appropriate normal distribution
	
	for (decl i = 0; i < cs; ++i)
	{
		decl state_i = vecindex(ms[][i]);
		mx[state_i][] = avMu[i] + mx[state_i][] .* avSD[i];
	}
	return mx;
}
/** Split vectorized parameters.
@returns array[3]: {probs[1][S], {mu_0[1][p],...}[S], {sd_0[1][p],...}[S] }.
*/
GaussianMix::MapToArgs(const vP, const cS, const cP)
{
	decl prob = vP[ : cS - 2]', amu, asd;
	amu = asd = new array[cS];
	
	for (decl i = 0, j = cS - 1; i < cS; ++i, j += 2 * cP)
	{
		amu[i] = vP[j : j + cP - 1]';
		asd[i] = vP[j + cP : j + 2 * cP - 1]';
	}
	prob ~= 1 - sumr(prob);
	return { prob, amu, asd };
}
/**	Map separate parameters to vectorized version.
@param vProb[S] probabilities of each state
@param avMu[S] array[S] with means for each state (for p-dimensional normal: matrix[1][p] with means for each state)
@param avSd[S] array[S] with standard deviations for each state (for p-dimensional normal: matrix[1][p] with standard
deviations for each state, so assuming diagonal variance matrix)
@returns vectorized parameters: S-1 probabilities, Means 1, SDs 1, Means 2, SDs 2, ...
*/
GaussianMix::MapToPar(const vProb, const avMu, const avSD)
{
	decl vp = vec(vProb[ : sizerc(vProb) - 2]);
	for (decl i = 0; i < sizeof(avMu); ++i)
	{
		vp |= vec(avMu[i]) | vec(avSD[i]);
	}
	return vp;
}
/**	Get the lower bound of the vectorized coefficients
*/
GaussianMix::GetLowerBound()
{
	decl vp = zeros(m_cS - 1, 1);
	for (decl i = 0; i < m_cS; ++i)
	{
		vp |= constant(-.Inf, m_cP, 1) | constant(1e-10, m_cP, 1);
	}
	return vp;
}
/**	Get the lower bound of the vectorized coefficients
*/
GaussianMix::GetUpperBound()
{
	decl vp = ones(m_cS - 1, 1);
	for (decl i = 0; i < m_cS; ++i)
	{
		vp |= constant(+.Inf, m_cP, 1) | constant(+.Inf, m_cP, 1);
	}
	return vp;
}
/**	Determines neutral starting values.
*/
GaussianMix::GetStartPar()
{
	decl prob = constant(1 / m_cS, 1, m_cS), amu, asd;
	amu = asd = new array[m_cS];

	// get distance measure of standardized data on [0,1] scale
	decl dist = sqrt(sumsqrr(standardize(m_mData)) / m_cP);
	dist -= min(dist);
	dist /= max(dist);
	
	// use this for initial partition
	for (decl i = 0, p_lo = 0, p_hi = 0; i < m_cS; ++i)
	{
		p_lo = p_hi;  p_hi = p_lo + prob[i];
		decl sel = vecindex(dist .>= p_lo .&& dist .<= p_hi);
		// need protection against selection getting too small?
		prob[i] = sizerc(sel) / sizer(m_mData);
		amu[i] = meanc(m_mData[sel][]);
		asd[i] = sqrt(varc(m_mData[sel][]));
	}
	return MapToPar(prob, amu, asd);
}
/**	Selects parameters for convergence check.
*/
GaussianMix::CheckPar(const vP)
{
	return vP;
}
/**	Returns the log of the multivariate normal distribution at vP.
*/
GaussianMix::GetLogGaussian(const vP)
{
	decl prob, amu, asd;
	[prob, amu, asd] = MapToArgs(vP, m_cS, m_cP);

	decl vlogdet = zeros(1, m_cS), mxtx = zeros(m_cN, m_cS);
	for (decl i = 0; i < m_cS; ++i)
	{
		vlogdet[i] = 2 * sumr(log(asd[i]));
		mxtx[][i] = sumsqrr((m_mData - amu[i]) ./ asd[i]);
	}
	return -0.5 * (mxtx + vlogdet + m_cP * log(M_2PI));
}
GaussianMix::UpdatePar(vP)
{
	decl mxtx = GetLogGaussian(vP), prob = vP[ : m_cS - 2]', mu, sd, vp_prev = vP;
	prob ~= 1 - sumr(prob);
	decl mxr_xtx = maxr(mxtx);
	decl posterior = prob .* exp(mxtx - mxr_xtx);
	posterior ./= sumr(posterior);

	prob = meanc(posterior);		   // new probabilities
	vP[ : m_cS - 2] = prob';
	posterior ./= prob;
	for (decl i = 0, j = m_cS - 1; i < m_cS; ++i, j += 2 * m_cP)
	{
		mu = meanc(posterior[][i] .* m_mData);
		sd = sqrt( meanc(posterior[][i] .* sqr(m_mData - mu)));
		vP[j : ] = mu';
		vP[j + m_cP : ] = sd';
	}
//println(vp_prev ~ vP);	
	return vP;
}
/**	Returns the loglikelihood.
*/
GaussianMix::LogLikAt(const vP)
{
	decl mxtx = GetLogGaussian(vP), prob = vP[ : m_cS - 2]';
	prob ~= 1 - sumr(prob);
	decl mxr_xtx = maxr(mxtx);
	decl mloglik = exp(mxtx - mxr_xtx) * prob';
	return double(meanc(mxr_xtx + log(mloglik)));
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
//	SetExtraUpdate(1);
}
MySwitching::Report(mResult, time)
{
	mResult[][0] -= max(mResult[][0]);
	println("%cf", {"%10.0f","%10.0f","%10.0f","%10.0f","%10.0f","%10.2f"}, "%c", {"iters","updates","logliks","failures","below","CPU(s)"},
		meanc(mResult[][2:]) ~ sumc(mResult[][1] .> MAX_WEAK_CONV) ~ .NaN ~ (timer() - time) / 100);
}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
RunGXM(const iLineSearchMode, const cM, const cN, const cS, const cP)
{
	decl vpi, amu, asd, mu;
	decl model, switching, mresult = zeros(cM, 5), vp_dgp, vp_0, vp, dfunc, retval;

	decl fn_random_args = [=]()
	{
		decl vprob = zeros(1, cS), amu, asd;
		amu = asd = new array[cS];
		for (decl i = 0; i < cS; ++i)
		{
			vprob[i] = ranu(1, 1) * (1 - sumr(vprob));
			amu[i] = i + rann(1, cP) / 10;
			asd[i] = 0.5 + ranu(1, cP);
		}
		vprob[cS - 1] = 1 - sumr(vprob[ : cS - 2]);
		return {vprob, amu, asd};
	};

	decl time = timer();

	// run all experiments on the same seed
	ranseed(-1);
	
//	parallel
	for (decl i = 0; i < cM; ++i)
	{
		// draw a design
		[vpi, amu, asd] = fn_random_args();
		vp_dgp = GaussianMix::MapToPar(vpi, amu, asd);
		// draw a sample and create model
		model = new GaussianMix(GaussianMix::Generate(cN, vpi, amu, asd), cS);
		switching = new MySwitching(iLineSearchMode);

		vp = vp_0 = model.GetStartPar();
//println(vpi, amu, asd, "%25.15g", model.LogLikAt(vp), vp_0 ~ vp_dgp ~ vp);		
//println("%25.15g", model.LogLikAt(vp), vp_0 ~ vp_dgp ~ vp ~ model.UpdatePar(vp));		
//exit(0);
		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt, model.UpdatePar, model.CheckPar,
			vp, model.GetLowerBound(), model.GetUpperBound(), 1e-12, 10000, -1);

//println("%25.15g", vp_0 ~ vp_dgp ~ vp);		
		mresult[i][] = dfunc ~ switching.GetIterInfo();

		delete model;
		delete switching;
	}
	print("Gaussian mixture model ", Switching::GetLineSearchName(iLineSearchMode), " M=", cM, " N=", cN, " S=", cS, " P=", cP);
	MySwitching::Report(mresult, time);

	return mresult;
}
///////////////////////////////////////////////////////////////////////

main()
{
	format(1000);

	decl als0 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_SQS3G, Switching::LS_PARABOLIC};
	decl mlogdets = <>;

	// 1. S=3,P=1: methods converge to different stationary points
	// 2. S=3,P=1: why does parabolic have failures?
	foreach (decl ls in als0)
		mlogdets ~=	RunGXM(ls, 100, 1000, 3, 1)[][0];

	decl mgap = fabs(mlogdets - maxr(mlogdets));
	// print the largest deviation from the maximum over all methods
	println("\nMaximum of the likelihoods in deviation from the maximum of those achieved by each method:",
		maxc(mgap));
	println("Premature convergence, deviations:");
	println((range(0, sizer(mgap) - 1)' ~ fabs(mgap))[vecindex(maxr(mgap) .> 1e-5)][] );
	println("Premature convergence, likelihoods:");
	println("%cf", {"%8.0f","%22.10g"}, (range(0, sizer(mgap) - 1)' ~ fabs(mlogdets))[vecindex(maxr(mgap) .> 1e-5)][] );
}