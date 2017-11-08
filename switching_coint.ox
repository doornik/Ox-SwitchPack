#include <oxstd.oxh>
#include <oxprob.oxh>

#import <database>
#import "switching"
#import "coint1"


///////////////////////////////////////////////////////////////////////
class MySwitching : Switching
{
	MySwitching(iLineSearchMode);
	static Report(mResult, time);
	static SaveTrace(aTrace, iRes, bBetaSwitching, iLineSearchMode);
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
	if (mResult == <>)
		return;
	mResult[][0] -= max(mResult[][0]);
	println("%cf", {"%10.0f","%10.0f","%10.0f","%10.0f","%10.0f","%10.2f"}, "%c", {"iters","updates","logliks","failures","below","CPU(s)"},
		meanc(mResult[][2:]) ~ sumc(mResult[][1] .> MAX_WEAK_CONV) ~ .NaN ~ (timer() - time) / 100);
}
MySwitching::SaveTrace(aTrace, iRes, bBetaSwitching, iLineSearchMode)
{
	decl cm = sizeof(aTrace), i;
	
	if (cm == 0)
		return;

	// flatten the trace
	decl mc_trace = nans(cm, 45);

	foreach (decl maxinfo in aTrace[i])
		mc_trace[i][] = maxinfo[ : 1][1]' ~ maxinfo[sizer(maxinfo) - 1][] ~ maxinfo[][1]';

	// save the trace
	decl sfile = "trace/Res" ~ sprint("%d", iRes, bBetaSwitching ? "_Beta" : "_AlphaBeta",
		"_LS", iLineSearchMode) ~ ".in7";
	
	savemat(sfile, mc_trace,
		{   "f0","f1","itcount","f_end","result","eps",
			"f_0",
			"f_1",   "f_2",   "f_3",   "f_4",   "f_5",   "f_6",   "f_7",   "f_8",   "f_9",
			"f_10",  "f_20",  "f_30",  "f_40",  "f_50",  "f_60",  "f_70",  "f_80",  "f_90",
			"f_100", "f_200", "f_300", "f_400", "f_500", "f_600", "f_700", "f_800", "f_900",
			"f_1000","f_2000","f_3000","f_4000","f_5000","f_6000","f_7000","f_8000","f_9000",
			"f_10000"
		});

	println("---- Results saved to: ", sfile);

}
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
RunCointExp(iLineSearchMode, const cM, const iRes=6, const bBetaSwitching=FALSE, dEps=1e-12, bTrace=FALSE)
{
	decl model, switching, mresult = zeros(cM, 5), atrace = new array[cM], vp, dfunc, retval,
		my, amg, amh, alpha0, beta0, chol0;
	decl time = timer();

	// load data
	decl db = new Database();
	db.Load("Danish_data_Juselius(2006).xlsx");
	my = db.GetVar({"Lm3r","Lyr","Dpy","Rm","Rb"});// Y_VAR: dependent variables
	delete db;

	decl H1 = unit(6)[][:2], H2 = unit(6)[][<0,5>], H3 = sumr(unit(6)[][<0,2,4>]);
	decl H6, H7, H8, A = <-1,0,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,1,0,0>;

	H6 = unit(6)[][0:1];	H6[2][1] = -1;
	H7 = unit(6)[][3:5];
	H7[4:5][2] = 1;
	H8 = unit(6)[][4:5];
	H8[2][1] = -500;
	H8[3][0] = 1;

	if (bBetaSwitching && iLineSearchMode > Switching::LS_NONE && iLineSearchMode < Switching::LS_USER2)
	{
		iLineSearchMode = Switching::LS_USER2;
	}

	switch_single (iRes)
	{
		case 5:	 amg = {};	amh = {H1, H2, unit(6)[][2 : ]};
		case 6:	 amg = {};	amh = {H1, H2, H3};
		case 7:	 amg = {unit(5)[][:3]};	amh = {H1, H2, H3};
		case 8:	 amg = {unit(5)[][1:]};	amh = {H1, H2, H3};
		case 9:  amg = {};	amh = {H6, H7, H8};
		case 10: amg = {A};	amh = {H6, H7, H8};
		case 11: amg = {unit(5)[][:2],unit(5)[][:3],unit(5)[][2:]};	amh = {H6, H7, H8};
	}

	println("\n\n", iRes, " ######### rank ", 3);

	// first run on actual data, then run on simulated data
	decl model0 = new Coint(my);
	vp = model0.SetRank(3);
	[alpha0, beta0]	= model0.MapToArgs(vp);
	vp = model0.Restrict(amg, amh);
	[alpha0, beta0]	= model0.MapToArgs(vp);
	chol0 = choleski(model0.OmegaAt(alpha0, beta0));

	switching = new MySwitching(iLineSearchMode);
	[vp] = switching.Estimate(model0.LogLikAt,
		bBetaSwitching ? model0.UpdatePar_BetaSwitching : model0.UpdatePar, model0.CheckPar,
		vp, <>, <>, dEps, 100000, 2,
		iLineSearchMode == Switching::LS_USER1 ? model0.LineSearchOpt : iLineSearchMode == Switching::LS_USER2 ? model0.LineSearch1Beta : model0.LineSearchQBeta);
 	delete switching;

	// run all experiments on the same seed
	ranseed(-1);
	
	parallel
	for (decl i = 0; i < cM; ++i)
	{
		model = new Coint(model0.GenerateVAR1(rann(sizer(my), sizec(my)) * chol0', alpha0, beta0));
		model.SetRank(3);
		vp = model.Restrict(amg, amh);
		
		switching = new MySwitching(iLineSearchMode);
		[vp, dfunc, retval] = switching.Estimate(model.LogLikAt,
			bBetaSwitching ? model.UpdatePar_BetaSwitching : model.UpdatePar, model.CheckPar,
			vp, <>, <>, dEps, 10000, -1,
			iLineSearchMode == Switching::LS_USER1 ? model.LineSearchOpt : iLineSearchMode == Switching::LS_USER2 ? model.LineSearch1Beta : model.LineSearchQBeta);

		mresult[i][] = dfunc ~ switching.GetIterInfo();
		if (bTrace)
			atrace[i] = switching.GetTrace();

		delete model;
		delete switching;
	}
	print("Cointegration ", Switching::GetLineSearchName(iLineSearchMode), " M=", cM, " Res=", iRes, " eps=", dEps);
	if (iLineSearchMode == Switching::LS_USER1)
		print(" Least squares line search");
	else if (iLineSearchMode == Switching::LS_USER2)
		print(" 1Step Beta line search");
	else if (iLineSearchMode == Switching::LS_USER3)
		print(" QStep Beta line search");
	print(bBetaSwitching ? ", beta switching" : ", alpha-beta switching");
	MySwitching::Report(mresult, time);
	if (bTrace)
		MySwitching::SaveTrace(atrace, iRes, bBetaSwitching, iLineSearchMode);

	return mresult;
}
///////////////////////////////////////////////////////////////////////

main()
{
	format(1000);

	// restrictions
    // 	in paper:          Aa    Ab    Bb    Cb    Ac    Dc    
    //  in code:            5     6     7     8     9    10         
    //  in trace:         Res5  Res6  Res7  Res8  Res9  Res10

	// Run in parallel (default) or with -rp1 for single process

	decl als0 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_1STEP1, Switching::LS_USER1, Switching::LS_USER2};
	decl als1 = {Switching::LS_NONE, Switching::LS_USER2};
	
	// initial estimation only with more detailed reports for the first table of results
//	foreach (decl ls in als0)
//		for (decl icase = 5; icase <= 11; ++icase)
//			RunCointExp(ls, 0, icase, FALSE);
//
//	foreach (decl ls in als1)
//		for (decl icase = 5; icase <= 10; ++icase)
//			RunCointExp(ls, 0, icase, TRUE);

	// simulations, case 6 and 10
	decl cm = 1000, btrace = FALSE, cases = {6, 10}; //{5, 6, 7, 8, 9, 10};

	// comment the next two lines in to create the longer trace	for Aa (figs 2,3)
//	btrace = TRUE;
//	cases = {5};
	
	Switching::SetTraceProgress(btrace);
	foreach (decl icase in cases)
	{
		decl als2 = {Switching::LS_NONE, Switching::LS_USER2, Switching::LS_USER3};
		foreach (decl ls in als2)
			RunCointExp(ls, cm, icase, TRUE, 1e-12, btrace);

		decl als3 = {Switching::LS_NONE, Switching::LS_1STEP, Switching::LS_1STD, Switching::LS_USER1, Switching::LS_USER2};
		foreach (decl ls in als3)
			RunCointExp(ls, cm, icase, FALSE, 1e-12, btrace);
	}
	foreach (decl icase in cases)
	{
		decl als3 = {Switching::LS_QSTEP, Switching::LS_USER2, Switching::LS_USER3, Switching::LS_BRENT, Switching::LS_BRENTSTD, Switching::LS_POWELL};
		foreach (decl ls in als3)
			RunCointExp(ls, cm, icase, FALSE, 1e-12, btrace);
	}
}