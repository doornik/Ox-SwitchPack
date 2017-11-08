#include <oxstd.oxh>
#include <oxdraw.oxh>
#import <maximize>
#import "switching"


static decl s_sFolder = "./";

enum
{	BETA_SWITCHING=0, AB_SWITCHING=1
};

setFolder(sFolder)
{
	s_sFolder = sFolder;
}						   
getFileName(iSwitching, iRes, iLineSearchMode)
{
	return s_sFolder ~  "Res" ~ sprint("%d", iRes, iSwitching == 0 ? "_Beta" : "_AlphaBeta",
		"_LS", iLineSearchMode) ~ ".in7";
}
getName(iSwitching)
{
	return {"Beta","AB"}[iSwitching];
}
getName_LS(iLinesearch)
{
	switch_single (iLinesearch)
	{
		case Switching::LS_USER1:		return "LOpt";
		case Switching::LS_USER2:		return "LBeta";
		default: 			return Switching::GetLineSearchName(iLinesearch);
	}
}

loadTrace(sFile, const bDropFailures=TRUE)
{
	decl mc_trace = loadmat(sFile);
	if (mc_trace == 0)
		oxrunerror(sFile ~ " not found");

	// {"f0", "f1", "itcount", "f_end", "result", "eps", "f_0","f_1","f_2", ...}
	decl cm = sizer(mc_trace);
	println("\nLoaded:   ", sFile);
	println("M:        ", cm);
	println("Failures: ", int(sumc(mc_trace[][4] .== MAX_FUNC_FAIL)));

	decl funcfail = mc_trace[][4] .== MAX_FUNC_FAIL;
	if (bDropFailures)
		mc_trace = deleteifr(mc_trace, funcfail);
	decl finit = mc_trace[][0] .< mc_trace[][1] .? mc_trace[][0] .: mc_trace[][1];
	decl f_end = mc_trace[][3];

	// get history and transform to 0-100 scale;
	decl mf = mc_trace[][6 : ]';
	mf = (mf - minc(mf)) ./ (maxc(mf) - minc(mf));

	decl counts = mc_trace[][2]';
	
	decl vx = <0:10,20:[10]100,200:[100]1000,2000:[1000]10000>;
	return {mf, vx, counts, f_end', funcfail', cm};
}

PlotHistogram(iArea, iRes, iLinesearch1, iSwitching1, iYmax)
{
	decl mf100a, vx, cnta, f_enda, faila;
	decl sres = sprint("Res", iRes);
	decl scmp1 = getName_LS(iLinesearch1);
	decl sbeta1 = getName(iSwitching1);

	[mf100a, vx, cnta, f_enda, faila] = loadTrace(getFileName(iSwitching1, iRes, iLinesearch1));

	decl barx = <1:16> / 4;
	DrawAdjust(ADJ_AREA_X, iArea, 0, 4.1);
	DrawAdjust(ADJ_AREA_Y, iArea, 0, iYmax);
	DrawAxisAuto(iArea, TRUE);
	DrawAdjust(ADJ_AXISSCALE, AXIS_LOG10);
	DrawAdjust(ADJ_AXISLABEL, 0, 250);
	DrawAxisAuto(iArea, FALSE);
	DrawAdjust(ADJ_AXISGRID, 1);
	DrawAdjust(ADJ_AXISLABEL, 1, 250);
	DrawHistogram(iArea, countr(cnta, 10 ^ barx), 0, 0.25, 3, 15);

	if (sumr(faila .== 1))
	{
		DrawLine(iArea, 4.05, 0, 4.05, sizerc(faila), 2);
		DrawAdjust(ADJ_COLOR, 2, 12);
	}
	DrawTitle(iArea, sprint(sbeta1, " ", scmp1));
	println(sres, " ", sbeta1, " ", scmp1);
	println("%c", {"median","mean","stddev","failures"}, quantiler(cnta) ~ meanr(cnta) ~ sqrt(varr(cnta)) ~ sumr(faila .== 1));
}

PlotConvergenceRate(iRes, iLinesearch1, iSwitching1, iLinesearch2, iSwitching2)
{
	decl mf100a, mf100a1, mf100b, mf100b1, vxa, vxb;
	decl sres = sprint("Res", iRes);
	decl scmp1 = getName_LS(iLinesearch1), scmp2 = getName_LS(iLinesearch2);
	decl stxt1 = getName_LS(iLinesearch1), stxt2 = getName_LS(iLinesearch2);
	decl sbeta1 = getName(iSwitching1), sbeta2 = getName(iSwitching2);

	SetDrawWindow(sres ~ "_" ~ sbeta1 ~ scmp1 ~ "_" ~ sbeta2 ~ scmp2);
	DrawAdjust(ADJ_AREA_X, 0, 0, 4.1);
	DrawAdjust(ADJ_AREA_Y, 0, 0, 1);
	DrawAxisAuto(0, TRUE);
	DrawAdjust(ADJ_AXISSCALE, AXIS_LOG10);
	DrawAdjust(ADJ_AXISLABEL, 0, 250);
	DrawAxisAuto(0, FALSE);
	DrawAdjust(ADJ_AXISLABEL, 1, 250);

	[mf100a, vxa] = loadTrace(getFileName(iSwitching1, iRes, iLinesearch1));
	// drop starting value, set missing values (at end) to 1
	mf100a = mf100a[][1 : ];  vxa = vxa[][1 : ];
	mf100a1	= mf100a .== .NaN .? 1 .: mf100a;
	
	DrawTitle(0, sres);
	DrawXMatrix(0, quantilec(mf100a1', <0.5>), "Median " ~ stxt1, log10(vxa), 0, 0, 3);
	DrawAdjust(ADJ_COLOR, 3, 7);
	DrawXMatrix(0, quantilec(mf100a1', <0.1>), "90% " ~ stxt1, log10(vxa), 0, 0, 3);
	DrawAdjust(ADJ_COLOR, 3, 6);
	DrawXMatrix(0, quantilec(mf100a1', <0.01>), "99% " ~ stxt1, log10(vxa), 0, 0, 3);
	DrawAdjust(ADJ_COLOR, 3, 5);

	[mf100b, vxb] = loadTrace(getFileName(iSwitching2, iRes, iLinesearch2));
	// drop starting value, set missing values (at end) to 1
	mf100b = mf100b[][1 : ];  vxb = vxb[][1 : ];
	mf100b1	= mf100b .== .NaN .? 1 .: mf100b;
	
	DrawXMatrix(0, quantilec(mf100b1', <0.5>), "Median " ~ stxt2, log10(vxb), 0, 0, 8);
	DrawAdjust(ADJ_COLOR, 8, 10);
	DrawXMatrix(0, quantilec(mf100b1', <0.1>), "90% " ~ stxt2, log10(vxb), 0, 0, 8);
	DrawAdjust(ADJ_COLOR, 8, 9);
	DrawXMatrix(0, quantilec(mf100b1', <0.01>), "99% " ~ stxt2, log10(vxb), 0, 0, 8);
	DrawAdjust(ADJ_COLOR, 8, 8);

	ShowDrawWindow();
}
PlotLogLik(iArea, iRes, iLinesearch1, iSwitching1, iLinesearch2, iSwitching2, bPlotDiff=TRUE)
{
	decl mf100a, mf100a1, mf100b, mf100b1, vx, cnta, cntb, f_enda, f_endb, faila, failb;
	decl sbeta1 = getName(iSwitching1), sbeta2 = getName(iSwitching2);
	decl scmp1 = getName_LS(iLinesearch1), scmp2 = getName_LS(iLinesearch2);
	decl ct = 119;

	[mf100a, vx, cnta, f_enda, faila] = loadTrace(getFileName(iSwitching1, iRes, iLinesearch1), FALSE);
	[mf100b, vx, cntb, f_endb, failb] = loadTrace(getFileName(iSwitching2, iRes, iLinesearch2), FALSE);

	// omit if either failed
	f_enda[vecindex(faila + failb)] = .NaN;
	f_endb[vecindex(faila + failb)] = .NaN;
//	f_endb = deleteifc(f_endb, faila + failb);
	f_enda *= (ct / 2);
	f_endb *= (ct / 2);

	if (bPlotDiff)
	{
		DrawTMatrix(iArea, f_enda - f_endb, sprint(sbeta1, scmp1, "-", sbeta2, scmp2));
		DrawAdjust(ADJ_INDEX, 1);
		if (any(faila .!= 0))
			DrawTMatrix(iArea, faila .!= 0 .? -1 .: .NaN, "", 1, 1, 1, 1, 3);
		if (any(failb .!= 0))
			DrawTMatrix(iArea, failb .!= 0 .? 1 .: .NaN, "", 1, 1, 1, 1, 3);
		if (max(fabs(f_enda - f_endb)) < 0.2)
			DrawAdjust(ADJ_AREA_Y, iArea, -0.201, 0.201);
	}
	else
		DrawTMatrix(iArea, f_enda | f_endb, {sbeta1 ~ scmp1, sbeta2 ~ scmp2});
}
Plot_I1_Danish(iRes)
{
	setFolder("trace/");
	
	//////////////////////////////////////
	PlotConvergenceRate( iRes, Switching::LS_USER2, AB_SWITCHING,   Switching::LS_1STEP, AB_SWITCHING);
	PlotConvergenceRate( iRes, Switching::LS_USER2, AB_SWITCHING,   Switching::LS_NONE, AB_SWITCHING);
	PlotConvergenceRate( iRes, Switching::LS_USER2, BETA_SWITCHING, Switching::LS_NONE, BETA_SWITCHING);

	decl fn_histall = [=](iRes, iYmax)
	{
		SetDrawWindow(sprint("Res", iRes, "_rate"));
		PlotHistogram(0, iRes, Switching::LS_NONE  , BETA_SWITCHING, iYmax);
		PlotHistogram(1, iRes, Switching::LS_USER2,  BETA_SWITCHING, iYmax);
		PlotHistogram(2, iRes, Switching::LS_NONE  , AB_SWITCHING ,  iYmax);
		PlotHistogram(3, iRes, Switching::LS_1STEP,  AB_SWITCHING ,  iYmax);
		PlotHistogram(4, iRes, Switching::LS_1STD,   AB_SWITCHING ,  iYmax);
		PlotHistogram(5, iRes, Switching::LS_USER2,  AB_SWITCHING ,  iYmax);
		PlotHistogram(6, iRes, Switching::LS_USER1,  AB_SWITCHING ,  iYmax);
		DrawAdjust(ADJ_AREAMATRIX, 2, 4);
		ShowDrawWindow();
		
		SetDrawWindow(sprint("Res", iRes, "_fdiff"));
		PlotLogLik(0, iRes, Switching::LS_NONE  , BETA_SWITCHING, Switching::LS_USER2, BETA_SWITCHING);
		PlotLogLik(1, iRes, Switching::LS_USER2,  BETA_SWITCHING, Switching::LS_USER2, AB_SWITCHING);
		PlotLogLik(2, iRes, Switching::LS_NONE  , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING);
		PlotLogLik(3, iRes, Switching::LS_1STEP,  AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING);
		PlotLogLik(4, iRes, Switching::LS_USER2 , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING);
		PlotLogLik(5, iRes, Switching::LS_USER1 , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING);
		DrawAdjust(ADJ_AREAMATRIX, 2, 3);
		ShowDrawWindow();
		
		SetDrawWindow(sprint("Res", iRes, "_loglik"));
		PlotLogLik(0, iRes, Switching::LS_NONE  , BETA_SWITCHING, Switching::LS_USER2, BETA_SWITCHING, FALSE);
		PlotLogLik(1, iRes, Switching::LS_USER2,  BETA_SWITCHING, Switching::LS_USER2, AB_SWITCHING, FALSE);
		PlotLogLik(2, iRes, Switching::LS_NONE  , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING, FALSE);
		PlotLogLik(3, iRes, Switching::LS_1STEP,  AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING, FALSE);
		PlotLogLik(4, iRes, Switching::LS_USER2 , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING, FALSE);
		PlotLogLik(5, iRes, Switching::LS_USER1 , AB_SWITCHING ,  Switching::LS_1STEP, AB_SWITCHING, FALSE);
		DrawAdjust(ADJ_AREAMATRIX, 2, 3);
		ShowDrawWindow();
	};
	fn_histall(iRes, 600);
}
main()
{
//	Plot_I1_Danish(6);
	Plot_I1_Danish(5);
}