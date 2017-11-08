#include <oxstd.oxh>
#include <oxfloat.oxh>

/** Brent (1973) algorithm to find maximum of a function of one variable without using derivatives.
@param Func function to maximize: Func(x) returns the function value
@param dP0 double, first starting value
@param dP1 double, second starting value
@param dP0 one coordinate
@param dF0 double, function value at dP0
@param dF1 double, function value at dP1
@param dA,dB doubles, defining interval to maximize over. Must have dP0,dP1 in (dA,dB), but this is not checked!
@param dTol convergence tolerance
@param mxIter int, maximum no of iterations	(default 1000)
@returns {argument, function value, no of function calls, no of iterations}
@notes avoids bracketing as this makes performance worse.
*/
MaxScalarBrent(const Func, const dP0, const dP1, const dF0, const dF1, const dA, const dB, const dTol, const mxIter)
{
	decl gold = 0.5 * (3 - sqrt(5)), eps = sqrt(DBL_EPSILON), a = dA, b = dB, cfunc = 0;
	decl tol, tol2, d = 0, e = 0, x, v, w, u, fx, fv, fw, fu, p, q, r, xm;

	// NB: code minimizes, so need to change sign on function value
	v = w = x = dP0;  fw = fv = fx = -dF0;
	u = dP1;		  fu = -dF1;

//println("\nBrent tol=", dTol, "%25.15g", dF0 ~ dF1);

	for (decl itno = 0; itno < mxIter; ++itno)
	{
		// Update a, b, v, w, and x
		if (fu <= fx)
		{
			if (u < x)  b = x;
			else		a = x;
			v = w; fv = fw;  w = x; fw = fx;  x = u; fx = fu;
		}
		else
		{
			if (u < x)	a = u;
			else		b = u;
			if (fu <= fw || w == x)
			{
				v = w; fv = fw;  w = u; fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u; fv = fu;
			}
		}
		// new midpoint
		xm = 0.5 * (a + b);
		tol = eps * fabs(x) + dTol;
		tol2 = 2.0 * tol;

		// check for convergence
		if (fabs(x - xm) <= tol2 - 0.5 * (b - a))
			return {x, -fx, cfunc, itno};

		p = q = r = 0;
		if (fabs(e) > tol)
		{
			// fit parabola
			r = (x - w) * (fx - fv);  		q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;	q = 2 * (q - r);

			if (q > 0)
				p = -p;
			else
				q = -q;

			r = e; e = d;
		}

		if (fabs(p) < fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x))
		{
			// parabolic interpolation step
			d = p / q;
			u = x + d;
			// f must not be evaluated too close to a or b
			if (u - a < tol2 || b - u < tol2)
			{
				d = x < xm ? tol : -tol;
			}
		}
		else
		{
			// golden-section step
			e = x < xm ? b - x : a - x;
			d = gold * e;
		}
		// function must not be evaluated too close to x
		u = x + (fabs(d) >= tol ? d : d > 0 ? tol : -tol);
		fu = -Func(u); ++cfunc;

//println("\tBrent itno=", itno, " xm=", xm, " fu=", "%25.15g", -fu, " fm=", "%25.15g", Func(xm));
	}

	return {x, -fx, cfunc, mxIter};
}
/** Powell (1963) algorithm to find maximum of a function of one variable without using derivatives.
@param Func function to maximize: Func(x) returns the function value
@param dP0 double, first starting value
@param dP1 double, second starting value
@param dP0 one coordinate
@param dF0 double, function value at dP0
@param dF1 double, function value at dP1
@param dA,dB doubles, defining interval to maximize over. Must have dP0,dP1 in (dA,dB), but this is not checked!
@param dTol	convergence tolerance
@param mxIter int, maximum no of iterations	(default 1000)
@returns {argument, function value, no of function calls, no of iterations}
@notes avoids bracketing as this makes performance worse.
*/
MaxScalarPowell(const Func, const dP0, const dP1, const dF0, const dF1, const dA, const dB, const dTol, const mxIter)
{
	decl num, crit, order, cfunc = 0, dtolf = DBL_EPSILON * 1000;
	decl x, f, x0, x1, x2, xs, f0, f1, f2, fs, xmax, imax, move, dx12, dx20, dx01;

	// find third point
	if (dF0 < dF1)
	{
		x0 = dP0; f0 = dF0; x1 = dP1; f1 = dF1; x2 = min(dP1 + (dP1 - dP0), dB * (1 - sqrt(DBL_EPSILON)));
		f2 = Func(x2); ++cfunc;
	}
	else
	{
		x0 = max(dP0 - (dP1 - dP0), dA * (1 + sqrt(DBL_EPSILON))); x1 = dP0; f1 = dF0; x2 = dP1; f2 = dF1;
		f0 = Func(x0); ++cfunc;
	}
	x = x0 | x1 | x2;
	f = f0 | f1 | f2;
//println("\nPowell tol=", dTol, "%25.15g", dF0 ~ dF1);

	for (decl itno = 0; itno < mxIter; ++itno)
	{
		// keep the coordinates ordered
		order = sortcindex(x);
		x = x[order];
		f = f[order];
		x0 = x[0]; f0 = f[0];
		x1 = x[1]; f1 = f[1];
		x2 = x[2]; f2 = f[2];

		// convergence 1: function values all approximately the same
		if (fabs(f0 - f2) <= dtolf * fabs(f1))
		{
			imax = int(maxcindex(f));
			return {x[imax], f[imax], cfunc, itno};
		}
		// predict maximum (or minimum) xs
		dx12 = x1 - x2;  dx20 = x2 - x0;  dx01 = x0 - x1;
		num = dx12 * f0 + dx20 * f1 + dx01 * f2;

		// because we keep them ordered: sign(dx12 * dx20 * dx01) = 1, so sign(num / (dx12 * dx20 * dx01)) == sign(num)
		if (num > 0)
		{
			// predicting a new maximum at xs
			xs = 0.5 * (dx12 * (x1 + x2) * f0 + dx20 * (x2 + x0) * f1 + dx01 * (x0 + x1) * f2) / num;
			// f must not be evaluated too close to a or b
			xs = setbounds(xs, 0.5 * (dA + x0), 0.5 * (dB + x2));

			imax = int(maxcindex(f));
			xmax = x[imax];
	
			// convergence 2: prediction is close to current best coordinate
			if (fabs(xs - xmax) <= dTol)
				return {xmax, f[imax], cfunc, itno};
				
			// replace and continue
			fs = Func(xs); ++cfunc;
			f |= fs;	  // include this to determine next move
			move = 0;

			switch_single (mincindex(f))
			{	case 0: order = <3,1,2>;	// replace x0
				case 1: order = <0,3,2>;	// replace x1
				case 2: order = <0,1,3>;	// replace x2
				case 3:						// deterioration, try to create a bracket and shrink
					if (imax == 2)
					{	// get a point beyond x2, also move lhs
						xs = min(x2 - dx12 / 2, dB);
						x1 -= dx12 / 2;	 f1 = Func(x1);	++cfunc;
						order = <1,2,3>;
					}
					else if (imax == 0)
					{	// get a point before x0, also move rhs
						xs = max(x0 + dx01 / 2, dA);
						x1 += dx01 / 2;	 f1 = Func(x1);	++cfunc;
						order = <3,0,1>;
					}
					else if (imax == 1)
					{	// we have a bracket, shrink it on one side
						if (fabs(dx01) > fabs(dx12))
						{
							xs = x0 - dx01 / 2;
							order = <3,1,2>;
						}
						else
						{
							xs = x2 + dx12 / 2;
							order = <0,1,3>;
						}
					}
					fs = Func(xs); ++cfunc;
			}
		}
		else
		{	// predicted step is minimum: wrong direction, so move bracket
			if (f2 < f0)				// move bracket left
			{
				xs = max(x1 + dx12, dA);
				order = <0,1,3>;
			}
			else                		// move bracket right
			{
				xs = min(x2 - dx01, dB);
				order = <1,2,3>;
			}
			fs = Func(xs); ++cfunc;
		}

		// now select the right elements and order them
		x = (x0 | x1 | x2 | xs)[order];
		f = (f0 | f1 | f2 | fs)[order];
	}
	return {xmax, f[imax], cfunc, mxIter};
}
