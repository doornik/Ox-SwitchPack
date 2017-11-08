#include <oxstd.oxh>
#import "maxscalar"


Func1(const vP)
{
	return -(sqr(vP) + 2 * exp(vP));
}
Func2(const vP)
{
	return -100.0 * sqr(vP^2 - 1.0) - (1.0 - vP^2)^2;
}
Func3(const vP)
{
	return -sqr(vP^2 - 1.0) - sqr(1.0 - vP^2);
}
Func4(const vP)
{
	return -100.0 * sqr(vP^3 - 1.0)^2 - (1.0 - vP^2)^2;
}
Func20(const vP)
{
	decl i = range(1, 20);
	return -double(sumr(sqr((2.0 * i - 5) ./ (vP - sqr(i)))));
}

Test(const func, const dP0, const dP1, const dA, const dB, const dTol)
{
	decl retval = MaxScalarBrent(func, dP0, dP1, func(dP0), func(dP1), dA, dB, dTol);
	print("Brent ", "%#25.15g", retval[0], "%#25.15g", retval[1], " cfunc=", retval[2]);
	retval = MaxScalarPowell(func, dP0, dP1, func(dP0), func(dP1), dA, dB, dTol);
	println(" Powell ", "%#25.15g", retval[0], "%#25.15g", retval[1], " cfunc=", retval[2]);
}

main()
{
	decl vp0, vp1;

	Test(Func1, 0.9, 1.1, -10, 10, 1e-10);
	Test(Func1, 0, 2, -10, 10, 1e-10);

	Test(Func2, -1.4, -1.3, -20, 20, 1e-10);
	Test(Func3, -2.0, -1.8, -20, 20, 1e-10);
	Test(Func4, 22.0, 16, -2000, 2000, 1e-10);

//	decl i, d, a, b, eps = 1e-1;
//	for (i = 1; i < 20; ++i)
//	{
//		d = double(i);
//		a = i^2, b = (1 + i)^2;
//		vp0 = 0.5 * (a + b);
//		vp1 = a + 0.75 * (b - a);
//		Test(Func20, vp0, vp1, a, b, eps);
//	}
}