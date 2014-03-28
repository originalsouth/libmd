// Returns constant zero (x must be larger than 1)
template<class number> number zero (ui i, number x)
{	switch (i)
	{	case  0 : return sqrt(x*x)-x;
		case  1 : return sqrt(x)*sqrt(x)-x;
		case  2 : return pow(pow(x,5.),.2)-x;
		case  3 : return exp(log(x))-x;
		case  4 : return log(exp(x))-x;
		case  5 : return log(pow(5.,x))/log(5.)-x;
		case  6 : return pow(5.,log(x)/log(5.))-x;
		case  7 : return sin(asin(1./x))-1./x;
		case  8 : return cos(acos(1./x))-1./x;
		case  9 : return tan(atan(x))-x;
		case 10 : return sinh(asinh(x))-x;
		case 11 : return cosh(acosh(x))-x;
		case 12 : return tanh(atanh(1./x))-1./x;
	}
	return 0;
}

template<class number> number test (ui i, number x, number y)
{	number z = x*x+y*y*y-x*y;
	switch (i)
	{	case  0 : return pow(z,3);
		case  1 : return pow(1.5,z);
		case  2 : return pow(z-5,z-4);
		case  3 : return exp(1./z);
		case  4 : return log(z);
		case  5 : return sin(z);
		case  6 : return cos(z);
		case  7 : return tan(z);
		case  8 : return sinh(z);
		case  9 : return cosh(z);
		case 10 : return tanh(z);
		case 11 : return asin(1./z);
		case 12 : return acos(1./z);
		case 13 : return atan(z);
		case 14 : return asinh(z);
		case 15 : return acosh(z);
		case 16 : return atanh(1./z);
	}
	return 42;
}

bool test_autodiff2_standard_expr()
{	duals<3> X(7,0), Y(3,1), Z(-2,2), F;
	ldf m = 0;
	for (ui t = 0; t < 14; t++)
	{	F = zero(t, X*X+Y*Y*Y-X*Y+1./Z);
		// Check if outcome is constant zero
		m = max(m, fabs(F.x));
		for (ui i = 0; i < 3; i++)
		{	m = max(m, fabs(F.dx[i]));
			for (ui j = 0; j < 3; j++)
				m = max(m, fabs(F.dxdy[i][j]));
		}
		if (m > 1e-10)
		{	printf("Part 1, test %d: expected zero, found %Lg\n", t, m);
			test_fail;
		}
	}
	
	ldf x = 3, y = -2, h = 1e-5;
	duals<2> G;
	m = 0;
	for (ui i = 0; i < 18; i++)
	{	G = test(i, duals<2>(x,0), duals<2>(y,1));
		// Compare outcome with finite difference
		m = max(m, fabs(G.x - test(i,x,y)));
		m = max(m, fabs(G.dx[0] - (test(i,x+h,y) - test(i,x-h,y)) / (2*h)));
		m = max(m, fabs(G.dx[1] - (test(i,x,y+h) - test(i,x,y-h)) / (2*h)));
		m = max(m, fabs(G.dxdy[0][0] - (test(i,x+h,y) - 2*test(i,x,y) + test(i,x-h,y)) / (h*h)));
		m = max(m, fabs(G.dxdy[1][1] - (test(i,x,y+h) - 2*test(i,x,y) + test(i,x,y-h)) / (h*h)));
		m = max(m, fabs(G.dxdy[0][1] - (test(i,x+h,y+h) - test(i,x+h,y-h) - test(i,x-h,y+h) + test(i,x-h,y-h)) / (4*h*h)));
		if (m > 1e-4)
		{	printf("part 2, test %d: expected a small number, found %Lg\n", i, m);
			test_fail;
		}
	}
	test_success;
}
