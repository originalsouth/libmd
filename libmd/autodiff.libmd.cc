// Single differentiation

struct var
{	ldf val, x;
	var() {}
	var (ldf f, ldf fx = 1) : val(f), x(fx) {}
	ldf diff() { return x; }
};

var operator+ (const var& F, const var& G) { return var(F.val+G.val, F.x+G.x); }
var operator- (const var& F, const var& G) { return var(F.val-G.val, F.x-G.x); }
var operator- (const var& F) { return var(-F.val, -F.x); }
var operator* (const var& F, const var& G) { return var(F.val*G.val, F.x*G.val + F.val*G.x); }
var operator/ (const var& F, const var& G) { return var(F.val/G.val, (F.x*G.val - F.val*G.x) / pow(G.val,2)); }

template<class number> var operator+ (const var& F, number a) { return var(F.val+a, F.x); }
template<class number> var operator- (const var& F, number a) { return var(F.val-a, F.x); }
template<class number> var operator* (const var& F, number a) { return var(F.val*a, F.x*a); }
template<class number> var operator/ (const var& F, number a) { return var(F.val/a, F.x/a); }

template<class number> var operator+ (number a, const var& F) { return F+a; }
template<class number> var operator- (number a, const var& F) { return (-F)+a; }
template<class number> var operator* (number a, const var& F) { return F*a; }
template<class number> var operator/ (number a, const var& F)
{ return var(a/F.val, -a*F.x/pow(F.val,2));
}

var sqrt (const var& F)
{ return var(sqrt(F.val), F.x / 2 / sqrt(F.val));
}
template<class number> var pow (const var& F, number n)
{ return var(pow(F.val, n), F.x * n * pow(F.val, n-1));
}
var exp (const var& F)
{ return var(exp(F.val), F.x * exp(F.val));
}
template<class number> var pow (number a, const var& G)
{ return var(pow(a, G.val), G.x * log(a) * pow(a, G.val));
}
var pow (const var& F, const var& G)
{ return var(pow(F.val, G.val), (G.x*log(F.val) + F.x*G.val/F.val) * pow(F.val, G.val));
}
var log (const var& F)
{ return var(log(F.val), F.x / F.val);
}
var sin (const var& F)
{ return var(sin(F.val), F.x * cos(F.val));
}
var cos (const var& F)
{ return var(cos(F.val), -F.x * sin(F.val));
}
var tan (const var& F) { return sin(F)/cos(F); }
var sinh (const var& F)
{ return var(sinh(F.val), F.x * cosh(F.val));
}
var cosh (const var& F)
{ return var(cosh(F.val), F.x * sinh(F.val));
}
var tanh (const var& F) { return sinh(F)/cosh(F); }



// Double differentiation

struct var2
{	ldf val, x, y, xy;
	var2() {}
	var2 (ldf f, ldf fx = 1, ldf fy = 1, ldf fxy = 0) : val(f), x(fx), y(fy), xy(fxy) {}
	ldf diff1() { return x; }
	ldf diff2() { return y; }
	ldf diff12() { return xy; }
};

var2 operator+ (const var2& F, const var2& G) { return var2(F.val+G.val, F.x+G.x, F.y+G.y, F.xy+G.xy); }
var2 operator- (const var2& F, const var2& G) { return var2(F.val-G.val, F.x-G.x, F.y-G.y, F.xy-G.xy); }
var2 operator- (const var2& F) { return var2(-F.val, -F.x, -F.y, -F.xy); }
var2 operator* (const var2& F, const var2& G)
{ return var2(F.val*G.val,
              F.val*G.x + F.x*G.val,
              F.val*G.y + F.y*G.val,
              F.val*G.xy + F.x*G.y + F.y*G.x + F.xy*G.val);
}
var2 operator/ (const var2& F, const var2& G)
{ return var2(F.val/G.val,
              (F.x*G.val - F.val*G.x) / pow(G.val,2),
              (F.y*G.val - F.val*G.y) / pow(G.val,2),
              (2*F.val*G.x*G.y - G.val*(F.val*G.xy + F.x*G.y + F.y*G.x - G.val*F.xy)) / pow(G.val,3));
}

template<class number> var2 operator+ (const var2& F, number a) { return var2(F.val+a, F.x, F.y, F.xy); }
template<class number> var2 operator- (const var2& F, number a) { return var2(F.val-a, F.x, F.y, F.xy); }
template<class number> var2 operator* (const var2& F, number a) { return var2(F.val*a, F.x*a, F.y*a, F.xy*a); }
template<class number> var2 operator/ (const var2& F, number a) { return var2(F.val/a, F.x/a, F.y/a, F.xy/a); }

template<class number> var2 operator+ (number a, const var2& F) { return F+a; }
template<class number> var2 operator- (number a, const var2& F) { return (-F)+a; }
template<class number> var2 operator* (number a, const var2& F) { return F*a; }
template<class number> var2 operator/ (number a, const var2& F)
{ return var2(a/F.val,
              -a*F.x/pow(F.val,2),
              -a*F.y/pow(F.val,2),
              a*(2*F.x*F.y-F.val*F.xy)/pow(F.val,3));
}

var2 sqrt (const var2& F)
{ return var2(sqrt(F.val),
              F.x / 2 / sqrt(F.val),
              F.y / 2 / sqrt(F.val),
              (F.xy - F.x*F.y/2/F.val) / 2 / sqrt(F.val));
}
template<class number> var2 pow (const var2& F, number n)
{ return var2(pow(F.val, n),
              F.x * n * pow(F.val, n-1),
              F.y * n * pow(F.val, n-1),
              n * (F.val*F.xy + (n-1)*F.x*F.y) * pow(F.val, n-2));
}
var2 exp (const var2& F)
{ ldf g = exp(F.val);
  return var2(g,
              F.x * g,
              F.y * g,
              (F.xy + F.x*F.y) * g);
}
template<class number> var2 pow (number a, const var2& G)
{ ldf z = pow(a, G.val), la = log(a);
  return var2(z,
              G.x * la * z,
              G.y * la * z,
              (G.xy + la*G.x*G.y) * la * z);
}
var2 pow (const var2& F, const var2& G)
{ ldf z = pow(F.val, G.val), lf = log(F.val);
  return var2(z,
              (G.x*lf + F.x*G.val/F.val) * z,
              (G.y*lf + F.y*G.val/F.val) * z,
              ((G.x*lf + F.x*G.val/F.val)*(G.y*lf + F.y*G.val/F.val) + G.xy*lf + (F.x*G.y+F.y*G.x+F.xy*G.val-F.x*F.y*G.val/F.val)/F.val) * z);
}
var2 log (const var2& F)
{ return var2(log(F.val),
              F.x / F.val,
              F.y / F.val,
              (F.xy*F.val - F.x*F.y) / pow(F.val,2));
}
var2 sin (const var2& F)
{ return var2(sin(F.val),
              F.x * cos(F.val),
              F.y * cos(F.val),
              F.xy * cos(F.val) - F.x*F.y * sin(F.val));
}
var2 cos (const var2& F)
{ return var2(cos(F.val),
              -F.x * sin(F.val),
              -F.y * sin(F.val),
              -F.xy * sin(F.val) - F.x*F.y * cos(F.val));
}
var2 tan (const var2& F) { return sin(F)/cos(F); }
var2 sinh (const var2& F)
{ return var2(sinh(F.val),
              F.x * cosh(F.val),
              F.y * cosh(F.val),
              F.xy * cosh(F.val) + F.x*F.y * sinh(F.val));
}
var2 cosh (const var2& F)
{ return var2(cosh(F.val),
              F.x * sinh(F.val),
              F.y * sinh(F.val),
              F.xy * sinh(F.val) + F.x*F.y * cosh(F.val));
}
var2 tanh (const var2& F) { return sinh(F)/cosh(F); }



typedef var (*functionpointer)(var[]);
typedef var2 (*functionpointer2)(var2[]);
typedef var (*potentialpointer)(var,var,vector<ldf> *);



// Example functions

template<ui dim, class number> number function1 (number X[])
{	return pow(X[0],3) * exp(-X[1]) * cos(X[dim-1]);
}

template<ui dim, class number> number function2 (number X[])
{	return X[0]*X[0]/X[1] + X[1] + 7;
}

template<class number> number LJ (number r, number rsq, vector<ldf> *parameters)
{
    (void) r;
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return 4.0*e*(pow(s,12)/pow(rsq,6)-pow(s,6)/pow(rsq,3));;
}

template<class number> number LJ2 (number r, number rsq, vector<ldf> *parameters)
{
    (void) rsq;
    const ldf e=parameters->at(0);
    const ldf s=parameters->at(1);
    return 4.0*e*(pow(s,12)/pow(r,12)-pow(s,6)/pow(r,6));;
}



// Example usage

template<ui dim> ldf df (functionpointer func, ui i, ldf X[])
{	var Y[dim];
	for (ui d = 0; d < dim; d++)
		Y[d] = var(X[d], i==d);
	return func(Y).diff();
}

ldf dr (potentialpointer func, ldf r, ldf rsq, vector<ldf> *parameters)
{	var R(r,1), Rsq(rsq, 2*r);
	return func(R, Rsq, parameters).diff();
}

template<ui dim> ldf d2f (functionpointer2 func, ui i, ui j, ldf X[])
{	var2 Y[dim];
	for (ui d = 0; d < dim; d++)
		Y[d] = var2(X[d], i==d, j==d, 0);
	return func(Y).diff12();
}



int test()
{	ldf X[] = {1,2,3};
	ldf Y[] = {5,7};
	printf("%Lf %Lf  %Lf %Lf %Lf\n", df<3>(&function1<3>,0,X), df<3>(&function1<3>,2,X), d2f<3>(&function1<3>,0,0,X), d2f<3>(&function1<3>,0,1,X), d2f<3>(&function1<3>,1,1,X));
	printf("%Lf %Lf  %Lf %Lf %Lf\n", df<2>(&function2<2>,0,Y), df<2>(&function2<2>,1,Y), d2f<2>(&function2<2>,0,0,Y), d2f<2>(&function2<2>,0,1,Y), d2f<2>(&function2<2>,1,1,Y));
	vector<ldf> P(2);
	P[0] = 3.0;
	P[1] = 2.0;
	ldf r = 4.0;
	printf("%Lf %Lf\n", dr(&LJ, r, r*r, &P), dr(&LJ2, r, r*r, &P));
	return 0;
}

