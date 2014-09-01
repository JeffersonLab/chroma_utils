#include <cmath>
#include "nr.h"
using namespace std;

//namespace NR{

namespace {
	inline void shft3(NR::DP &a, NR::DP &b, NR::DP &c, const NR::DP d)
	{
		a=b;
		b=c;
		c=d;
	}
}



void NR::mnbrak(NR::DP &ax, NR::DP &bx, NR::DP &cx, NR::DP &fa, NR::DP &fb, NR::DP &fc,
	const ScalFunc &func)
//void NR::mnbrak(NR::DP &ax, NR::DP &bx, NR::DP &cx, NR::DP &fa, NR::DP &fb, NR::DP &fc,
//	NR::DP func(const NR::DP))
{
	const NR::DP GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
	NR::DP ulim,u,r,q,fu;

	fa=func(ax);
	fb=func(bx);
	if (fb > fa) {
		SWAP(ax,bx);
		SWAP(fb,fa);
	}
	cx=bx+GOLD*(bx-ax);
	fc=func(cx);
	while (fb > fc) {
		r=(bx-ax)*(fb-fc);
		q=(bx-cx)*(fb-fa);
		u=bx-((bx-cx)*q-(bx-ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=bx+GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu=func(u);
			if (fu < fc) {
				ax=bx;
				bx=u;
				fa=fb;
				fb=fu;
				return;
			} else if (fu > fb) {
				cx=u;
				fc=fu;
				return;
			}
			u=cx+GOLD*(cx-bx);
			fu=func(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu=func(u);
			if (fu < fc) {
				shft3(bx,cx,u,cx+GOLD*(cx-bx));
				shft3(fb,fc,fu,func(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u=ulim;
			fu=func(u);
		} else {
			u=cx+GOLD*(cx-bx);
			fu=func(u);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
}

void NR::mnbrak(NR::DP &ax, NR::DP &bx, NR::DP &cx, NR::DP &fa, NR::DP &fb, NR::DP &fc,
		NR::DP func(const NR::DP)){
  PtrScalFunc p(func);
  return NR::mnbrak(ax, bx, cx,fa, fb, fc, p);
}

//}
