"""Exact certificates for large rectangular holes in a square.

Companion to Section 6 and Appendix D of manuscript.tex, entitled
"Polya inequalities for displaced square holes and perforated squares".
Run with Python 3.8+; only the standard library is required.
All proof decisions are rational or integer. No floating point is used.
The optional --output PATH saves a deterministic result record; no file
is read, and default invocation only prints the record.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json

P=F(333,106)
H=1<<24


def require(c,s):
    if not c: raise RuntimeError(s)


def pi_certificate():
    def atan_bounds(d):
        s=sum((1 if k%2==0 else -1)*F(1,(2*k+1)*d**(2*k+1)) for k in range(12))
        return s,s+F(1,25*d**25)
    a,b=atan_bounds(5)
    c,d=atan_bounds(239)
    require(P<16*a-4*d<16*b-4*c<F(355,113),'pi enclosure')
    return dict(method='Machin identity with twelve-term rational alternating enclosures',
                lower='333/106',upper='355/113',passed=True)


def up(z):
    require(z>=0,'negative radical')
    return isqrt(z.numerator*H*H//z.denominator)+1


def down(z):
    return up(z)-1


@lru_cache(maxsize=8192)
def midpoint(r,a):
    if a==0:return 0
    out=0
    k=0
    while True:
        z=r-F((2*k+1)**2,4)/(a*a)
        if z<0:return out
        out+=up(z)
        k+=1


@lru_cache(maxsize=8192)
def riesz(z):
    return sum(up(z*z-k*k) for k in range(1,z.numerator//z.denominator+1))


@lru_cache(maxsize=8192)
def riesz_normalized(z):
    return sum(up(1-F(k*k)/(z*z)) for k in range(1,z.numerator//z.denominator+1))


def midpoint_refined():
    stack=[(F(1,4),F(9),0)]
    tested=accepted=depthmax=0
    width=F(0)
    least=None
    while stack:
        l,u,d=stack.pop()
        tested+=1
        depthmax=max(depthmax,d)
        margin=P*l/4+F(4,25)*F(down(l),H)-F(midpoint(u,F(1)),H)
        if margin>=0:
            accepted+=1;width+=u-l
            least=margin if least is None else min(least,margin)
        else:
            require(d<30,f'midpoint failed {l} {u}')
            m=(l+u)/2
            stack.extend([(l,m,d+1),(m,u,d+1)])
    require(width==F(35,4),'midpoint coverage')
    require(F(2,27)/3<=F(4,25)**2,'midpoint tail')
    require((tested,accepted,depthmax)==(245,123,12),'midpoint statistics')
    return dict(bound='4/25',y_range='[1/2,3]',tail_start=3,tested=tested,
                accepted=accepted,depth=depthmax,covered_squared_width=str(width),least_margin=str(least))


def L_certificate(start,stop,target):
    tested=accepted=depthmax=0
    covered=F(0)
    least=None
    for n in range(start,stop):
        stack=[(F(n),F(n+1),0)]
        while stack:
            l,u,d=stack.pop()
            tested+=1;depthmax=max(depthmax,d)
            margin=11*P*l*l/4-11*F(riesz(u),H)-n-12*target*u
            if margin>=0:
                accepted+=1;covered+=u-l
                least=margin if least is None else min(least,margin)
            else:
                require(d<30,f'L failed {l} {u} {target}')
                m=(l+u)/2
                stack.extend([(l,m,d+1),(m,u,d+1)])
    require(covered==stop-start,'L interval coverage')
    # L(z) >= 3/8 - (11/12)*sqrt(2/27)/sqrt(z).
    require(F(3,8)>target,'L tail sign')
    require(stop*(F(3,8)-target)**2 >= F(11,12)**2*F(2,27),'L tail')
    if (start,stop,target)==(1,7,F(7,25)):
        require((tested,accepted,depthmax)==(434,220,9),'global L statistics')
    elif (start,stop,target)==(4,21,F(8,25)):
        require((tested,accepted,depthmax)==(1995,1006,10),'high L statistics')
    else:
        raise RuntimeError('Unexpected L certificate parameters')
    return dict(z_range=f'[{start},{stop})',lower_bound=str(target),tail_start=stop,
                tested=tested,accepted=accepted,depth=depthmax,covered_width=str(covered),least_margin=str(least))


def horizontal_middle():
    # H(x,a) <= 4/15 for 2<=x<=4, 0<=a<=1/2.
    G=K=1<<14
    stack=[(4*G,16*G,0,K//2,0)]
    tested=accepted=depthmax=covered=0
    least=None
    while stack:
        rl,ru,al,au,d=stack.pop()
        tested+=1;depthmax=max(depthmax,d)
        l,u,a,b=F(rl,G),F(ru,G),F(al,K),F(au,K)
        margin=P*l/4+F(4,15)*F(down(l),H)-F(midpoint(u,b)+midpoint(u,1-a),H)
        if margin>=0:
            accepted+=1;covered+=(ru-rl)*(au-al)
            least=margin if least is None else min(least,margin)
            continue
        require(ru-rl>1 or au-al>1,f'horizontal failed {l} {u} {a} {b}')
        splitr=(ru-rl)*K>=2*(au-al)*ru
        if (splitr and ru-rl>1) or au-al<=1:
            m=(rl+ru)//2
            stack.extend([(rl,m,al,au,d+1),(m,ru,al,au,d+1)])
        else:
            m=(al+au)//2
            stack.extend([(rl,ru,al,m,d+1),(rl,ru,m,au,d+1)])
    require(covered==12*G*(K//2),'horizontal coverage')
    require((tested,accepted,depthmax)==(5395,2698,18),'horizontal statistics')
    return dict(x_range='[2,4]',a_range='[0,1/2]',upper_bound='4/15',
                tested=tested,accepted=accepted,depth=depthmax,covered_grid_area=covered,
                expected_grid_area=12*G*(K//2),least_margin=str(least))


def triangle_double_area(xl,xu,zl,zu):
    left=max(0,min(xu,zl)-xl)*(zu-zl)*2
    lo=max(xl,zl)
    hi=min(xu,zu)
    if hi>lo:left+=(hi-lo)*(2*zu-hi-lo)
    return left


def center_triangle():
    # h_center(x) <= K(z,x) for 1<=x<=2, x<=z<4.
    # floor(z)=n on each half-open initial z interval.
    G=1<<16
    tested=accepted=skipped=depthmax=covered=0
    least=None
    per_band=[]
    for n in range(1,4):
        stack=[(G,2*G,n*G,(n+1)*G,0)]
        area=0
        while stack:
            xl,xu,zl,zu,d=stack.pop()
            clipped=triangle_double_area(xl,xu,zl,zu)
            if clipped==0:
                skipped+=1
                continue
            tested+=1;depthmax=max(depthmax,d)
            a,b,c,e=F(xl,G),F(xu,G),F(zl,G),F(zu,G)
            lower=P*(c/4+11*a/48)
            upper=(1-a/(12*e))*F(riesz_normalized(e),H)+F(n)/(12*c)+2*F(up(1-1/(b*b)),H)
            margin=lower-upper
            if margin>=0:
                accepted+=1;area+=clipped;covered+=clipped
                least=margin if least is None else min(least,margin)
                continue
            require(xu-xl>1 or zu-zl>1,f'center triangle failed {a} {b} {c} {e}')
            if (xu-xl>=zu-zl and xu-xl>1) or zu-zl<=1:
                m=(xl+xu)//2
                stack.extend([(xl,m,zl,zu,d+1),(m,xu,zl,zu,d+1)])
            else:
                m=(zl+zu)//2
                stack.extend([(xl,xu,zl,m,d+1),(xl,xu,m,zu,d+1)])
        expected=G*G if n==1 else 2*G*G
        require(area==expected,f'triangle band {n} coverage')
        per_band.append(dict(floor_z=n,twice_clipped_grid_area=area,expected=expected))
    require(covered==5*G*G,'total triangle coverage')
    require((tested,accepted,skipped,depthmax)==(9280,4620,43,18),'triangle statistics')
    return dict(x_range='[1,2]',z_range='[x,4)',grid_denominator=G,
                tested=tested,accepted=accepted,skipped_empty_boxes=skipped,depth=depthmax,
                twice_clipped_grid_area=covered,coverage_by_band=per_band,least_margin=str(least))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result={'pi_certificate':pi_certificate()}
    for name,fn in [('midpoint_refined',midpoint_refined),
                    ('L_global',lambda:L_certificate(1,7,F(7,25))),
                    ('L_high',lambda:L_certificate(4,21,F(8,25))),
                    ('horizontal_middle',horizontal_middle),
                    ('center_triangle',center_triangle)]:
        result[name]=fn()
    result['all_checks_passed']=True
    print(json.dumps(result,indent=2),flush=True)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
