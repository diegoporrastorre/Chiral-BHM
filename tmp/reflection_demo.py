import numpy as np
from scipy.sparse import diags, bmat, eye
from scipy.sparse.linalg import expm_multiply
import sys
sys.path.append('C:/Users/diego/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/Lib/site-packages')
from PIL import Image, ImageDraw, ImageFont

N=24; kappa=3.; delta=.1; t=np.linspace(0,60,1201)
def evolve(n, initial, delta):
    s=diags([np.ones(n-1)],[1],shape=(n,n),format='csr'); ident=eye(n,format='csr')
    a=bmat([[(1-kappa/2)*ident+2*s,delta*ident],[-delta*ident,-(1+kappa/2)*ident-2*s.T]],format='csr')
    y=np.zeros(2*n); y[initial]=1
    return expm_multiply(a,y,start=0,stop=60,num=len(t),traceA=-kappa*n)
y=evolve(N,17,delta); z=evolve(N,17,0)
pad=180; big=evolve(N+2*pad,pad+17,delta)
bulk=np.concatenate([big[:,pad:pad+N],big[:,N+2*pad+pad:N+2*pad+pad+N]],axis=1)
diff=y-bulk
# Independent time-step refinement by fourth-order Runge-Kutta.
def rhs(v):
    p,q=v[:N],v[N:]; dp=-.5*p+.1*q; dq=-2.5*q-.1*p
    dp[:-1]+=2*p[1:]; dq[1:]-=2*q[:-1]
    return np.r_[dp,dq]
v=np.zeros(2*N);v[17]=1;dt=.005
for _ in range(12000):
    a=rhs(v);b=rhs(v+dt*a/2);c=rhs(v+dt*b/2);d=rhs(v+dt*c);v+=dt*(a+2*b+2*c+d)/6
print('RK4 relative final error',np.linalg.norm(v-y[-1])/np.linalg.norm(y[-1]))
norm=np.linalg.norm(y,axis=1); norm0=np.linalg.norm(z,axis=1)
print('late fitted amplitude growth',np.polyfit(t[t>=45],np.log(norm[t>=45]),1)[0])
for time in [0,5,10,15,20,30,40,60]:
    i=round(time/.05);print(time,norm[i],norm0[i],np.linalg.norm(diff[i]))

im=Image.new('RGB',(1560,1110),'#fbfcfe');draw=ImageDraw.Draw(im)
fontpath='C:/Windows/Fonts/arial.ttf'
def f(n):return ImageFont.truetype(fontpath,n)
def txt(x,y,s,size=22,fill='#172438'):draw.text((x,y),s,font=f(size),fill=fill)
txt(55,25,'A localized kick in a detuned, open Bogoliubov chain',32)
txt(55,72,'N = 24   |   J = gs = gc = 1   |   phi = pi/2   |   Delta = 0.1   |   kappa = 3   |   p18(0) = 1, q(0) = 0',21)
colors=np.array([[17,25,54],[37,71,135],[36,148,166],[244,197,74],[254,247,215]])
def color(val):
    u=np.clip((val+5)/23,0,1)*4;k=np.minimum(u.astype(int),3);a=(u-k)[...,None]
    return ((1-a)*colors[k]+a*colors[k+1]).astype('uint8')
def heat(x,y,data,title):
    w,h=640,320
    arr=color(np.log10(np.maximum(abs(data[::2]),1e-30)))
    panel=Image.fromarray(arr).resize((w,h),Image.Resampling.NEAREST)
    im.paste(panel,(x,y));txt(x,y-34,title,23)
    for time in [0,10,20,30,40,50,60]:
        yy=y+time/60*h;draw.line((x-5,yy,x,yy),fill='black');txt(x-38,yy-9,str(time),17)
    for site in [1,6,12,18,24]:txt(x+(site-.5)/24*w-8,y+h+7,str(site),17)
    txt(x+270,y+h+30,'site j',18)
heat(85,160,y[:,:N],'Left-moving channel |p|')
heat(855,160,y[:,N:],'Right-moving channel |q|')
txt(15,270,'t',23)
heat(85,590,diff[:,N:],'Boundary-induced return: |q open - q long|')
txt(15,700,'t',23)
# Time histories on log scale.
x,y0,w,h=855,590,640,320
txt(x,y0-34,'Total perturbation amplitude',23)
for power in [-3,0,3,6,9,12,15,18]:
    yy=y0+h-(power+3)/21*h;draw.line((x,yy,x+w,yy),fill='#d8dee8');txt(x-52,yy-9,'10^'+str(power),15)
for time in [0,10,20,30,40,50,60]:
    xx=x+time/60*w;draw.line((xx,y0,xx,y0+h),fill='#e3e7ee');txt(xx-8,y0+h+7,str(time),17)
for values,col in [(norm,'#c34829'),(norm0,'#275db2')]:
    pts=[(x+tt/60*w,y0+h-(np.clip(np.log10(max(val,1e-30)),-3,18)+3)/21*h) for tt,val in zip(t,values)]
    draw.line(pts,fill=col,width=3)
txt(x+20,y0+15,'Delta = 0.1: feedback and eventual growth',18,'#c34829')
txt(x+20,y0+43,'Delta = 0: large transient, then decay',18,'#275db2')
txt(x+270,y0+h+30,'time t J',18)
# Shared color scale for the three maps.
bar=np.linspace(-5,18,500)[None,:];im.paste(Image.fromarray(color(bar)).resize((500,20)),(85,997))
for value in [-5,0,5,10,15,18]:txt(85+(value+5)/23*500-15,1022,'10^'+str(value),17)
txt(620,994,'Amplitude scale (same in all three maps)',19)
txt(85,1060,'Unit initial amplitude is a normalization; multiply every field by an arbitrarily small initial kick.',18)
im.save('tmp/reflection_demo.png')

