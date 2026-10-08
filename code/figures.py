"""Render documentation figures from the report and saved benchmark records.

Outputs default to ignored results/latest/figures/. See assets/figures/README.md.
"""
from pathlib import Path
import argparse
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/latest/figures'
FORMATS=('svg','png')


def figure_data():
    rows=re.findall(r'^\| (\d+) \| (\d+) \| (\d+) \| ([\d.]+) \| ([^\n]+) \|$',
                    (ROOT/'REPORT.md').read_text(),re.M)
    frontier=[dict(budget=int(budget),degree=int(degree),N=int(N),rate=float(rate))
              for budget,degree,N,rate,status in rows if int(budget)<=32]
    if [row['budget'] for row in frontier]!=[2,4,8,16,32]:
        raise ValueError('Expected the five certified degree-budget rows in REPORT.md')
    return dict(frontier=frontier,
                benchmark=json.loads((ROOT/'results/benchmark_newton.json').read_text()))


DATA=figure_data()
BG='#F7F4EC'; INK='#25364A'; COPPER='#BC7046'; TEAL='#3E7D7C'
MUTED='#646B71'; GRID='#D8D8D1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,
    'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
    'axes.edgecolor':GRID,'axes.facecolor':BG,'figure.facecolor':BG,
    'savefig.facecolor':BG,'mathtext.fontset':'dejavuserif',
    'svg.fonttype':'path','svg.hashsalt':'ramanujan-series-figures-v1',
    'axes.spines.top':False,'axes.spines.right':False})


def base(title):
    fig=plt.figure(figsize=(16,9))
    fig.text(.055,.93,title,fontsize=27,fontfamily='DejaVu Serif')
    return fig


def save(fig,name):
    OUTPUT.mkdir(parents=True,exist_ok=True)
    for suffix in FORMATS:
        metadata={'Date':None,'Creator':'ramanujan-series/code/figures.py'} if suffix=='svg' else None
        path=OUTPUT/(name+'.'+suffix)
        fig.savefig(path,dpi=160,metadata=metadata)
        if suffix=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def torus_paths():
    fig=base('Infinite journeys. Two basic directions.')
    cases=[(1,0,INK,'Around the hole'),(0,1,COPPER,'Around the tube'),(2,3,TEAL,'Combining both directions')]
    uu,vv=np.meshgrid(np.linspace(0,2*np.pi,45),np.linspace(0,2*np.pi,25))
    R,r=2.05,.78
    xx=(R+r*np.cos(vv))*np.cos(uu); yy=(R+r*np.cos(vv))*np.sin(uu); zz=r*np.sin(vv)
    tt=np.linspace(0,2*np.pi,1400)
    for i,(m,n,color,label) in enumerate(cases):
        left=.04+i*.32
        fig.text(left+.155,.837,label,ha='center',fontsize=16,weight='bold',color=color)
        fig.text(left+.155,.791,f'Winding numbers ({m}, {n})',ha='center',fontsize=12,color=MUTED)
        ax=fig.add_axes([left,.49,.31,.29],projection='3d')
        ax.plot_wireframe(xx,yy,zz,color=INK,linewidth=.5,alpha=.14,rstride=2,cstride=3)
        u=m*tt;v=n*tt
        ax.plot((R+r*np.cos(v))*np.cos(u),(R+r*np.cos(v))*np.sin(u),r*np.sin(v),color=color,lw=3)
        ax.scatter([R+r],[0],[0],s=44,color=color,edgecolors=BG,linewidths=1.6,depthshade=False)
        ax.set_box_aspect([1,1,.5],zoom=1.4);ax.view_init(28,-58);ax.set_axis_off()
        ax.set(xlim=(-3,3),ylim=(-3,3),zlim=(-1,1))
        grid=fig.add_axes([left+.055,.16,.22,.285])
        grid.set_aspect('equal');grid.set_xlim(-.35,3.5);grid.set_ylim(-.35,3.5)
        for q in range(4):
            grid.plot([0,3],[q,q],color=GRID,lw=.8,zorder=0)
            grid.plot([q,q],[0,3],color=GRID,lw=.8,zorder=0)
        for a in range(4):
            for b in range(4):grid.plot(a,b,'o',color=GRID,ms=3)
        grid.plot([0,m],[0,n],color=color,lw=3,zorder=3)
        grid.annotate('',xy=(m*.67,n*.67),xytext=(m*.45,n*.45),
                      arrowprops=dict(arrowstyle='-|>',color=color,lw=2,mutation_scale=16))
        grid.plot(0,0,'o',color=color,ms=8,zorder=4)
        grid.plot(m,n,'o',mfc=BG,mec=color,mew=2.5,ms=10,zorder=4)
        grid.text(m+.16,n+.12,f'({m}, {n})',fontsize=12,color=color)
        grid.text(-.13,-.25,'(0, 0)',fontsize=10,ha='right',color=MUTED)
        grid.text(3.2,-.1,'m',fontsize=12,color=MUTED)
        grid.text(-.12,3.25,'n',fontsize=12,color=MUTED)
        grid.set_axis_off()
    fig.text(.5,.063,'Filled and open endpoints coincide when the lattice is wrapped into a torus.',ha='center',fontsize=14,color=MUTED)
    save(fig,'torus-and-lattice')


def frontier():
    fig=base('Convergence at five degree budgets')
    ax=fig.add_axes([.11,.20,.82,.60])
    points=DATA['frontier'];xs=np.arange(len(points));ys=[p['rate'] for p in points]
    ax.set_axisbelow(True);ax.grid(axis='y',color=GRID,lw=.9)
    ax.set_ylim(0,410);ax.set_xlim(-.42,4.5)
    ax.set_yticks([0,100,200,300,400]);ax.set_xticks(xs,[str(p['budget']) for p in points])
    ax.tick_params(axis='both',length=0,pad=11)
    ax.spines['left'].set_visible(False)
    ax.set_ylabel('Digits per term (asymptotic)',fontsize=15,labelpad=15)
    ax.set_xlabel('Allowed joint algebraic degree  ·  doubling scale',fontsize=15,labelpad=17)
    for x,y in zip(xs,ys):
        color=COPPER if x==4 else INK
        ax.plot([x,x],[0,y],color=color,alpha=.23,lw=2)
        ax.scatter([x],[y],s=125 if x==4 else 75,color=color,zorder=3)
        ax.annotate(f'{y:.2f}',(x,y),xytext=(0,15),textcoords='offset points',ha='center',fontsize=18,color=color,weight='bold')
    save(fig,'convergence-frontier')


def symmetry():
    fig=base('Two descriptions, one parameter')
    ax=fig.add_axes([.04,.14,.5,.67]);ax.set_xlim(0,10);ax.set_ylim(0,8);ax.axis('off')
    ax.text(1.65,7.6,r'$t_i$',fontsize=20,ha='center',color=INK)
    ax.text(4.35,7.6,r'$4096/t_i$',fontsize=20,ha='center',color=COPPER)
    ax.text(8.2,7.6,r'$x_i$',fontsize=20,ha='center',color=TEAL)
    for i,y in enumerate([6.65,4.95,3.25,1.55],1):
        ax.plot([1.65,4.35],[y,y],color=GRID,lw=1.5,zorder=0)
        for x,color,label in [(1.65,INK,str(i)),(4.35,COPPER,str(i))]:
            ax.add_patch(Circle((x,y),.35,facecolor=BG,edgecolor=color,lw=2.4))
            ax.text(x,y,label,ha='center',va='center',color=color,fontsize=13)
        for start,rad in [((2.05,y+.08),-.2),((4.75,y-.08),.12)]:
            ax.add_patch(FancyArrowPatch(start,(7.72,y),arrowstyle='-|>',connectionstyle=f'arc3,rad={rad}',
                                        mutation_scale=14,color=GRID,lw=1.5,zorder=0))
        ax.add_patch(Circle((8.2,y),.43,facecolor=TEAL,edgecolor=TEAL))
        ax.text(8.2,y,str(i),color='white',ha='center',va='center',fontsize=13)
    ax.text(3,.34,'8 modular values',ha='center',fontsize=15)
    ax.text(8.2,.34,'4 parameters',ha='center',fontsize=15,color=TEAL)
    fig.text(.6,.73,r'$x(t)=\frac{256t}{(t+64)^2}$',fontsize=25)
    fig.text(.6,.645,r'$x(4096/t)=x(t)$',fontsize=22,color=TEAL)
    fig.text(.6,.515,'The whole identity has to follow.',fontsize=17,weight='bold')
    fig.text(.6,.463,r'For $N\equiv5\;(\mathrm{mod}\,8)$:',fontsize=15)
    fig.text(.6,.404,r'$\sqrt{N}\ \mapsto\ -\sqrt{N}$',fontsize=19,color=COPPER)
    fig.text(.6,.347,r'$\sqrt{1-x}\ \mapsto\ -\sqrt{1-x}$',fontsize=19,color=COPPER)
    fig.text(.6,.283,r'$B=\sqrt{N}\sqrt{1-x}$ stays fixed.',fontsize=17,color=TEAL)
    fig.text(.6,.222,'A needs a separate CM descent argument.',fontsize=13,color=MUTED)
    save(fig,'symmetry-and-degree')


def benchmark():
    fig=base('Fewer terms is only part of the story.')
    handles=[Line2D([0],[0],lw=9,color=INK),Line2D([0],[0],lw=9,color=COPPER),
             Line2D([0],[0],marker='o',color='none',markerfacecolor=TEAL,markersize=7)]
    fig.legend(handles,['Coefficient setup','Summation + error bound','Median total'],
               loc='upper left',bbox_to_anchor=(.055,.862),ncol=3,frameon=False,fontsize=12,columnspacing=2.1)
    methods=['ramanujan-1103','chudnovsky','cm-degree-32']
    labels=['Ramanujan (1103)','Chudnovsky','CM degree 32']
    for j,(digits,scale,unit,limit,ticks) in enumerate([(1000,1000,'milliseconds',14,[0,3,6,9,12]),(10000,1,'seconds',.85,[0,.2,.4,.6,.8])]):
        ax=fig.add_axes([.19+j*.445,.215,.29,.46])
        rows={r['method']:r for r in DATA['benchmark']['rows'] if r['target_absolute_error_digits']==digits}
        ys=np.array([2,1,0]);ax.set_xlim(0,limit);ax.set_ylim(-.6,2.6)
        ax.set_yticks(ys,[labels[i]+'\n'+f"{rows[m]['terms']:,} terms" for i,m in enumerate(methods)])
        ax.set_xticks(ticks);ax.tick_params(axis='both',length=0,pad=8,labelsize=11)
        ax.set_axisbelow(True);ax.grid(axis='x',color=GRID,lw=.8);ax.spines['left'].set_visible(False)
        ax.set_xlabel('Time ('+unit+')  ·  lower is faster',fontsize=12,labelpad=12)
        ax.set_title(f'{digits:,}-digit error target',loc='left',fontsize=18,pad=30,fontfamily='DejaVu Serif')
        for y,m in zip(ys,methods):
            r=rows[m];setup=r['median_setup_seconds']*scale;total=r['median_total_seconds']*scale;add=r['median_summation_and_bound_seconds']*scale
            ax.barh(y,setup,height=.36,color=INK,zorder=3)
            ax.barh(y,add,left=setup,height=.36,color=COPPER,zorder=3)
            ax.plot(total,y,'o',color=TEAL,ms=6,zorder=4)
            text=f'{total:.2f} ms' if j==0 else f'{total:.3f} s'
            ax.text(max(total,setup+add)+limit*.035,y,text,va='center',fontsize=12,weight='bold')
    save(fig,'terms-and-time')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    parser.add_argument('--format',nargs='+',choices=('svg','png'),default=FORMATS)
    args=parser.parse_args()
    OUTPUT=args.output_dir
    FORMATS=args.format
    torus_paths();frontier();symmetry();benchmark()
    print('Rendered four figures in:',OUTPUT)
