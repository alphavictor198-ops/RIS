#!/usr/bin/env python3
"""RIS-Assisted Satellite Link Simulator — Full Parameter Edition"""
import sys,subprocess
def _e(p):
    try:__import__(p.replace("-","_"))
    except ImportError:subprocess.check_call([sys.executable,"-m","pip","install",p,"-q"])
for _p in ["PyQt6","numpy","matplotlib"]:_e(_p)
import numpy as np,matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FC
from matplotlib.figure import Figure
import matplotlib.patches as mp
from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt,QTimer
C=3e8
# ─── Theme ───
BG="#f4f4f6";CD="#ffffff";PN="#eaeaee";LN="#d0d0d8"
CY="#0891b2";TX="#1e1e2e";DM="#6b7280"
RD="#dc2626";GN="#16a34a";BL="#2563eb";AM="#d97706"
CSS=f"""
QMainWindow,QWidget{{background:{BG};color:{TX};font-family:'Segoe UI',sans-serif;font-size:13px;}}
QPushButton{{background:{CY};color:#fff;border:none;padding:9px 0;border-radius:6px;font-weight:700;font-size:13px;}}
QPushButton:hover{{background:#0e7490;}}
QLabel{{color:{TX};background:transparent;}}
QDoubleSpinBox,QSpinBox{{background:#fff;color:{BL};border:1px solid {LN};border-radius:5px;padding:3px 6px;font-size:12px;font-weight:700;min-width:75px;}}
QDoubleSpinBox:focus,QSpinBox:focus{{border-color:{CY};}}
QDoubleSpinBox::up-button,QDoubleSpinBox::down-button,QSpinBox::up-button,QSpinBox::down-button{{background:{PN};border:none;width:16px;border-radius:2px;}}
QDoubleSpinBox::up-arrow,QSpinBox::up-arrow{{image:none;width:0;height:0;border-left:3px solid transparent;border-right:3px solid transparent;border-bottom:4px solid {CY};margin:auto;}}
QDoubleSpinBox::down-arrow,QSpinBox::down-arrow{{image:none;width:0;height:0;border-left:3px solid transparent;border-right:3px solid transparent;border-top:4px solid {CY};margin:auto;}}
QComboBox{{background:#fff;color:{BL};border:1px solid {LN};border-radius:5px;padding:3px 6px;font-size:12px;font-weight:700;min-width:75px;}}
QComboBox QAbstractItemView{{background:#fff;color:{BL};selection-background-color:{CY};selection-color:#fff;}}
QScrollArea{{border:none;background:{PN};}}
"""
# ─── Physics ───
def slant_range(alt_km,el):
    R=6371e3;h=alt_km*1e3;t=np.deg2rad(el)
    return -R*np.sin(t)+np.sqrt((R*np.sin(t))**2+2*R*h+h*h)

def rain_atten(fg,rr,el):
    _f=[1,2,4,8,10,15,20,25,30,50,100]
    _k=[2.6e-5,1.54e-4,6.5e-4,4.54e-3,1.01e-2,3.67e-2,7.51e-2,.124,.187,.536,1.31]
    _a=[.97,1.07,1.12,1.33,1.26,1.15,1.10,1.06,1.02,.86,.74]
    k=np.interp(fg,_f,_k);a=np.interp(fg,_f,_a)
    return k*rr**a*3.0/np.sin(np.deg2rad(max(el,5)))

def pq_loss(bits):
    if bits==0:return 1.0
    return np.sinc(1.0/2**bits)**2

def compute(ptx,fghz,bw_mhz,alt,el,gtx,grx,nf,nris,espac,eta,pbits,kric,rr,perr,rdist):
    f=fghz*1e9;lam=C/f;d=slant_range(alt,el)
    Lfs=20*np.log10(4*np.pi*d*f/C)
    Latm=0.5/np.sin(np.deg2rad(max(el,1)))
    Lrain=rain_atten(fghz,rr,el) if rr>0 else 0
    Nfl=-174+10*np.log10(bw_mhz*1e6)+nf
    pno=ptx+gtx+grx-Lfs-Latm-Lrain
    de=espac*lam;A=de**2;eq=pq_loss(pbits)
    Ns=max(int(np.sqrt(nris)),1)
    bwdeg=np.rad2deg(lam/(Ns*de))
    Lpt=min(12*(perr/max(bwdeg,.01))**2,30)
    Kl=10**(kric/10);rc=Kl/(Kl+1)
    pris=(ptx+gtx+grx+20*np.log10(nris*A)+20*np.log10(eta)
          +10*np.log10(eq)+10*np.log10(rc)-Lpt
          -10*np.log10(64*np.pi**3)-20*np.log10(d)-20*np.log10(rdist)
          -Latm-Lrain)
    pd=10**(pno/10);pr=10**(pris/10)
    pw=10*np.log10((np.sqrt(pd)+np.sqrt(pr))**2)
    mu=3.986e14;vorb=np.sqrt(mu/((6371+alt)*1e3))
    dop=f*vorb*np.cos(np.deg2rad(el))/C
    return dict(lam=lam,d=d,Lfs=Lfs,Latm=Latm,Lrain=Lrain,Lpt=Lpt,
                eqdb=10*np.log10(eq) if eq>0 else -30,Nfl=Nfl,
                pno=pno,pw=pw,gain=pw-pno,snr_no=pno-Nfl,snr_w=pw-Nfl,
                dop=dop,vorb=vorb,bwdeg=bwdeg)
# ─── UI helpers ───
def _sep():
    f=QFrame();f.setFrameShape(QFrame.Shape.HLine)
    f.setStyleSheet(f"background:{LN};max-height:1px;border:none;");return f
def _lb(t,s=12,c=TX,b=False):
    w=QLabel(t);w.setStyleSheet(f"color:{c};font-size:{s}px;font-weight:{'700' if b else '400'};background:transparent;");return w

class FS(QWidget):
    def __init__(s,lab,lo,hi,v,st=1.,dc=1,cb=None):
        super().__init__();s.setFixedHeight(26)
        r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch()
        s.b=QDoubleSpinBox();s.b.setRange(lo,hi);s.b.setValue(v)
        s.b.setSingleStep(st);s.b.setDecimals(dc)
        s.b.setAlignment(Qt.AlignmentFlag.AlignRight);r.addWidget(s.b)
        if cb:s.b.valueChanged.connect(cb)
    def value(s):return s.b.value()

class IS(QWidget):
    def __init__(s,lab,lo,hi,v,cb=None):
        super().__init__();s.setFixedHeight(26)
        r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch()
        s.b=QSpinBox();s.b.setRange(lo,hi);s.b.setValue(v)
        s.b.setAlignment(Qt.AlignmentFlag.AlignRight);r.addWidget(s.b)
        if cb:s.b.valueChanged.connect(cb)
    def value(s):return s.b.value()

class CB(QWidget):
    def __init__(s,lab,items,idx=0,cb=None):
        super().__init__();s.setFixedHeight(26)
        r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch()
        s.b=QComboBox();s.b.addItems(items);s.b.setCurrentIndex(idx);s.b.setFixedWidth(85)
        r.addWidget(s.b)
        if cb:s.b.currentIndexChanged.connect(cb)
    def text(s):return s.b.currentText()

class Card(QWidget):
    def __init__(s,lab,col=CY):
        super().__init__();s.setFixedHeight(44)
        s.setStyleSheet(f"background:{CD};border-radius:5px;border:1px solid {LN};border-left:3px solid {col};")
        l=QVBoxLayout(s);l.setContentsMargins(8,2,8,2);l.setSpacing(0)
        s.v=QLabel("--");s.v.setStyleSheet(f"color:{col};font-size:13px;font-weight:700;background:transparent;border:none;")
        s.l=QLabel(lab);s.l.setStyleSheet(f"color:{DM};font-size:8px;background:transparent;border:none;")
        l.addWidget(s.v);l.addWidget(s.l)
    def set(s,t):s.v.setText(str(t))

# ─── Drawing ───
def _sat(ax,cx,cy):
    ax.add_patch(mp.FancyBboxPatch((cx-.2,cy-.1),.4,.2,boxstyle="round,pad=0.02",fc="#2563eb",ec="#3b82f6",lw=1.4,zorder=10))
    for dx in [-.58,.2]:
        ax.add_patch(plt.Rectangle((cx+dx,cy-.05),.34,.1,fc="#1d4ed8",ec="#3b82f6",lw=.8,zorder=10))
    ax.plot([cx,cx],[cy+.1,cy+.22],color="#60a5fa",lw=1.2,zorder=11)
    ax.plot(cx,cy+.24,'o',color="#60a5fa",ms=2.5,zorder=11)
    ax.text(cx,cy+.38,"CubeSat",ha="center",fontsize=9,color="#1d4ed8",fontweight="bold")

def _gs(ax,cx,cy):
    ax.plot([cx,cx],[cy,cy+.4],color="#475569",lw=2,zorder=10)
    ax.add_patch(plt.Polygon([(cx-.25,cy+.65),(cx+.03,cy+.33),(cx+.03,cy+.65)],closed=True,fc="#94a3b8",ec="#475569",lw=1.2,zorder=10))
    ax.plot([cx-.14,cx+.14],[cy,cy],color="#475569",lw=2.5,zorder=10)
    ax.text(cx,cy-.2,"Ground Station",ha="center",fontsize=8,color="#475569",fontweight="bold")

def _ris(ax,cx,cy,n):
    co=min(int(np.sqrt(n)),8);ro=min(n//max(co,1),8)
    w=co*.07;h=ro*.07
    ax.add_patch(mp.FancyBboxPatch((cx-w/2-.04,cy-h/2-.04),w+.08,h+.08,boxstyle="round,pad=0.02",fc="#ecfeff",ec=CY,lw=1.2,zorder=10))
    for r in range(ro):
        for c in range(co):
            ax.add_patch(plt.Rectangle((cx-w/2+c*.07,cy-h/2+r*.07),.05,.05,fc="#06b6d4",ec=CY,lw=.3,alpha=.7,zorder=11))
    ax.plot([cx,cx],[cy-h/2-.04,cy-h/2-.25],color="#64748b",lw=1.5,zorder=9)
    ax.text(cx,cy+h/2+.18,f"RIS ({n} elem.)",ha="center",fontsize=8,color=CY,fontweight="bold")

def draw_dia(ax,r,n,el,alt):
    ax.clear();ax.set_xlim(-.5,10.5);ax.set_ylim(-.2,6.8)
    ax.set_aspect("equal");ax.axis("off");ax.set_facecolor("#ffffff")
    sx,sy=5,5.8;gx,gy=7.5,.9;rx,ry=2.5,1.2
    ax.axhline(.4,color="#d1d5db",lw=1);ax.fill_between([-.5,10.5],-.2,.4,color="#f9fafb")
    _sat(ax,sx,sy);_gs(ax,gx,gy);_ris(ax,rx,ry,n)
    ax.annotate("",xy=(gx-.1,gy+.7),xytext=(sx+.1,sy-.15),arrowprops=dict(arrowstyle="-|>",color=RD,lw=1.5,ls="--",connectionstyle="arc3,rad=.03"),zorder=8)
    ax.text((sx+gx)/2+.35,(sy+gy)/2+.2,"Direct Path",fontsize=8,color=RD,fontweight="bold",ha="center",bbox=dict(fc="white",ec="none",alpha=.85,pad=1.5))
    ax.annotate("",xy=(rx+.1,ry+.35),xytext=(sx-.15,sy-.15),arrowprops=dict(arrowstyle="-|>",color=CY,lw=1.3,connectionstyle="arc3,rad=.02"),zorder=8)
    ax.text((sx+rx)/2-.4,(sy+ry)/2+.2,"h_sr",fontsize=8,color=CY,fontweight="bold",ha="center",bbox=dict(fc="white",ec="none",alpha=.85,pad=1.5))
    ax.annotate("",xy=(gx-.3,gy+.35),xytext=(rx+.4,ry-.05),arrowprops=dict(arrowstyle="-|>",color=BL,lw=1.3,connectionstyle="arc3,rad=.02"),zorder=8)
    ax.text((rx+gx)/2,ry-.3,"h_rg",fontsize=8,color=BL,fontweight="bold",ha="center",bbox=dict(fc="white",ec="none",alpha=.85,pad=1.5))
    ax.text(gx+.35,gy+.2,f"{el:.0f}°",fontsize=9,color=AM,fontweight="bold")
    ax.text(.3,3.3,f"Alt: {alt:.0f} km",fontsize=8,color=DM,rotation=90,ha="center")
    if r:
        parts=[f"FSPL: {r['Lfs']:.1f} dB",f"Atm: {r['Latm']:.1f} dB"]
        if r['Lrain']>0.01:parts.append(f"Rain: {r['Lrain']:.1f} dB")
        if abs(r['Lpt'])>0.01:parts.append(f"Point: {r['Lpt']:.1f} dB")
        if abs(r['eqdb'])>0.01:parts.append(f"Quant: {r['eqdb']:.1f} dB")
        parts.append(f"λ={r['lam']*100:.2f} cm")
        parts.append(f"Doppler: {r['dop']/1e3:.1f} kHz")
        ax.text(5.2,-.05,"  |  ".join(parts),fontsize=7,color=DM,ha="center",family="monospace")
    ax.set_title("System Architecture  —  RIS-Assisted Satellite Link",fontsize=11,fontweight="bold",color=TX,pad=6)

def draw_bars(ax,r):
    ax.clear();ax.set_facecolor("#ffffff")
    for s in ax.spines.values():s.set_color("#e5e5ea")
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.tick_params(colors=TX,labelsize=10)
    if r is None:
        ax.text(.5,.5,"Adjust parameters — auto-updates",transform=ax.transAxes,ha="center",va="center",color=DM,fontsize=13,style="italic");return
    pno=r["pno"];pw=r["pw"]
    xmin=min(pno,pw)-10;xmax=max(pno,pw)+35
    y=[1,0];vals=[pno,pw]
    ax.barh(y,[v-xmin for v in vals],left=xmin,height=.4,color=[RD,BL],alpha=.82,zorder=5)
    ax.text(pno+1,1,f"  {pno:.1f} dBm  (SNR {r['snr_no']:.1f} dB)",va="center",ha="left",fontsize=10,fontweight="bold",color=RD,zorder=6)
    ax.text(pw+1,0,f"  {pw:.1f} dBm  (SNR {r['snr_w']:.1f} dB)",va="center",ha="left",fontsize=10,fontweight="bold",color=BL,zorder=6)
    mx=max(pno,pw)
    ax.text(mx+1,.5,f"  ▲ +{r['gain']:.1f} dB gain",va="center",ha="left",fontsize=11,fontweight="bold",color=GN,zorder=6,bbox=dict(fc="#f0fdf4",ec=GN,alpha=.9,pad=3,boxstyle="round,pad=0.3"))
    ax.set_yticks(y);ax.set_yticklabels(["Without RIS","With RIS"],fontsize=10,fontweight="600")
    ax.set_xlim(xmin,xmax);ax.set_xlabel("Received Power (dBm)",fontsize=10,color=DM)
    ax.grid(True,axis="x",color="#e5e5ea",lw=.7,ls="--");ax.set_axisbelow(True)
    ax.set_title("Received Power Comparison",fontsize=11,fontweight="bold",color=TX,pad=6)

# ─── Main Window ───
class Win(QMainWindow):
    def __init__(s):
        super().__init__()
        s.setWindowTitle("RIS Satellite Link Simulator");s.setMinimumSize(1150,720);s.setStyleSheet(CSS)
        root=QWidget();s.setCentralWidget(root)
        ml=QHBoxLayout(root);ml.setSpacing(0);ml.setContentsMargins(0,0,0,0)
        # Side panel with scroll
        scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setFixedWidth(280)
        scroll.setStyleSheet(f"background:{PN};border-right:1px solid {LN};")
        side=QWidget();sl=QVBoxLayout(side);sl.setContentsMargins(12,10,12,10);sl.setSpacing(4)
        up=lambda _:s._au()
        sl.addWidget(_lb("Satellite Link Simulator",13,CY,True));sl.addWidget(_sep())
        # Link
        sl.addWidget(_lb("▸ Link Parameters",10,DM,True))
        s.p_freq=FS("Carrier Freq (GHz)",.1,200,2,.5,1,up)
        s.p_bw=FS("Bandwidth (MHz)",.1,5000,10,1,1,up)
        s.p_ptx=FS("Tx Power (dBm)",0,60,33,1,1,up)
        for w in[s.p_freq,s.p_bw,s.p_ptx]:sl.addWidget(w)
        sl.addWidget(_sep())
        # Antennas
        sl.addWidget(_lb("▸ Antennas",10,DM,True))
        s.p_gtx=FS("Tx Gain (dBi)",0,40,6,1,1,up)
        s.p_grx=FS("Rx Gain (dBi)",0,40,10,1,1,up)
        s.p_nf=FS("Noise Figure (dB)",0,15,3,.5,1,up)
        for w in[s.p_gtx,s.p_grx,s.p_nf]:sl.addWidget(w)
        sl.addWidget(_sep())
        # Satellite
        sl.addWidget(_lb("▸ Satellite Orbit",10,DM,True))
        s.p_alt=FS("Altitude (km)",100,2000,550,50,0,up)
        s.p_el=FS("Elevation Angle (°)",5,90,10,1,1,up)
        for w in[s.p_alt,s.p_el]:sl.addWidget(w)
        sl.addWidget(_sep())
        # RIS
        sl.addWidget(_lb("▸ RIS Configuration",10,DM,True))
        s.p_nris=IS("Elements (N)",4,4096,100,up)
        s.p_esp=FS("Spacing (λ)",.1,1,.5,.1,1,up)
        s.p_eta=FS("Reflect. Efficiency (η)",.1,1,.9,.05,2,up)
        s.p_pb=CB("Phase Bits",["∞ (cont.)","2-bit","4-bit"],0,up)
        s.p_rdist=FS("RIS→GS Dist (m)",1,500,50,5,0,up)
        for w in[s.p_nris,s.p_esp,s.p_eta,s.p_pb,s.p_rdist]:sl.addWidget(w)
        sl.addWidget(_sep())
        # Channel / Environment
        sl.addWidget(_lb("▸ Channel & Environment",10,DM,True))
        s.p_kric=FS("Rician K-factor (dB)",0,30,15,1,1,up)
        s.p_rain=FS("Rain Rate (mm/h)",0,150,0,5,0,up)
        s.p_perr=FS("Pointing Error (°)",0,10,0,.1,1,up)
        for w in[s.p_kric,s.p_rain,s.p_perr]:sl.addWidget(w)
        sl.addWidget(_sep())
        # Result cards
        s.c_no=Card("Rx Power — No RIS",RD)
        s.c_w=Card("Rx Power — With RIS",BL)
        s.c_g=Card("Power Gain",GN)
        s.c_snr=Card("SNR Improvement",AM)
        s.c_nf=Card("Noise Floor",DM)
        s.c_dop=Card("Doppler Shift","#7c3aed")
        for c in[s.c_no,s.c_w,s.c_g,s.c_snr,s.c_nf,s.c_dop]:sl.addWidget(c)
        sl.addStretch()
        scroll.setWidget(side);ml.addWidget(scroll)
        # Right: diagram + bars
        right=QWidget();right.setStyleSheet(f"background:{CD};")
        rl=QVBoxLayout(right);rl.setContentsMargins(8,6,8,6);rl.setSpacing(4)
        s.fdia=Figure(figsize=(10,4),facecolor="#fff");s.axd=s.fdia.add_subplot(111)
        s.cvd=FC(s.fdia);s.cvd.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Expanding)
        rl.addWidget(s.cvd,stretch=58)
        s.fbar=Figure(figsize=(10,2),facecolor="#fff");s.axb=s.fbar.add_subplot(111)
        s.cvb=FC(s.fbar);s.cvb.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Expanding)
        rl.addWidget(s.cvb,stretch=42)
        ml.addWidget(right)
        s._tm=QTimer();s._tm.setSingleShot(True);s._tm.setInterval(120);s._tm.timeout.connect(s._run)
        QTimer.singleShot(50,s._run)

    def _au(s):s._tm.start()

    def _pb(s):
        t=s.p_pb.text()
        if "∞" in t:return 0
        return int(t.split("-")[0])

    def _run(s):
        r=compute(s.p_ptx.value(),s.p_freq.value(),s.p_bw.value(),
                  s.p_alt.value(),s.p_el.value(),s.p_gtx.value(),s.p_grx.value(),
                  s.p_nf.value(),s.p_nris.value(),s.p_esp.value(),s.p_eta.value(),
                  s._pb(),s.p_kric.value(),s.p_rain.value(),s.p_perr.value(),
                  s.p_rdist.value())
        s.c_no.set(f"{r['pno']:.1f} dBm");s.c_w.set(f"{r['pw']:.1f} dBm")
        s.c_g.set(f"+{r['gain']:.1f} dB");s.c_snr.set(f"+{r['snr_w']-r['snr_no']:.1f} dB")
        s.c_nf.set(f"{r['Nfl']:.1f} dBm");s.c_dop.set(f"{r['dop']/1e3:.1f} kHz")
        draw_dia(s.axd,r,s.p_nris.value(),s.p_el.value(),s.p_alt.value())
        s.fdia.tight_layout();s.cvd.draw()
        draw_bars(s.axb,r);s.fbar.tight_layout();s.cvb.draw()

if __name__=="__main__":
    app=QApplication(sys.argv);w=Win();w.show();sys.exit(app.exec())
