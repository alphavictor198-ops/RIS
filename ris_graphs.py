#!/usr/bin/env python3
"""RIS Satellite Link — Interactive Graphs with Radio Buttons + Editable Parameters"""
import sys,subprocess
for p in["PyQt6","numpy","matplotlib","scipy"]:
    try:__import__(p.replace("-","_"))
    except:subprocess.check_call([sys.executable,"-m","pip","install",p,"-q"])
import numpy as np
from scipy.special import erfc
import matplotlib;matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FC
from matplotlib.figure import Figure
from PyQt6.QtWidgets import *
from PyQt6.QtCore import Qt,QTimer
C=3e8
BG="#f4f4f6";CD="#ffffff";PN="#eaeaee";LN="#d0d0d8"
CY="#0891b2";TX="#1e1e2e";DM="#6b7280"
RD="#dc2626";GN="#16a34a";BL="#2563eb";AM="#d97706"
CSS=f"""
QMainWindow,QWidget{{background:{BG};color:{TX};font-family:'Segoe UI',sans-serif;font-size:13px;}}
QLabel{{color:{TX};background:transparent;}}
QRadioButton{{color:{TX};font-size:12px;font-weight:600;padding:4px 0;background:transparent;}}
QRadioButton::indicator{{width:14px;height:14px;border-radius:7px;border:2px solid {LN};background:#fff;}}
QRadioButton::indicator:checked{{background:{CY};border-color:{CY};}}
QRadioButton:hover{{color:{CY};}}
QDoubleSpinBox,QSpinBox{{background:#fff;color:{BL};border:1px solid {LN};border-radius:4px;padding:3px 6px;font-size:12px;font-weight:700;min-width:72px;}}
QDoubleSpinBox:focus,QSpinBox:focus{{border-color:{CY};}}
QDoubleSpinBox::up-button,QDoubleSpinBox::down-button,QSpinBox::up-button,QSpinBox::down-button{{background:{PN};border:none;width:16px;border-radius:2px;}}
QDoubleSpinBox::up-arrow,QSpinBox::up-arrow{{image:none;width:0;height:0;border-left:3px solid transparent;border-right:3px solid transparent;border-bottom:4px solid {CY};margin:auto;}}
QDoubleSpinBox::down-arrow,QSpinBox::down-arrow{{image:none;width:0;height:0;border-left:3px solid transparent;border-right:3px solid transparent;border-top:4px solid {CY};margin:auto;}}
QComboBox{{background:#fff;color:{BL};border:1px solid {LN};border-radius:4px;padding:3px 6px;font-size:12px;font-weight:700;min-width:72px;}}
QComboBox QAbstractItemView{{background:#fff;color:{BL};selection-background-color:{CY};selection-color:#fff;}}
QScrollArea{{border:none;background:{PN};}}
"""
def slant_range(a,e):
    R=6371e3;h=a*1e3;t=np.deg2rad(e);return -R*np.sin(t)+np.sqrt((R*np.sin(t))**2+2*R*h+h*h)
def rain_atten(fg,rr,e):
    _f=[1,2,4,8,10,15,20,25,30,50,100]
    _k=[2.6e-5,1.54e-4,6.5e-4,4.54e-3,1.01e-2,3.67e-2,7.51e-2,.124,.187,.536,1.31]
    _a=[.97,1.07,1.12,1.33,1.26,1.15,1.10,1.06,1.02,.86,.74]
    return np.interp(fg,_f,_k)*rr**np.interp(fg,_f,_a)*3.0/np.sin(np.deg2rad(max(e,5)))
def pql(b):return 1.0 if b==0 else np.sinc(1.0/2**b)**2
def comp(ptx,fghz,bw,alt,el,gtx,grx,nf,nris,esp,eta,pb,kric,rr,pe,rd):
    f=fghz*1e9;lam=C/f;d=slant_range(alt,el)
    Lfs=20*np.log10(4*np.pi*d*f/C);La=0.5/np.sin(np.deg2rad(max(el,1)))
    Lr=rain_atten(fghz,rr,el) if rr>0 else 0;Nfl=-174+10*np.log10(bw*1e6)+nf
    pno=ptx+gtx+grx-Lfs-La-Lr;de=esp*lam;A=de**2;eq=pql(pb)
    Ns=max(int(np.sqrt(nris)),1);bwd=np.rad2deg(lam/(Ns*de))
    Lp=min(12*(pe/max(bwd,.01))**2,30);Kl=10**(kric/10);rc=Kl/(Kl+1)
    pr=(ptx+gtx+grx+20*np.log10(nris*A)+20*np.log10(eta)+10*np.log10(eq)
        +10*np.log10(rc)-Lp-10*np.log10(64*np.pi**3)-20*np.log10(d)-20*np.log10(rd)-La-Lr)
    pw=10*np.log10((np.sqrt(10**(pno/10))+np.sqrt(10**(pr/10)))**2)
    vo=np.sqrt(3.986e14/((6371+alt)*1e3));dp=f*vo*np.cos(np.deg2rad(el))/C
    return dict(lam=lam,d=d,Lfs=Lfs,Nfl=Nfl,pno=pno,pw=pw,gain=pw-pno,sno=pno-Nfl,sw=pw-Nfl,dp=dp,vo=vo)
# ─── Defaults per graph (ensure both lines visible) ───
DEFAULTS=dict(ptx=33,fghz=2.0,bw=10,alt=550,el=10,gtx=6,grx=10,nf=3,
              nris=100,esp=0.5,eta=0.9,pb=0,kric=15,rr=0,pe=0,rd=50)
# ─── UI ───
def _sep():
    f=QFrame();f.setFrameShape(QFrame.Shape.HLine);f.setStyleSheet(f"background:{LN};max-height:1px;border:none;");return f
def _lb(t,s=12,c=TX,b=False):
    w=QLabel(t);w.setStyleSheet(f"color:{c};font-size:{s}px;font-weight:{'700' if b else '400'};background:transparent;");return w
class FS(QWidget):
    def __init__(s,lab,lo,hi,v,st=1.,dc=1,cb=None):
        super().__init__();s.setFixedHeight(26);s._df=v;r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch();s.b=QDoubleSpinBox()
        s.b.setRange(lo,hi);s.b.setValue(v);s.b.setSingleStep(st);s.b.setDecimals(dc)
        s.b.setAlignment(Qt.AlignmentFlag.AlignRight);r.addWidget(s.b)
        if cb:s.b.valueChanged.connect(cb)
    def value(s):return s.b.value()
    def reset(s):s.b.blockSignals(True);s.b.setValue(s._df);s.b.blockSignals(False)
class IS(QWidget):
    def __init__(s,lab,lo,hi,v,cb=None):
        super().__init__();s.setFixedHeight(26);s._df=v;r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch();s.b=QSpinBox()
        s.b.setRange(lo,hi);s.b.setValue(v);s.b.setAlignment(Qt.AlignmentFlag.AlignRight);r.addWidget(s.b)
        if cb:s.b.valueChanged.connect(cb)
    def value(s):return s.b.value()
    def reset(s):s.b.blockSignals(True);s.b.setValue(s._df);s.b.blockSignals(False)
class CB(QWidget):
    def __init__(s,lab,items,idx=0,cb=None):
        super().__init__();s.setFixedHeight(26);s._df=idx;r=QHBoxLayout(s);r.setContentsMargins(0,0,0,0);r.setSpacing(4)
        r.addWidget(_lb(lab,11,DM));r.addStretch();s.b=QComboBox();s.b.addItems(items)
        s.b.setCurrentIndex(idx);s.b.setFixedWidth(85);r.addWidget(s.b)
        if cb:s.b.currentIndexChanged.connect(cb)
    def text(s):return s.b.currentText()
    def reset(s):s.b.blockSignals(True);s.b.setCurrentIndex(s._df);s.b.blockSignals(False)
# ─── Graphs (FIXED axis scales) ───
def sty(ax):
    ax.set_facecolor("#fafafa");ax.tick_params(colors=TX,labelsize=10)
    ax.grid(True,color="#e5e5ea",lw=.6,ls="--");ax.set_axisbelow(True)
    for s in ax.spines.values():s.set_color("#d1d5db")
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
def g_ber(ax,P):
    ax.clear();sty(ax);r=comp(**P);g=r['gain']
    snr=np.linspace(-10,25,300);sl=10**(snr/10)
    ax.semilogy(snr,np.clip(.5*erfc(np.sqrt(sl)),1e-8,1),color=RD,lw=2.2,label='Without RIS')
    ax.semilogy(snr,np.clip(.5*erfc(np.sqrt(sl*10**(g/10))),1e-8,1),color=CY,lw=2.2,label='With RIS')
    ax.axvline(r['sno'],color=RD,ls=':',lw=1.2,alpha=.6)
    ax.axvline(r['sw'],color=CY,ls=':',lw=1.2,alpha=.6)
    ax.axhline(1e-3,color=AM,ls=':',lw=1);ax.text(24,1.5e-3,'BER=10⁻³',color=AM,fontsize=9,ha='right')
    ax.set_xlim(-10,25);ax.set_ylim(1e-7,1)
    ax.set_xlabel('SNR (dB)',fontsize=11);ax.set_ylabel('Bit Error Rate',fontsize=11)
    ax.set_title(f'BER vs SNR  (RIS gain ≈ {g:.1f} dB)',fontsize=13,fontweight='bold',color=TX,pad=10)
    ax.legend(fontsize=10,facecolor=CD,edgecolor=LN)
def g_snr_n(ax,P):
    ax.clear();sty(ax)
    ns=np.array([4,16,36,64,100,200,400,600,800,1000,1500,2000]);sno=[];sw=[]
    for n in ns:pp=P.copy();pp['nris']=int(n);r=comp(**pp);sno.append(r['sno']);sw.append(r['sw'])
    ax.plot(ns,sno,color=RD,lw=2.2,marker='o',ms=5,label='Without RIS')
    ax.plot(ns,sw,color=CY,lw=2.2,marker='s',ms=5,label='With RIS')
    ax.fill_between(ns,sno,sw,color=CY,alpha=.08)
    ax.axvline(P['nris'],color=AM,ls='--',lw=1.2,label=f"Current N={P['nris']}")
    ax.set_xlim(0,2100);ax.set_ylim(-40,40)
    ax.set_xlabel('Number of RIS Elements (N)',fontsize=11);ax.set_ylabel('SNR (dB)',fontsize=11)
    ax.set_title('SNR vs RIS Elements',fontsize=13,fontweight='bold',color=TX,pad=10)
    ax.legend(fontsize=10,facecolor=CD,edgecolor=LN)
def g_pl(ax,P):
    ax.clear();sty(ax)
    fs=np.linspace(0.5,100,300);pn=[];pw=[]
    for f in fs:pp=P.copy();pp['fghz']=f;r=comp(**pp);pn.append(-r['pno']);pw.append(-r['pw'])
    ax.plot(fs,pn,color=RD,lw=2.2,label='Total loss (No RIS)')
    ax.plot(fs,pw,color=CY,lw=2.2,label='Total loss (With RIS)')
    ax.fill_between(fs,pn,pw,color=CY,alpha=.08)
    ax.axvline(P['fghz'],color=AM,ls='--',lw=1.2,label=f"Current f={P['fghz']} GHz")
    ax.set_xlim(0,100);ax.set_ylim(100,220)
    ax.set_xlabel('Carrier Frequency (GHz)',fontsize=11);ax.set_ylabel('Effective Path Loss (dB)',fontsize=11)
    ax.set_title('Path Loss vs Frequency',fontsize=13,fontweight='bold',color=TX,pad=10)
    ax.legend(fontsize=10,facecolor=CD,edgecolor=LN,loc='upper left')
def g_dop(ax,P):
    ax.clear();sty(ax)
    r0=comp(**P);vo=r0['vo'];f=P['fghz']*1e9
    t=np.linspace(0,600,500);el=5+85*np.sin(np.pi*t/600)
    dp=f*vo*np.cos(np.deg2rad(el))/C/1e3
    ax.plot(t,dp,color=RD,lw=2.2,label='Doppler (No RIS)')
    ax.plot(t,dp*.95,color=CY,lw=2,ls='--',label='Doppler (RIS phase-tracked)')
    ax.fill_between(t,dp,dp*.95,color=CY,alpha=.08)
    ax2=ax.twinx();ax2.plot(t,el,color=AM,lw=1.2,ls=':',alpha=.5,label='Elevation')
    ax2.set_ylabel('Elevation (°)',color=AM,fontsize=10);ax2.tick_params(colors=AM,labelsize=9)
    ax2.set_ylim(0,100);ax2.spines['right'].set_color(AM)
    h1,l1=ax.get_legend_handles_labels();h2,l2=ax2.get_legend_handles_labels()
    ax.legend(h1+h2,l1+l2,fontsize=9,facecolor=CD,edgecolor=LN,loc='upper right')
    ax.set_xlim(0,600);ax.set_ylim(0,None)
    ax.set_xlabel('Time (s)',fontsize=11);ax.set_ylabel('Doppler Shift (kHz)',fontsize=11)
    ax.set_title('Doppler Shift During Satellite Pass',fontsize=13,fontweight='bold',color=TX,pad=10)
def g_rain(ax,P):
    ax.clear();sty(ax)
    rrs=np.linspace(0,150,200);sno=[];sw=[]
    for rr in rrs:pp=P.copy();pp['rr']=rr;r=comp(**pp);sno.append(r['sno']);sw.append(r['sw'])
    ax.plot(rrs,sno,color=RD,lw=2.2,label='Without RIS')
    ax.plot(rrs,sw,color=CY,lw=2.2,label='With RIS')
    ax.fill_between(rrs,sno,sw,color=CY,alpha=.08)
    ax.axhline(0,color=AM,ls=':',lw=1);ax.text(145,.8,'SNR=0 dB',color=AM,fontsize=9,ha='right')
    ax.set_xlim(0,150);ax.set_ylim(-60,30)
    ax.set_xlabel('Rain Rate (mm/h)',fontsize=11);ax.set_ylabel('SNR (dB)',fontsize=11)
    ax.set_title('Weather/Rain Impact on SNR',fontsize=13,fontweight='bold',color=TX,pad=10)
    ax.legend(fontsize=10,facecolor=CD,edgecolor=LN)
GRAPHS=[("BER vs SNR",g_ber),("SNR vs RIS Elements",g_snr_n),("Path Loss vs Frequency",g_pl),
        ("Doppler Shift vs Time",g_dop),("Weather/Rain Impact",g_rain)]
# ─── Window ───
class Win(QMainWindow):
    def __init__(s):
        super().__init__();s.setWindowTitle("RIS Link — Analysis Graphs");s.setMinimumSize(1150,720);s.setStyleSheet(CSS)
        s._prev_idx=-1
        root=QWidget();s.setCentralWidget(root)
        ml=QHBoxLayout(root);ml.setSpacing(0);ml.setContentsMargins(0,0,0,0)
        scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setFixedWidth(280)
        scroll.setStyleSheet(f"background:{PN};border-right:1px solid {LN};")
        side=QWidget();sl=QVBoxLayout(side);sl.setContentsMargins(12,10,12,10);sl.setSpacing(4)
        up=lambda _:s._au();gswitch=lambda _:s._on_graph_switch()
        sl.addWidget(_lb("Graph Analysis",13,CY,True));sl.addWidget(_sep())
        sl.addWidget(_lb("▸ Select Graph",10,DM,True))
        s.radios=[]
        for i,(nm,_) in enumerate(GRAPHS):
            rb=QRadioButton(nm);rb.setChecked(i==0);rb.toggled.connect(gswitch);s.radios.append(rb);sl.addWidget(rb)
        sl.addWidget(_sep());sl.addWidget(_lb("▸ Link Parameters",10,DM,True))
        s.p_freq=FS("Carrier Freq (GHz)",.1,200,2,.5,1,up)
        s.p_bw=FS("Bandwidth (MHz)",.1,5000,10,1,1,up)
        s.p_ptx=FS("Tx Power (dBm)",0,60,33,1,1,up)
        for w in[s.p_freq,s.p_bw,s.p_ptx]:sl.addWidget(w)
        sl.addWidget(_sep());sl.addWidget(_lb("▸ Antennas",10,DM,True))
        s.p_gtx=FS("Tx Gain (dBi)",0,40,6,1,1,up)
        s.p_grx=FS("Rx Gain (dBi)",0,40,10,1,1,up)
        s.p_nf=FS("Noise Figure (dB)",0,15,3,.5,1,up)
        for w in[s.p_gtx,s.p_grx,s.p_nf]:sl.addWidget(w)
        sl.addWidget(_sep());sl.addWidget(_lb("▸ Satellite Orbit",10,DM,True))
        s.p_alt=FS("Altitude (km)",100,2000,550,50,0,up)
        s.p_el=FS("Elevation Angle (°)",5,90,10,1,1,up)
        for w in[s.p_alt,s.p_el]:sl.addWidget(w)
        sl.addWidget(_sep());sl.addWidget(_lb("▸ RIS Configuration",10,DM,True))
        s.p_nris=IS("Elements (N)",4,4096,100,up)
        s.p_esp=FS("Spacing (λ)",.1,1,.5,.1,1,up)
        s.p_eta=FS("Reflect. Efficiency (η)",.1,1,.9,.05,2,up)
        s.p_pb=CB("Phase Bits",["∞ (cont.)","2-bit","4-bit"],0,up)
        s.p_rd=FS("RIS→GS Dist (m)",1,500,50,5,0,up)
        for w in[s.p_nris,s.p_esp,s.p_eta,s.p_pb,s.p_rd]:sl.addWidget(w)
        sl.addWidget(_sep());sl.addWidget(_lb("▸ Channel & Environment",10,DM,True))
        s.p_kric=FS("Rician K-factor (dB)",0,30,15,1,1,up)
        s.p_rain=FS("Rain Rate (mm/h)",0,150,0,5,0,up)
        s.p_perr=FS("Pointing Error (°)",0,10,0,.1,1,up)
        for w in[s.p_kric,s.p_rain,s.p_perr]:sl.addWidget(w)
        s._all_spins=[s.p_freq,s.p_bw,s.p_ptx,s.p_gtx,s.p_grx,s.p_nf,s.p_alt,s.p_el,
                      s.p_nris,s.p_esp,s.p_eta,s.p_pb,s.p_rd,s.p_kric,s.p_rain,s.p_perr]
        sl.addStretch();scroll.setWidget(side);ml.addWidget(scroll)
        right=QWidget();right.setStyleSheet(f"background:{CD};")
        rl=QVBoxLayout(right);rl.setContentsMargins(10,8,10,8)
        s.fig=Figure(figsize=(10,7),facecolor="#fff");s.ax=s.fig.add_subplot(111)
        s.cv=FC(s.fig);s.cv.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Expanding)
        rl.addWidget(s.cv);ml.addWidget(right)
        s._tm=QTimer();s._tm.setSingleShot(True);s._tm.setInterval(150);s._tm.timeout.connect(s._run)
        QTimer.singleShot(50,s._run)
    def _au(s):s._tm.start()
    def _on_graph_switch(s):
        idx=next((i for i,r in enumerate(s.radios) if r.isChecked()),0)
        if idx!=s._prev_idx:
            for w in s._all_spins:w.reset()
            s._prev_idx=idx
        s._tm.start()
    def _pb(s):
        t=s.p_pb.text();return 0 if "∞" in t else int(t.split("-")[0])
    def _P(s):
        return dict(ptx=s.p_ptx.value(),fghz=s.p_freq.value(),bw=s.p_bw.value(),
            alt=s.p_alt.value(),el=s.p_el.value(),gtx=s.p_gtx.value(),grx=s.p_grx.value(),
            nf=s.p_nf.value(),nris=s.p_nris.value(),esp=s.p_esp.value(),eta=s.p_eta.value(),
            pb=s._pb(),kric=s.p_kric.value(),rr=s.p_rain.value(),pe=s.p_perr.value(),rd=s.p_rd.value())
    def _run(s):
        idx=next((i for i,r in enumerate(s.radios) if r.isChecked()),0)
        s._prev_idx=idx;GRAPHS[idx][1](s.ax,s._P());s.fig.tight_layout();s.cv.draw()
if __name__=="__main__":
    app=QApplication(sys.argv);w=Win();w.show();sys.exit(app.exec())
