#!/usr/bin/env python3
"""Coordinate-based preview: reference structure and extracted first shell."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]

def xyz(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    rows = [s.split() for s in lines[2:] if s.strip()]
    assert len(rows) == int(lines[0])
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])

def render():
    plt.rcParams.update({"font.family":"Arial","font.size":8,"svg.fonttype":"none","pdf.fonttype":42})
    fig = plt.figure(figsize=(7.1,3.65))
    colors = {"W":"#267A98","O":"#C94D4D","P":"#D6A950"}
    sizes = {"W":130,"O":50,"P":90}
    src = ROOT/"references"
    for panel, (filename,title) in enumerate([
            ("alpha_PW12O40_crystallographic.xyz",r"Complete $[\mathrm{PW}_{12}\mathrm{O}_{40}]^{3-}$ reference"),
            ("W1O6_first_shell_DIAGNOSTIC.xyz",r"Extracted W–O$_6$ first shell")]):
        elements, coords = xyz(src/filename)
        ax = fig.add_subplot(1,2,panel+1,projection="3d")
        for i, el in enumerate(elements):
            for j in range(i+1,len(elements)):
                pair = {el,elements[j]}
                distance = np.linalg.norm(coords[i]-coords[j])
                if (pair=={"W","O"} and distance<2.6) or (pair=={"P","O"} and distance<1.8):
                    ax.plot(*coords[[i,j]].T,color="#7B8388",lw=.65,zorder=1)
        for el in ("O","P","W"):
            selected = np.array([e==el for e in elements])
            if selected.any():
                ax.scatter(*coords[selected].T,s=sizes[el]*(1.9 if panel==1 else 1),
                    color=colors[el],edgecolor="white",linewidth=.4,depthshade=True,zorder=3)
        ax.view_init(elev=18,azim=30 if panel==0 else 30)
        center=(coords.min(0)+coords.max(0))/2
        half=max(np.ptp(coords,axis=0))*.55
        ax.set(xlim=(center[0]-half,center[0]+half),ylim=(center[1]-half,center[1]+half),zlim=(center[2]-half,center[2]+half))
        ax.set_box_aspect((1,1,1))
        ax.set_axis_off()
        ax.set_title(title,fontsize=9,pad=2)
        ax.text2D(.02,1.01,"ab"[panel],transform=ax.transAxes,weight="bold",fontsize=11)
        if panel==1:
            ax.text2D(.5,-.015,r"W–O: 1.714, 1.918 × 2, 1.942 × 2, 2.347 $\AA$",ha="center",transform=ax.transAxes,fontsize=7)
    for index,el in enumerate(("W","O","P")):
        fig.text(.40+index*.08,.105,"●",color=colors[el],fontsize=11,ha="center")
        fig.text(.418+index*.08,.105,el,fontsize=8)
    fig.text(.5,.035,"Coordinates from published crystallography; the seven-atom fragment is a geometric diagnostic only.",ha="center",fontsize=7)
    fig.subplots_adjust(left=.015,right=.985,bottom=.2,top=.89,wspace=.02)
    out = ROOT/"preview"
    out.mkdir(exist_ok=True)
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"W_O_models.{ext}",dpi=400)
    plt.close(fig)
    (out/"figure_contract.json").write_text(json.dumps({
        "claim":"The complete reference has O6 first-shell W coordination, but includes other W atoms beyond this shell.",
        "archetype":"coordinate-based structural comparison", "backend":"Python matplotlib",
        "source":"references/*.xyz", "is_spectrum":False,
        "not_supported_single_atom_catalyst":True,"size_inches":[7.1,3.65],"exports":["PNG","PDF","SVG"]},indent=2),encoding="utf-8")

if __name__=="__main__":
    render()
