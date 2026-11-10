# Figure contract

Core question: How do calculated W, Co and Fe spectra differ among explicit, relaxed coordination/adsorption models? The figures test candidate-dependent responses; no desired direction of spectral change is assumed.

Archetype: quantitative grid. Primary output: three vertical rows, W L3 / Co K / Fe K, each overlaying three structural states. The same state has the same colour/line style in every panel. Additional FT, k3chi and WT panels expose the structural basis of the comparison.

Backend: Python/matplotlib only. White background; Arial when installed with documented fallback; editable PDF/SVG and high-resolution PNG/TIFF. No manually generated peaks, smoothing, added noise, or individual energy alignment.

Source contract: completed native FEFF output, validated hashes and full intended absorber-site averages. TXT contains axis units, source identity, native energy references, signed chi and complex transforms. Complete native and processed outputs remain available beside the publication figures.

Scientific risks: static neutral-cell models are not OCP/1.6V electrode states; FT peak R is not a directly fitted bond length; FT/WT are derived from the same chi; a first-shell distance cutoff does not determine bond order; element-specific spectral change does not establish an exclusive catalytic site. DFT/QE charge density is not transferred to FEFF by this workflow.

Statistics: deterministic calculations, no experimental replicates or artificial error bars. Site averaging is not a statistical uncertainty estimate. Any future disorder, charge-state mixture or potential population must have independently specified weights and provenance.

QA: verify actual W coverage10100-10300eV before plotting, no missing comparison states, no false full-site averages, no axis/label overlap, actual font recorded, TXT numeric values consistent with source files, and adverse candidate outcomes retained in audits.
