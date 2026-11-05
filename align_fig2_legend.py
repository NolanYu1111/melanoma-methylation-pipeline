import re

fig2_exact_legend_md = (
    "**Fig. 2 Multi-panel Kaplan-Meier overall survival curves across top prognostic promoter CpG loci (n = 395).** "
    "Survival curves comparing high (blue solid line and shaded 95% confidence interval) versus low (red solid line and shaded 95% confidence interval) promoter methylation groups stratified at median beta-value. "
    "(A) *ZCCHC11* (`cg00754357`): High methylation extended overall survival (median: 9.76 vs. 4.49 yr; log-rank p = 4.52 × 10⁻⁵). Silencing this oncogene restores let-7 microRNA maturation. "
    "(B) *KMO* (`cg00606312`): High methylation prolonged survival (median: 8.92 vs. 5.09 yr; log-rank p = 1.00 × 10⁻⁴). Silencing *KMO* prevents immunosuppressive tryptophan degradation. "
    "(C) *EPHA10* (`cg00667938`): High methylation improved survival (median: 8.60 vs. 5.28 yr; log-rank p = 3.70 × 10⁻³) by dampening oncogenic receptor tyrosine kinase signaling. "
    "(D) *RGS21* (`cg00158654`): High methylation reduced survival (median: 5.23 vs. 9.37 yr; log-rank p = 5.20 × 10⁻³). Hypermethylation silences this GPCR-dampening tumor suppressor. "
    "(E) *CYMP* (`cg00076797`): High methylation shortened survival (median: 5.28 vs. 9.25 yr; log-rank p = 8.00 × 10⁻³) by silencing tumor suppressor regulation. "
    "(F) *AMPD2* (`cg00177290`): High methylation reduced survival (median: 5.23 vs. 9.37 yr; log-rank p = 2.78 × 10⁻²). Silencing *AMPD2* rewires purine pools, driving GTP accumulation."
)

fig2_exact_legend_html = (
    "<strong>Fig. 2</strong> Multi-panel Kaplan-Meier overall survival curves across top prognostic promoter CpG loci (n = 395). "
    "Survival curves comparing high (blue solid line and shaded 95% confidence interval) versus low (red solid line and shaded 95% confidence interval) promoter methylation groups stratified at median beta-value. "
    "(A) <em>ZCCHC11</em> (<code>cg00754357</code>): High methylation extended overall survival (median: 9.76 vs. 4.49 yr; log-rank <em>p</em> = 4.52 &times; 10<sup>-5</sup>). Silencing this oncogene restores <em>let-7</em> microRNA maturation. "
    "(B) <em>KMO</em> (<code>cg00606312</code>): High methylation prolonged survival (median: 8.92 vs. 5.09 yr; log-rank <em>p</em> = 1.00 &times; 10<sup>-4</sup>). Silencing <em>KMO</em> prevents immunosuppressive tryptophan degradation. "
    "(C) <em>EPHA10</em> (<code>cg00667938</code>): High methylation improved survival (median: 8.60 vs. 5.28 yr; log-rank <em>p</em> = 3.70 &times; 10<sup>-3</sup>) by dampening oncogenic receptor tyrosine kinase signaling. "
    "(D) <em>RGS21</em> (<code>cg00158654</code>): High methylation reduced survival (median: 5.23 vs. 9.37 yr; log-rank <em>p</em> = 5.20 &times; 10<sup>-3</sup>). Hypermethylation silences this GPCR-dampening tumor suppressor. "
    "(E) <em>CYMP</em> (<code>cg00076797</code>): High methylation shortened survival (median: 5.28 vs. 9.25 yr; log-rank <em>p</em> = 8.00 &times; 10<sup>-3</sup>) by silencing tumor suppressor regulation. "
    "(F) <em>AMPD2</em> (<code>cg00177290</code>): High methylation reduced survival (median: 5.23 vs. 9.37 yr; log-rank <em>p</em> = 2.78 &times; 10<sup>-2</sup>). Silencing <em>AMPD2</em> rewires purine pools, driving GTP accumulation."
)

# 1. Update manuscript.md
with open('manuscript.md', 'r', encoding='utf-8') as f:
    md = f.read()

md = re.sub(r'\*\*Fig\. 2 Multi-panel Kaplan-Meier overall survival curves.*?(?=\n\n|\n\*\*Key|\n---)', fig2_exact_legend_md, md, flags=re.DOTALL)

with open('manuscript.md', 'w', encoding='utf-8') as f:
    f.write(md)

# 2. Update manuscript.html
with open('manuscript.html', 'r', encoding='utf-8') as f:
    html = f.read()

html = re.sub(r'<div class="fig-caption"><strong>Fig\. 2</strong>.*?</div>', f'<div class="fig-caption">{fig2_exact_legend_html}</div>', html, flags=re.DOTALL)

with open('manuscript.html', 'w', encoding='utf-8') as f:
    f.write(html)

print("Aligned Figure 2 legend exactly across manuscript.md and manuscript.html!")
