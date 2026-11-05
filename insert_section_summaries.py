import re

# 1. Update manuscript.md
with open('manuscript.md', 'r', encoding='utf-8') as f:
    md = f.read()

# Insert summary 1
s1 = "\n\n**Key Results Achieved in this Section:** Established a quality-controlled data reduction pipeline that resolved longitudinal clinical record redundancy in TCGA-SKCM (recovering 395 fully matched patient profiles from 470 cases) and isolated 214 high-confidence promoter CpG sites for downstream prognostic and diagnostic modeling.\n"
md = md.replace("The analytical workflow and data filtering stages are summarized in Fig. 1 and Table 1.", "The analytical workflow and data filtering stages are summarized in Fig. 1 and Table 1." + s1)

# Insert summary 2
s2 = "\n\n**Key Results Achieved in this Section:** Identified a 24-CpG prognostic survival panel in melanoma exhibiting distinct epigenetic directionality: 15 loci function as silenced oncogenic drivers (protective hypermethylation), while 9 loci function as silenced tumor suppressors (hazardous hypermethylation).\n"
md = md.replace("hazardous silencing).\n\n**Table 3", "hazardous silencing)." + s2 + "\n**Table 3")

# Insert summary 3
s3 = "\n\n**Key Results Achieved in this Section:** Demonstrated complete epigenetic segregation between normal controls and melanoma tumors, proving near-universal promoter hypermethylation of tumor suppressors (mean $\\beta > 0.75$; $> 95\\%$ methylated) and profound promoter hypomethylation of oncogenes (mean $\\beta < 0.15$; $> 87\\%$ demethylated) in melanoma.\n"
md = md.replace("relative to methylated normal controls ($\\beta > 0.75$).\n\n### Diagnostic", "relative to methylated normal controls ($\\beta > 0.75$).\n" + s3 + "\n### Diagnostic")

# Insert summary 4
s4 = "\n\n**Key Results Achieved in this Section:** Built and validated six diagnostic machine learning classifiers on the primary training set (*n* = 790), achieving up to 92.15% 10-fold cross-validation accuracy (AUC-ROC = 0.9722) and identifying *AMPD2* (`cg00177290`) and *WNT2B* (`cg00255837`) as the top diagnostic descriptors.\n"
md = md.replace("without overfitting.\n\n### Independent External", "without overfitting.\n" + s4 + "\n### Independent External")

# Insert summary 5
s5 = "\n\n**Key Results Achieved in this Section:** Confirmed robust external generalizability of our diagnostic classifiers on unseen clinical tissue (GSE105191, *n* = 40), achieving 95.00% external validation accuracy, 100.00% sensitivity, and AUC-ROC = 0.9800 with zero overfitting.\n"
md = md.replace("Score = 0.0278).\n\n### Visualizing", "Score = 0.0278).\n" + s5 + "\n### Visualizing")

# Insert summary 6
s6 = "\n\n**Key Results Achieved in this Section:** Visually confirmed the clinical impact of epigenetic directionality across top prognostic loci, demonstrating that hypermethylation of oncogenes (*ZCCHC11*, *KMO*, *EPHA10*) significantly extends overall survival (up to 9.76 yr), whereas hypermethylation of tumor suppressors (*RGS21*, *CYMP*, *AMPD2*) significantly shortens overall survival (down to 5.23 yr).\n"
md = md.replace("rewiring purine pools.\n\n---", "rewiring purine pools.\n" + s6 + "\n---")

with open('manuscript.md', 'w', encoding='utf-8') as f:
    f.write(md)
print("Updated manuscript.md with section results summaries!")

# 2. Update manuscript.html
with open('manuscript.html', 'r', encoding='utf-8') as f:
    html = f.read()

hs1 = '<p><strong>Key Results Achieved in this Section:</strong> Established a quality-controlled data reduction pipeline that resolved longitudinal clinical record redundancy in TCGA-SKCM (recovering 395 fully matched patient profiles from 470 cases) and isolated 214 high-confidence promoter CpG sites for downstream prognostic and diagnostic modeling.</p>'
html = html.replace('<p>The analytical workflow and data filtering stages are summarized in Fig. 1 and Table 1.</p>', '<p>The analytical workflow and data filtering stages are summarized in Fig. 1 and Table 1.</p>\n        ' + hs1)

hs2 = '<p><strong>Key Results Achieved in this Section:</strong> Identified a 24-CpG prognostic survival panel in melanoma exhibiting distinct epigenetic directionality: 15 loci function as silenced oncogenic drivers (protective hypermethylation), while 9 loci function as silenced tumor suppressors (hazardous hypermethylation).</p>'
html = html.replace('<p>Log-rank testing identified 24 promoter CpG sites significantly associated with overall survival (<em>p</em> &lt; 0.05; Table 3). The loci segregated into 15 oncogenic loci where hypermethylation correlated with longer survival (protective silencing) and 9 tumor-suppressor loci where hypermethylation correlated with shorter survival (hazardous silencing).</p>', '<p>Log-rank testing identified 24 promoter CpG sites significantly associated with overall survival (<em>p</em> &lt; 0.05; Table 3). The loci segregated into 15 oncogenic loci where hypermethylation correlated with longer survival (protective silencing) and 9 tumor-suppressor loci where hypermethylation correlated with shorter survival (hazardous silencing).</p>\n        ' + hs2)

hs3 = '<p><strong>Key Results Achieved in this Section:</strong> Demonstrated complete epigenetic segregation between normal controls and melanoma tumors, proving near-universal promoter hypermethylation of tumor suppressors (mean &beta; &gt; 0.75; &gt; 95% methylated) and profound promoter hypomethylation of oncogenes (mean &beta; &lt; 0.15; &gt; 87% demethylated) in melanoma.</p>'
html = re.sub(r'(<div class="fig-caption"><strong>Fig\. 3</strong>.*?</div>\s*</div>)', r'\1\n\n    <div class="two-column-body">\n        ' + hs3 + '\n    </div>', html, flags=re.DOTALL)

hs4 = '<p><strong>Key Results Achieved in this Section:</strong> Built and validated six diagnostic machine learning classifiers on the primary training set (<em>n</em> = 790), achieving up to 92.15% 10-fold cross-validation accuracy (AUC-ROC = 0.9722) and identifying <em>AMPD2</em> (<code>cg00177290</code>) and <em>WNT2B</em> (<code>cg00255837</code>) as the top diagnostic descriptors.</p>'
html = re.sub(r'(<div class="fig-caption"><strong>Fig\. 5</strong>.*?</div>\s*</div>)', r'\1\n\n    <div class="two-column-body">\n        ' + hs4 + '\n    </div>', html, flags=re.DOTALL)

hs5 = '<p><strong>Key Results Achieved in this Section:</strong> Confirmed robust external generalizability of our diagnostic classifiers on unseen clinical tissue (GSE105191, <em>n</em> = 40), achieving 95.00% external validation accuracy, 100.00% sensitivity, and AUC-ROC = 0.9800 with zero overfitting.</p>'
html = re.sub(r'(<div class="fig-caption"><strong>Fig\. 6</strong>.*?</div>\s*</div>)', r'\1\n\n    <div class="two-column-body">\n        ' + hs5 + '\n    </div>', html, flags=re.DOTALL)

hs6 = '<p><strong>Key Results Achieved in this Section:</strong> Visually confirmed the clinical impact of epigenetic directionality across top prognostic loci, demonstrating that hypermethylation of oncogenes (<em>ZCCHC11</em>, <em>KMO</em>, <em>EPHA10</em>) significantly extends overall survival (up to 9.76 yr), whereas hypermethylation of tumor suppressors (<em>RGS21</em>, <em>CYMP</em>, <em>AMPD2</em>) significantly shortens overall survival (down to 5.23 yr).</p>'
html = re.sub(r'(<div class="fig-caption"><strong>Fig\. 2</strong>.*?</div>\s*</div>)', r'\1\n\n    <div class="two-column-body">\n        ' + hs6 + '\n    </div>', html, flags=re.DOTALL)

with open('manuscript.html', 'w', encoding='utf-8') as f:
    f.write(html)
print("Updated manuscript.html with section results summaries!")
