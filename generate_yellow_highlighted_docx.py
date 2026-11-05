import sys
import os
import re
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls, qn

def enable_track_changes(doc):
    """Enable Track Changes in Word document settings."""
    settings = doc.settings.element
    track_rev = settings.find(qn('w:trackRevisions'))
    if track_rev is None:
        track_rev = parse_xml(r'<w:trackRevisions %s/>' % nsdecls('w'))
        settings.append(track_rev)

def set_table_wolters_kluwer_borders(table):
    """Apply Wolters Kluwer table borders: top line 1.5pt, header bottom line 1.0pt, bottom line 1.5pt, no vertical lines."""
    tblPr = table._tbl.tblPr
    borders_xml = parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        f'<w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        f'<w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        f'<w:insideH w:val="single" w:sz="4" w:space="0" w:color="E0E0E0"/>'
        f'<w:left w:val="none"/>'
        f'<w:right w:val="none"/>'
        f'<w:insideV w:val="none"/>'
        f'</w:tblBorders>'
    )
    tblPr.append(borders_xml)

def set_header_row_border(row):
    """Add thick bottom border under table header row."""
    trPr = row._tr.get_or_add_trPr()
    for cell in row.cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcBorders = parse_xml(
            f'<w:tcBorders {nsdecls("w")}>'
            f'<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
            f'</w:tcBorders>'
        )
        tcPr.append(tcBorders)

def add_hyperlink(paragraph, text, url, color="0000FF", underline=True):
    """Add a real Word XML hyperlink element to a paragraph."""
    part = paragraph.part
    r_id = part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True
    )
    
    hyperlink = parse_xml(f'<w:hyperlink {nsdecls("w", "r")} r:id="{r_id}" w:history="1"/>')
    
    rPr_xml = '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="22"/>'
    if color:
        rPr_xml += f'<w:color w:val="{color}"/>'
    if underline:
        rPr_xml += '<w:u w:val="single"/>'
    rPr_xml += '</w:rPr>'
    
    escaped_text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    r_xml = f'<w:r {nsdecls("w")}>{rPr_xml}<w:t xml:space="preserve">{escaped_text}</w:t></w:r>'
    
    hyperlink.append(parse_xml(r_xml))
    paragraph._p.append(hyperlink)

def parse_inline_markdown(paragraph, text, default_font_name='Times New Roman', default_font_size=11, highlight_yellow=False):
    text = text.replace('<sup>', '').replace('</sup>', '').replace('&nbsp;', ' ')
    text = re.sub(r'\$([^$]+)\$', r'\1', text)
    text = text.replace(r'\beta', 'β').replace(r'\Delta', 'Δ').replace(r'\kappa', 'κ').replace(r'\in', '∈').replace(r'\ge', '≥').replace(r'\le', '≤').replace(r'\text{', '').replace(r'}', '').replace(r'\\', '')
    
    pattern = r'(\[.*?\]\(https?://[^\s\)]+\)|https?://[^\s\)]+|\*\*\*.*?\*\*\*|\*\*.*?\*\*|\*.*?\*|__.*?__|__.*?__|_.*?_)'
    tokens = re.split(pattern, text)
    
    for token in tokens:
        if not token:
            continue
            
        md_link_match = re.match(r'^\[(.*?)\]\((https?://[^\s\)]+)\)$', token)
        if md_link_match:
            link_text = md_link_match.group(1)
            link_url = md_link_match.group(2)
            add_hyperlink(paragraph, link_text, link_url)
            continue
            
        if token.startswith('http://') or token.startswith('https://'):
            clean_url = token.rstrip('.,;:')
            add_hyperlink(paragraph, clean_url, clean_url)
            continue
        
        is_bold = False
        is_italic = False
        clean_text = token
        
        if token.startswith('***') and token.endswith('***'):
            is_bold = True
            is_italic = True
            clean_text = token[3:-3]
        elif token.startswith('**') and token.endswith('**'):
            is_bold = True
            clean_text = token[2:-2]
        elif token.startswith('__') and token.endswith('__'):
            is_bold = True
            clean_text = token[2:-2]
        elif token.startswith('*') and token.endswith('*'):
            is_italic = True
            clean_text = token[1:-1]
        elif token.startswith('_') and token.endswith('_'):
            is_italic = True
            clean_text = token[1:-1]
            
        run = paragraph.add_run(clean_text)
        run.font.name = default_font_name
        run.font.size = Pt(default_font_size)
        run.bold = is_bold
        run.italic = is_italic
        if highlight_yellow:
            run.font.highlight_color = WD_COLOR_INDEX.YELLOW

def split_table_row(row_str):
    cleaned = row_str.strip()
    if cleaned.startswith('|'):
        cleaned = cleaned[1:]
    if cleaned.endswith('|'):
        cleaned = cleaned[:-1]
    parts = cleaned.split(' | ')
    return [p.strip() for p in parts]

def add_tracked_revision(paragraph, deleted_text, inserted_text, author="Editor", rev_id=1):
    """Add a tracked deletion and insertion in Word XML format."""
    p_elem = paragraph._p
    
    if deleted_text:
        del_xml = parse_xml(
            f'<w:del {nsdecls("w")} w:id="{rev_id}" w:author="{author}" w:date="2026-08-23T12:00:00Z">'
            f'<w:r>'
            f'<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="22"/></w:rPr>'
            f'<w:delText xml:space="preserve">{deleted_text}</w:delText>'
            f'</w:r>'
            f'</w:del>'
        )
        p_elem.append(del_xml)
        
    if inserted_text:
        ins_xml = parse_xml(
            f'<w:ins {nsdecls("w")} w:id="{rev_id+1}" w:author="{author}" w:date="2026-08-23T12:00:00Z">'
            f'<w:r>'
            f'<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="22"/></w:rPr>'
            f'<w:t xml:space="preserve">{inserted_text}</w:t>'
            f'</w:r>'
            f'</w:ins>'
        )
        p_elem.append(ins_xml)

def convert_md_to_perfect_docx(md_path, docx_path):
    print(f"Reading Markdown from {md_path}...")
    with open(md_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    doc = Document()
    enable_track_changes(doc)
    
    # Configure Page Setup (Single Column, 1 in margins, Times New Roman 11pt, Double Spacing)
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(1.0)
    sec.bottom_margin = Inches(1.0)
    sec.left_margin = Inches(1.0)
    sec.right_margin = Inches(1.0)
    
    # Add Article Title (Times New Roman 16pt Bold)
    p_title = doc.add_paragraph()
    p_title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_title.paragraph_format.line_spacing = 2.0
    p_title.paragraph_format.space_before = Pt(0)
    p_title.paragraph_format.space_after = Pt(12)
    r_title = p_title.add_run("Diagnostics and prognostics of melanoma using methylation profiles")
    r_title.font.name = 'Times New Roman'
    r_title.font.size = Pt(16)
    r_title.bold = True
    
    # Add Authors Line (Times New Roman 11pt)
    p_author = doc.add_paragraph()
    p_author.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p_author.paragraph_format.line_spacing = 2.0
    p_author.paragraph_format.space_before = Pt(0)
    p_author.paragraph_format.space_after = Pt(12)
    
    authors_data = [
        ("Nolan Yu", "1", ", "),
        ("Valentina L. Kouznetsova", "1,2", ", "),
        ("Santosh Kesari", "3", ", and "),
        ("Igor F. Tsigelny", "1,2,4*", "")
    ]
    for name, sup, sep in authors_data:
        r_name = p_author.add_run(name)
        r_name.font.name = 'Times New Roman'
        r_name.font.size = Pt(11)
        r_sup = p_author.add_run(sup)
        r_sup.font.name = 'Times New Roman'
        r_sup.font.size = Pt(11)
        r_sup.font.superscript = True
        r_sep = p_author.add_run(sep)
        r_sep.font.name = 'Times New Roman'
        r_sep.font.size = Pt(11)
        
    # Add Affiliations Block
    p_aff = doc.add_paragraph()
    p_aff.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p_aff.paragraph_format.line_spacing = 1.5
    p_aff.paragraph_format.space_before = Pt(0)
    p_aff.paragraph_format.space_after = Pt(12)
    
    affiliations = [
        ("1", "Biana Institute, San Diego, CA, United States"),
        ("2", "San Diego Supercomputer Center, University of California, San Diego, La Jolla, CA, United States"),
        ("3", "Lundquist Institute, Torrance, CA, United States"),
        ("4", "Department of Neurosciences, University of California San Diego, La Jolla, CA, United States")
    ]
    for idx, (sup, text) in enumerate(affiliations):
        r_sup = p_aff.add_run(sup)
        r_sup.font.name = 'Times New Roman'
        r_sup.font.size = Pt(10)
        r_sup.font.superscript = True
        r_text = p_aff.add_run(text + ("\n" if idx < len(affiliations) - 1 else ""))
        r_text.font.name = 'Times New Roman'
        r_text.font.size = Pt(10)
        
    # Add Correspondence Block
    p_corr = doc.add_paragraph()
    p_corr.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p_corr.paragraph_format.line_spacing = 1.5
    p_corr.paragraph_format.space_before = Pt(0)
    p_corr.paragraph_format.space_after = Pt(18)
    
    r_ast = p_corr.add_run("*")
    r_ast.font.name = 'Times New Roman'
    r_ast.font.size = Pt(10)
    r_ast.font.superscript = True
    
    r_lbl = p_corr.add_run("Correspondence: ")
    r_lbl.font.name = 'Times New Roman'
    r_lbl.font.size = Pt(10)
    r_lbl.bold = True
    
    r_name = p_corr.add_run("Igor F. Tsigelny, ")
    r_name.font.name = 'Times New Roman'
    r_name.font.size = Pt(10)
    
    add_hyperlink(p_corr, "itsigeln@ucsd.edu", "mailto:itsigeln@ucsd.edu")
    
    in_code_block = False
    i = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        if line == '---':
            i += 1
            continue
            
        # Skip top metadata lines from Markdown as rendered above
        if (line.startswith('# ') or 
            line.startswith('**Nolan Yu') or 
            line.startswith('Nolan Yu') or 
            'Biana Institute' in line or 
            'San Diego Supercomputer' in line or 
            'Lundquist Institute' in line or 
            'Department of Neurosciences' in line or 
            'Correspondence:' in line or 
            line.startswith('<sup>') or 
            line.startswith('¹') or line.startswith('²') or line.startswith('³') or line.startswith('⁴')):
            i += 1
            continue
            
        if line.startswith('```'):
            in_code_block = not in_code_block
            i += 1
            continue
            
        if in_code_block:
            raw_line = lines[i].rstrip('\n')
            p = doc.add_paragraph()
            p.paragraph_format.first_line_indent = Inches(0)
            p.paragraph_format.left_indent = Inches(0.3)
            p.paragraph_format.line_spacing = 1.15
            p.paragraph_format.space_after = Pt(2)
            run = p.add_run(raw_line)
            run.font.name = 'Consolas'
            run.font.size = Pt(8.5)
            i += 1
            continue
            
        # Display Math Equations (\[ ... \] or $$ ... $$)
        if line.startswith(r'\[') or line.startswith('$$'):
            math_lines = []
            if line.endswith(r'\]') or (len(line) > 2 and line.endswith('$$')):
                math_content = line.strip(r'\[\]$ ')
            else:
                math_lines.append(line.strip(r'\[\]$ '))
                i += 1
                while i < len(lines):
                    m_line = lines[i].strip()
                    if m_line.endswith(r'\]') or m_line.endswith('$$'):
                        math_lines.append(m_line.strip(r'\[\]$ '))
                        break
                    math_lines.append(m_line)
                    i += 1
                math_content = " ".join(math_lines)
            
            p_eq = doc.add_paragraph()
            p_eq.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_eq.paragraph_format.line_spacing = 2.0
            p_eq.paragraph_format.space_before = Pt(12)
            p_eq.paragraph_format.space_after = Pt(12)
            
            norm = math_content.replace('\\', '').replace(' ', '')
            
            if "t=" in norm or "Delta" in norm:
                r_eq = p_eq.add_run("t = |ΔBeta| / √([RMSD_Diseased² / n] + [RMSD_Normal² / n])")
            elif "RMSD" in norm and "ge" in norm:
                r_eq = p_eq.add_run("RMSD_Diseased ≥ 0.5 × |ΔBeta|    and    RMSD_Normal ≥ 0.5 × |ΔBeta|")
            elif "IG(" in norm or "entropy" in norm:
                r_eq = p_eq.add_run("IG(Y, X) = H(Y) - H(Y | X) = ∑ P(x, y) log₂( P(x, y) / [ P(x) P(y) ] )")
            elif "SplitInfo" in norm:
                r_eq = p_eq.add_run("SplitInfo(X) = - ∑ (|S_i| / |S|) log₂(|S_i| / |S|),   GainRatio(Y, X) = IG(Y, X) / SplitInfo(X)")
            else:
                clean_m = math_content.replace(r'\frac', '').replace(r'\sqrt', '√').replace(r'\text{', '').replace(r'}', '').replace(r'\quad', '   ').replace(r'\ge', '≥')
                r_eq = p_eq.add_run(clean_m)
                
            r_eq.font.name = 'Times New Roman'
            r_eq.font.size = Pt(11)
            r_eq.italic = True
            i += 1
            continue
            
        # Section Headings (H2)
        if line.startswith('## '):
            heading_text = line[3:].strip()
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.line_spacing = 2.0
            run = p.add_run(heading_text)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(12)
            run.bold = True
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.keep_with_next = True
            i += 1
            continue
            
        # Subsection Headings (H3)
        if line.startswith('### '):
            heading_text = line[4:].strip()
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.line_spacing = 2.0
            run = p.add_run(heading_text)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(11)
            run.bold = True
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.keep_with_next = True
            i += 1
            continue
            
        # Table Titles / Captions
        if line.startswith('**Table '):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.keep_with_next = True
            parse_inline_markdown(p, line, default_font_name='Times New Roman', default_font_size=10)
            i += 1
            continue
            
        # Tables (Wolters Kluwer Horizontal Border Style)
        if line.startswith('|'):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                table_lines.append(lines[i].strip())
                i += 1
                
            if len(table_lines) >= 3:
                headers = split_table_row(table_lines[0])
                rows = []
                for row_line in table_lines[2:]:
                    row_data = split_table_row(row_line)
                    rows.append(row_data)
                    
                table = doc.add_table(rows=1, cols=len(headers))
                table.alignment = WD_TABLE_ALIGNMENT.CENTER
                set_table_wolters_kluwer_borders(table)
                
                hdr_trPr = table.rows[0]._tr.get_or_add_trPr()
                hdr_trPr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
                hdr_trPr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
                set_header_row_border(table.rows[0])
                
                hdr_cells = table.rows[0].cells
                for idx, h_text in enumerate(headers):
                    p = hdr_cells[idx].paragraphs[0]
                    p.text = ""
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER if idx > 0 else WD_ALIGN_PARAGRAPH.LEFT
                    p.paragraph_format.line_spacing = 1.15
                    p.paragraph_format.space_after = Pt(3)
                    p.paragraph_format.space_before = Pt(3)
                    clean_h = h_text.replace(r'\(p\)', '*p*').replace(r'\(', '').replace(r'\)', '').replace(r'\text{', '').replace(r'}', '').replace(r'\\', '')
                    parse_inline_markdown(p, clean_h, default_font_name='Times New Roman', default_font_size=9.5)
                    for r in p.runs:
                        r.bold = True
                        
                for r_data in rows:
                    row_obj = table.add_row()
                    trPr = row_obj._tr.get_or_add_trPr()
                    trPr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
                    row_cells = row_obj.cells
                    for idx in range(len(headers)):
                        val = r_data[idx] if idx < len(r_data) else ""
                        p = row_cells[idx].paragraphs[0]
                        p.text = ""
                        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if (len(val) < 12 or idx in [0, 2, 3, 4, 5]) else WD_ALIGN_PARAGRAPH.LEFT
                        p.paragraph_format.line_spacing = 1.15
                        p.paragraph_format.space_after = Pt(3)
                        p.paragraph_format.space_before = Pt(3)
                        
                        is_gene_symbol = (idx == 1 and headers[0] == "#") or val.startswith('*')
                        clean_val = val.replace('*', '').replace(r'\(', '').replace(r'\)', '').replace(r'\text{', '').replace(r'}', '').replace(r'\\', '').strip()
                        
                        if is_gene_symbol:
                            r = p.add_run(clean_val)
                            r.font.name = 'Times New Roman'
                            r.font.size = Pt(9.5)
                            r.italic = True
                        else:
                            parse_inline_markdown(p, clean_val, default_font_name='Times New Roman', default_font_size=9.5)
                            
                num_cols = len(headers)
                if num_cols == 4:
                    col_widths = [Inches(2.2), Inches(1.0), Inches(0.8), Inches(2.5)]
                elif num_cols == 7:
                    col_widths = [Inches(1.5), Inches(1.1), Inches(0.7), Inches(0.8), Inches(0.8), Inches(0.8), Inches(0.8)]
                elif num_cols == 8:
                    col_widths = [Inches(0.4), Inches(0.9), Inches(0.9), Inches(0.7), Inches(0.8), Inches(0.8), Inches(0.8), Inches(1.2)]
                elif num_cols == 5:
                    col_widths = [Inches(0.5), Inches(1.6), Inches(1.1), Inches(1.6), Inches(1.7)]
                elif num_cols == 6:
                    col_widths = [Inches(1.6), Inches(1.0), Inches(1.0), Inches(1.0), Inches(1.0), Inches(0.9)]
                else:
                    col_widths = [Inches(6.5 / num_cols)] * num_cols

                for row_obj_item in table.rows:
                    for c_idx, w_val in enumerate(col_widths):
                        if c_idx < len(row_obj_item.cells):
                            row_obj_item.cells[c_idx].width = w_val

                p_space = doc.add_paragraph()
                p_space.paragraph_format.line_spacing = 1.0
                p_space.paragraph_format.space_after = Pt(6)
            continue
            
        # Figures & Images
        if line.startswith('!['):
            img_match = re.match(r'!\[(.*?)\]\((.*?)\)', line)
            if img_match:
                alt_text = img_match.group(1)
                img_path = img_match.group(2)
                
                if os.path.exists(img_path):
                    p = doc.add_paragraph()
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    p.paragraph_format.space_before = Pt(12)
                    p.paragraph_format.space_after = Pt(4)
                    p.paragraph_format.line_spacing = 1.0
                    r = p.add_run()
                    r.add_picture(img_path, width=Inches(6.0))
                    
                    p_cap = doc.add_paragraph()
                    p_cap.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
                    p_cap.paragraph_format.line_spacing = 1.5
                    p_cap.paragraph_format.space_before = Pt(2)
                    p_cap.paragraph_format.space_after = Pt(12)
                    parse_inline_markdown(p_cap, alt_text, default_font_name='Times New Roman', default_font_size=10)
            i += 1
            continue

        # Lists / Bullet Points
        if line.startswith('* ') or line.startswith('- ') or re.match(r'^\d+\.\s', line):
            p = doc.add_paragraph()
            p.paragraph_format.line_spacing = 2.0
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.left_indent = Inches(0.25)
            parse_inline_markdown(p, line, default_font_name='Times New Roman', default_font_size=11)
            i += 1
            continue

        # Regular Body Paragraph (Strict Double Spacing, 2 Intervals Between Lines)
        if line:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 2.0 # Strict double spacing!
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.first_line_indent = Inches(0.5)
            
            # Tracked Revisions & Yellow Highlights Preserved
            if "In oncogenesis, promoter hypermethylation" in line:
                prefix = "In oncogenesis, promoter hypermethylation is a well-established mechanism for silencing tumor-suppressor genes, acting as one of the \"two hits\" in Knudson's hypothesis for gene inactivation. Conversely, "
                yellow_orig = "promoter hypomethylation can lead to the aberrant activation of oncogenes or cancer-testis antigens."
                new_rewrite_next = " (i.e., promoter hypomethylation can lead to the aberrant transcription and pathological activation of oncogenes or cancer-testis antigens)."
                suffix = " measurement of these methylation levels is quantified using \"beta values\" (β-values), which range from 0 (completely unmethylated) to 1 (fully methylated), representing the proportion of methylated alleles at a given locus across a cell population [6]."
                
                p_run1 = p.add_run(prefix)
                p_run1.font.name = 'Times New Roman'
                p_run1.font.size = Pt(11)
                
                p_run_yellow = p.add_run(yellow_orig)
                p_run_yellow.font.name = 'Times New Roman'
                p_run_yellow.font.size = Pt(11)
                p_run_yellow.font.highlight_color = WD_COLOR_INDEX.YELLOW
                
                p_run_rewrite = p.add_run(new_rewrite_next + " ")
                p_run_rewrite.font.name = 'Times New Roman'
                p_run_rewrite.font.size = Pt(11)
                
                add_tracked_revision(p, "Usually", "The", rev_id=10)
                
                p_run2 = p.add_run(suffix)
                p_run2.font.name = 'Times New Roman'
                p_run2.font.size = Pt(11)
                
            elif "In the current article we try to elucidate" in line:
                add_tracked_revision(p, "In the current article we try to elucidate if ", "")
                p_run1 = p.add_run("DNA methylation profiles ")
                p_run1.font.name = 'Times New Roman'
                p_run1.font.size = Pt(11)
                add_tracked_revision(p, "", "can ")
                p_run2 = p.add_run("be utilized to develop a robust, statistically significant ")
                p_run2.font.name = 'Times New Roman'
                p_run2.font.size = Pt(11)
                add_tracked_revision(p, "", "diagnostic and ")
                p_run3 = p.add_run("prognostic panel for cutaneous melanoma ")
                p_run3.font.name = 'Times New Roman'
                p_run3.font.size = Pt(11)
                add_tracked_revision(p, "survival", "")
                p_run4 = p.add_run(", and ")
                p_run4.font.name = 'Times New Roman'
                p_run4.font.size = Pt(11)
                add_tracked_revision(p, "", "how ")
                p_run5 = p.add_run("the identified epigenetic alterations ")
                p_run5.font.name = 'Times New Roman'
                p_run5.font.size = Pt(11)
                add_tracked_revision(p, "", "can ")
                p_run6 = p.add_run("be therapeutically targeted")
                p_run6.font.name = 'Times New Roman'
                p_run6.font.size = Pt(11)
                add_tracked_revision(p, "", ".")
                add_tracked_revision(p, " using translationally viable, pathway-level interventions given that CRISPR-dCas9 epigenome editing is not yet approved for clinical practice?", "")
                
            elif "To eliminate multi-line record redundancy" in line or "To resolve a critical clinical data processing bug" in line:
                orig_dedup_yellow = (
                    "To resolve a critical clinical data processing bug that caused silent dropping of patient records during database merging, we developed a deduplication protocol. "
                    "The raw clinical.tsv file from the GDC portal contained 5,207 rows representing treatment history lines for 470 unique patients. "
                    "We collapsed these rows to one unique row per patient case ID. The longest follow-up or survival time was recovered by extracting the maximum of daystodeath (for deceased patients) and daystolastfollowup (for censored patients). "
                    "The patient's vital status was extracted from the primary entry."
                )
                new_dedup_next = (
                    "\n\n[REWRITTEN AGGREGATION PROTOCOL]: To eliminate multi-line record redundancy and prevent sample-selection bias during multi-database merging, a standardized clinical aggregation protocol was established. "
                    "The raw clinical metadata (clinical.tsv) from the GDC portal contained 5,207 longitudinal treatment records corresponding to 470 unique patient cases. "
                    "We implemented a group-by aggregation by unique case ID (cases.case_id), recovering the maximum available follow-up or overall survival duration across recorded days_to_death (for deceased patients) and days_to_last_follow_up (for censored patients), while retaining primary vital status annotations."
                )
                dedup_suffix = (
                    " We constructed a patient-by-CpG methylation matrix using processed tumor tissue beta values for the 214 candidate promoter probes. "
                    "Finally, we matched and joined the deduplicated clinical dataset with the molecular methylation matrix using the unique patient case ID, producing a final validated clinical-methylation cohort of 395 patients."
                )
                
                p_yellow = p.add_run(orig_dedup_yellow)
                p_yellow.font.name = 'Times New Roman'
                p_yellow.font.size = Pt(11)
                p_yellow.font.highlight_color = WD_COLOR_INDEX.YELLOW
                
                p_rewrite = p.add_run(new_dedup_next)
                p_rewrite.font.name = 'Times New Roman'
                p_rewrite.font.size = Pt(11)
                
                p_suf = p.add_run(dedup_suffix)
                p_suf.font.name = 'Times New Roman'
                p_suf.font.size = Pt(11)
                
            elif "To evaluate whether the identified promoter CpG sites" in line:
                p_run1 = p.add_run("To evaluate whether the identified promoter CpG sites can ")
                p_run1.font.name = 'Times New Roman'
                p_run1.font.size = Pt(11)
                del_t = "function as an accurate diagnostic machine learning classifier for melanoma detection, we formulated a binary classification task comparing melanoma tumor tissue methylation profiles against"
                ins_t = "serve as an accurate diagnostic machine-learning classifier for melanoma detection, we formulated a binary classification task that compares melanoma tumor tissue methylation profiles with"
                add_tracked_revision(p, del_t, ins_t, rev_id=30)
                p_run2 = p.add_run(" normal tissue profiles across the 24 significant CpG descriptors. Input features consisted of quantitative promoter CpG methylation beta values (β-values ranging from 0.0 to 1.0).")
                p_run2.font.name = 'Times New Roman'
                p_run2.font.size = Pt(11)
                
            elif "Finally, we performed pathway mapping to classify" in line or "literature curation and pathway mapping" in line:
                p_run1 = p.add_run("For each of the 214 promoter CpG sites, the 395 matched patients were stratified into high-expression and low-expression groups at the median beta value. To prevent statistical artifacts in sparsely populated groups, we excluded sites containing fewer than 5 patients in either group. We modeled cumulative overall survival (OS) distributions using the Kaplan-Meier estimator and performed independent log-rank tests to identify loci significantly associated with patient survival (p < 0.05). Finally, we performed ")
                p_run1.font.name = 'Times New Roman'
                p_run1.font.size = Pt(11)
                add_tracked_revision(p, "literature curation and ", "", rev_id=50)
                p_run2 = p.add_run("pathway mapping to classify the significant genes as oncogenes or tumor suppressors and map their interactions to identify targetable downstream pathway nodes.")
                p_run2.font.name = 'Times New Roman'
                p_run2.font.size = Pt(11)
                
            elif "Across 10-fold cross-validation on the primary training cohort" in line or "learning performance metrics. Write much more text" in line:
                yellow_prompt = "learning performance metrics. Write much more text explaining all these results. Anyway – these are your main results of the paper."
                new_results_next = (
                    " Across 10-fold cross-validation on the primary training cohort (n = 790), all models demonstrated strong predictive performance (Table 5, Figure 5). "
                    "Logistic Regression and Naïve Bayes achieved 92.15% and 92.03% cross-validation accuracy, respectively, with an AUC-ROC of 0.9723. "
                    "The Multi-Layer Perceptron (MLP) and SPegasos SVM followed closely with 91.77% and 91.65% cross-validation accuracy."
                )
                p_yellow = p.add_run(yellow_prompt)
                p_yellow.font.name = 'Times New Roman'
                p_yellow.font.size = Pt(11)
                p_yellow.font.highlight_color = WD_COLOR_INDEX.YELLOW
                
                p_next = p.add_run(new_results_next)
                p_next.font.name = 'Times New Roman'
                p_next.font.size = Pt(11)
                
            else:
                parse_inline_markdown(p, line, default_font_name='Times New Roman', default_font_size=11)
                
        i += 1
        
    print(f"Saving Document to {docx_path}...")
    doc.save(docx_path)
    print("DOCX saved successfully!")

if __name__ == '__main__':
    md_file = r"C:\Users\ethan\.gemini\antigravity\scratch\melanoma_methylation_paper\manuscript.md"
    docx_file = r"C:\Users\ethan\.gemini\antigravity\scratch\melanoma_methylation_paper\manuscript.docx"
    convert_md_to_perfect_docx(md_file, docx_file)
