import os
import re
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import parse_xml, OxmlElement
from docx.oxml.ns import nsdecls, qn

def enable_track_changes(doc):
    """Enable Track Changes in Word document settings."""
    settings = doc.settings.element
    # Check if trackRevisions already exists
    track_rev = settings.find(qn('w:trackRevisions'))
    if track_rev is None:
        track_rev = parse_xml(r'<w:trackRevisions %s/>' % nsdecls('w'))
        settings.append(track_rev)

def add_tracked_paragraph(doc, text_runs, align=WD_ALIGN_PARAGRAPH.LEFT, double_spaced=True, rev_id_start=1, author="Nolan Yu", date_str="2026-08-08T12:00:00Z"):
    """
    text_runs is a list of tuples/dicts:
    ('plain', 'text', bold=False, italic=False)
    ('ins', 'inserted text', bold=False, italic=False)
    ('del', 'deleted text', bold=False, italic=False)
    """
    p = doc.add_paragraph()
    p.alignment = align
    
    p_format = p.paragraph_format
    if double_spaced:
        p_format.line_spacing = 2.0
        p_format.space_after = Pt(0)
    else:
        p_format.line_spacing = 1.15
        p_format.space_after = Pt(6)
        
    rev_id = rev_id_start
    
    for item in text_runs:
        run_type = item[0]
        text = item[1]
        bold = item[2] if len(item) > 2 else False
        italic = item[3] if len(item) > 3 else False
        
        if run_type == 'plain':
            run = p.add_run(text)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(12)
            run.bold = bold
            run.italic = italic
            
        elif run_type == 'ins':
            rPr_xml = '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/>'
            if bold:
                rPr_xml += '<w:b/>'
            if italic:
                rPr_xml += '<w:i/>'
            rPr_xml += '</w:rPr>'
            
            # Escape special XML chars
            escaped_text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            
            ins_xml = (
                f'<w:ins {nsdecls("w")} w:id="{rev_id}" w:author="{author}" w:date="{date_str}">'
                f'<w:r>{rPr_xml}<w:t xml:space="preserve">{escaped_text}</w:t></w:r>'
                f'</w:ins>'
            )
            p._p.append(parse_xml(ins_xml))
            rev_id += 1
            
        elif run_type == 'del':
            rPr_xml = '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/>'
            if bold:
                rPr_xml += '<w:b/>'
            if italic:
                rPr_xml += '<w:i/>'
            rPr_xml += '</w:rPr>'
            
            escaped_text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            
            del_xml = (
                f'<w:del {nsdecls("w")} w:id="{rev_id}" w:author="{author}" w:date="{date_str}">'
                f'<w:r>{rPr_xml}<w:delText xml:space="preserve">{escaped_text}</w:delText></w:r>'
                f'</w:del>'
            )
            p._p.append(parse_xml(del_xml))
            rev_id += 1

    return rev_id

print("Helper functions defined successfully!")
