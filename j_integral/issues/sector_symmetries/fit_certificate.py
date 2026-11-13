"""Reuse the existing finite-field full-DE total-derivative certificate search."""
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[3]
mask=sys.argv[1] if len(sys.argv)>1 else '1110101'
assert len(mask)==7 and set(mask)<=set('01')
p=root/'j_integral/tests/e6/analysis/compare_fi_j.py'
text=p.read_text().replace("OUT=CASE/'analysis';",f"OUT=ROOT/'j_integral/issues/sector_symmetries/e6_{mask}';OUT.mkdir(parents=True,exist_ok=True);")
text=text.replace('mask=(0,1,1,0,1,0,1);','mask='+repr(tuple(map(int,mask)))+';')
exec(compile(text,str(p),'exec'),{'__file__':str(p),'__name__':'__main__'})
