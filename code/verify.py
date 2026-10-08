"""Fast exact/ball checks of the local certificates. --full reruns completeness.
Run from any directory. This is a verifier using FLINT, not formal proof software.
"""
import argparse,subprocess,sys,json,hashlib,os
from pathlib import Path
from flint import fmpz_poly
from run_paths import RUN_ENV,fresh_run_directory
ROOT=Path(__file__).resolve().parents[1]

def run(name,*args):
    subprocess.run([sys.executable,str(ROOT/'code'/name),*args],check=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
    directory=fresh_run_directory()
    os.environ[RUN_ENV]=str(directory)
    print('Fresh run outputs:',directory,flush=True)
    manifest=ROOT/'MANIFEST.json'
    if manifest.exists():
        m=json.loads(manifest.read_text())
        for path,digest in m['sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
        print('Input/source manifest PASS',len(m['sha256']),'files',flush=True)
    c=json.loads((ROOT/'results'/'candidate32.json').read_text());p=fmpz_poly(c['p'])
    fs=p.factor()[1];assert len(fs)==1 and fs[0][0].degree()==32 and fs[0][1]==1
    r=c['radical'];num=fmpz_poly(r['numerator']);den=fmpz_poly(r['denominator'])
    assert den%p and (num*num-fmpz_poly([c['N'],-c['N']])*den*den)%p==0
    print('Degree-32 x/B certificate PASS',flush=True)
    run('check_maps.py');run('rank16.py');run('check_level1.py')
    run('structure.py');run('check_phi17.py');run('identity32.py')
    if args.full:
        run('reconstruct.py');run('catalogue.py','--bound','64');run('extend32.py')
        run('select_higher.py');run('higher_candidates.py')
    print('PASS: requested checks completed. The report identifies the mathematical inputs.',flush=True)

if __name__=='__main__':main()
