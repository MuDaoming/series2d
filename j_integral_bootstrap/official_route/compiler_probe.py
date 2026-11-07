"""Fail early on compute-image compiler/linker/runtime dependency mismatches."""
import argparse,os,subprocess,tempfile,json,socket
from pathlib import Path

def probe(prefix):
    prefix=Path(prefix).resolve()
    compiler=os.environ.get('CXX','g++')
    subprocess.run(['cc','--version'],check=True,stdout=subprocess.DEVNULL)
    version=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]
    with tempfile.TemporaryDirectory(prefix='ab-compiler-probe-') as directory:
        source=Path(directory)/'probe.cpp';binary=Path(directory)/'probe'
        source.write_text('#include <filesystem>\n#include <flint/flint.h>\n#include <iostream>\nint main(){flint_rand_t r;flint_randinit(r);flint_randclear(r);std::cout<<std::filesystem::current_path()<<" OK\\n";}\n')
        subprocess.run([compiler,'-std=c++17',str(source),'-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lflint','-lgmp','-o',str(binary)],check=True)
        result=subprocess.check_output([str(binary)],text=True).strip()
    return dict(host=socket.gethostname(),compiler=compiler,version=version,result=result)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prefix',required=True)
    print(json.dumps(probe(parser.parse_args().prefix)),flush=True)
