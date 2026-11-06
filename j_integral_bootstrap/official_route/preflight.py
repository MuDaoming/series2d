"""One-sided rejection of an AB ansatz using the original native ProductionFitter.

An incompatible prefix proves the whole supplied set incompatible. A compatible
prefix does not justify accepting any envelope or omitting final holdouts.
"""
import os
from pathlib import Path
from time import perf_counter


def preflight(ansatz, samples, *, max_points=32, kernel_directory=None):
    """Return status incompatible/compatible/undetermined, never silently fall back.

    compatible means a unique fit to all supplied samples was obtained. Otherwise
    a consistent bounded prefix is undetermined. Only the original fitter's
    InconsistentSampleError is interpreted as mathematical incompatibility.
    Engine failures remain undetermined with an error. No fixed-entry identity
    is imposed here: callers must still run their full joint constrained fit.
    """
    from analytic_bootstrap.upstream.fitting import ProductionFitter
    from analytic_bootstrap.upstream.native_fit import NativeUnavailable
    from analytic_bootstrap.sampling import InconsistentSampleError
    if type(max_points) is not int or max_points<1:raise ValueError('positive max_points required')
    values=list(samples);start=perf_counter();fitter=None;used=0
    directory=kernel_directory or os.environ.get('ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR')
    report=dict(status='undetermined',input_samples=len(values),used_samples=0,
                kernel_directory=str(Path(directory).resolve()) if directory else 'upstream default',
                unique=False)
    old=os.environ.get('ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR')
    try:
        if kernel_directory is not None:os.environ['ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR']=str(Path(kernel_directory).resolve())
        # stream explicitly refuses missing/unsupported native libraries.
        fitter=ProductionFitter(ansatz,backend='stream')
        report['method']=fitter.method
        for sample in values[:max_points]:
            used+=1;result=fitter.add(sample)
            report['rank']=result['rank'];report['dimension']=result['dimension']
        unique=fitter.result()['parameters'] is not None
        report['unique']=unique
        if used==len(values) and unique:report['status']='compatible'
        else:report['reason']='sample prefix consistent; full joint fit still required'
    except InconsistentSampleError as e:
        report.update(status='incompatible',reason=str(e))
    except (NativeUnavailable,OSError,RuntimeError,ArithmeticError,ValueError) as e:
        report.update(status='undetermined',error=repr(e))
    finally:
        if fitter is not None:fitter.close()
        if old is None:os.environ.pop('ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR',None)
        else:os.environ['ANALYTIC_BOOTSTRAP_UPSTREAM_KERNEL_DIR']=old
        report.update(used_samples=used,seconds=perf_counter()-start)
    return report
