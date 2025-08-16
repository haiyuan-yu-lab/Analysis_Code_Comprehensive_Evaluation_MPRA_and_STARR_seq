# === Standard Library ===
import os
from subprocess import Popen, PIPE, STDOUT, run
import os
import shutil
import subprocess
from pathlib import Path

# === Bioinformatics Tools ===
import pybedtools

def gzip_file(file):
    """
    Compress a file using gzip.
    
    Args:
        file (str): Path to the file to be compressed.
    """
    
    os.system("gzip {}".format(file))

def gunzip_file(file_gz, file=None, keep=True):
    
    """
    Decompress a gzip file.

    Args:
        file_gz (str): Path to the gzipped file.
        file (str, optional): Output path for decompressed content. If not provided, decompresses in place.
        keep (bool): If True and output file is specified, preserves the original .gz file.
    """
    if keep and file:
        os.system(f"gunzip < {file_gz} > {file}")
    else:
        os.system(f"gunzip {file_gz}")

def safe_bedsort(input_file, output_file):
    """
    Sort a BED file by chromosome and start coordinate.

    Args:
        input_file (str): Path to input BED file.
        output_file (str): Path to write the sorted BED file.
    """
    
    os.system(f"sort -k1,1 -k2,2n {input_file} > {output_file}")

def safe_bedsort_optimized(input_file, output_file, tmpdir=None, parallel=None, mem=None):
    """
    Sort a BED file by chr and start (1,1; 2,2n) using GNU sort.
    - Forces C locale for speed/consistency
    - Lets you control temp dir, threads, and memory
    - Raises on error and surfaces stderr
    """
    input_file = str(input_file)
    output_file = str(output_file)
    tmpdir = tmpdir or os.environ.get("TMPDIR") or "/tmp"
    parallel = parallel or os.cpu_count() or 1
    mem = mem or "50%"  # GNU sort memory per process

    # Basic checks
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input not found: {input_file}")
    if os.path.getsize(input_file) == 0:
        # Still write an empty file, then return
        Path(output_file).write_bytes(b"")
        return

    # Ensure output dir exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(tmpdir).mkdir(parents=True, exist_ok=True)

    # Environment: C locale speeds up sort a lot
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    env["TMPDIR"] = tmpdir

    # If your input might be gzipped, handle it transparently
    is_gz = input_file.endswith(".gz")

    # Build command
    if is_gz:
        # Prefer pigz if available (faster), else gzip -dc
        zcat = "pigz" if shutil.which("pigz") else "gzip"
        cat_args = ["-dc"] if zcat != "pigz" else ["-dc"]
        # zcat input | sort ... > output
        sort_cmd = [
            "sort", "-k1,1", "-k2,2n",
            "--parallel", str(parallel),
            "-S", str(mem),
            "-T", tmpdir,
        ]
        # Use a shell pipeline safely via subprocess with pipes
        with open(output_file, "wb") as fout:
            p1 = subprocess.Popen([zcat] + cat_args + [input_file], stdout=subprocess.PIPE, env=env)
            p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=fout, stderr=subprocess.PIPE, env=env)
            p1.stdout.close()  # allow p1 to get SIGPIPE if p2 exits
            _, err = p2.communicate()
            rc = p2.returncode
            if rc != 0:
                raise RuntimeError(f"sort failed (rc={rc}): {err.decode('utf-8', 'ignore')}")
    else:
        # sort input > output (no PIPE for stdout; write directly to file)
        sort_cmd = [
            "sort", "-k1,1", "-k2,2n",
            "--parallel", str(parallel),
            "-S", str(mem),
            "-T", tmpdir,
            input_file,
        ]
        with open(output_file, "wb") as fout:
            res = subprocess.run(sort_cmd, stdout=fout, stderr=subprocess.PIPE, env=env)
        if res.returncode != 0:
            raise RuntimeError(f"sort failed (rc={res.returncode}): {res.stderr.decode('utf-8','ignore')}")


def safe_remove(file):
    """
    Remove a file if it exists.

    Args:
        file (str): Path to the file to be removed.
    """
    
    if os.path.exists(file):
        os.remove(file)

def set_pybedtools_tmp_dir(tmp_dir):
    """
    Set a temporary directory for pybedtools.

    Args:
        tmp_dir (str): Path to temporary directory.
    """
    
    pybedtools.helpers.set_tempdir(tmp_dir)

def set_dir(path):
    """
    Create a directory if it doesn't already exist.

    Args:
        path (str): Directory path to create.
    """
    
    if not os.path.exists(path):
        os.mkdir(path)

def get_row_number(file):
    """
    Count the number of lines (rows) in a file.

    Args:
        file (str): Path to the file.

    Returns:
        int: Number of lines in the file.
    """
    
    results = run(['wc', '-l', file], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    return int(results.stdout.split()[0])

def get_row_number_zipped(file):
    
    """
    Count the number of lines in a file, including gzipped files.

    Args:
        file (str): Path to the (possibly gzipped) file.

    Returns:
        int: Number of lines in the file.
    """
    
    cmd = f"zcat {file} | wc -l" if file.endswith('.gz') else f"wc -l {file}"
    ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    output = ps.communicate()[0].decode("utf-8")
    return int(output.strip().split()[0])
