import os
import zipfile
import magic
import gzip
from mkdir import mkdir

def file_type(f):
    file_type = magic.from_file(f, mime=True)
    suffix = os.path.splitext(f)[1]
    if file_type == 'text/plain' and suffix in ['.fastq', '.fq']:
        return 'FQ'
    if file_type == 'text/plain' and suffix in ['.fasta', '.fa', 'fna']:
        return 'FA'
    if file_type == 'application/x-gzip' and suffix in ['.gzip', '.gz', '.zip']:
        return 'gz'
    return False

def un_gz(gz_file_name, out_file_name):
    g_file = gzip.GzipFile(gz_file_name)
    open(out_file_name, 'wb').write(g_file.read())
    g_file.close()

def un_zip(file_name, out_file_dir):
    zip_file = zipfile.ZipFile(file_name)
    mkdir(out_file_dir)
    for names in zip_file.namelist():
        zip_file.extract(names, out_file_dir)
    zip_file.close()

