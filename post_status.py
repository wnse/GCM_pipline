import os
import shutil
import urllib.parse
import urllib.request
from urllib.error import URLError
import logging
import json
import argparse
import datetime
import time

def post_url(taskID, status, url='http://localhost/task/createSubNodeReport/'):
    url_tmp = urllib.parse.urljoin(url, f'{taskID}/{status}')
    logging.info(url_tmp)
    keep_request = True
    try_request_time = 1
    while keep_request:
        try:
            request = urllib.request.Request(url_tmp)
            respose = json.loads(urllib.request.urlopen(request).read().decode('utf-8'))
            if respose['status']:
                logging.error(f'post url {url_tmp} ')
                logging.error(f'respose {respose}')
            else:
                logging.info(f'post url {url_tmp} ')
                logging.info(f'respose {respose}')
            keep_request = False
            return respose['status']
        except URLError as e:
            logging.info(f'{e} try {try_request_time}')
            try_request_time += 1
            if try_request_time > 5:
                keep_request = False
            time.sleep(10)
        except Exception as e:
            keep_request = False
            logging.error(f'post url {e}')

def copy_file(out_file_list, dest_dir):
    for i in out_file_list:
        if os.path.isfile(i):
            try:
                des_file = shutil.copy(i, dest_dir)
                logging.info(des_file)
            except Exception as e:
                logging.error(f'copy file {e}')
        elif os.path.isdir(i):
            try:
                des_file = shutil.copytree(i, os.path.join(dest_dir, os.path.split(i)[1]))
                logging.info(des_file)
            except Exception as e:
                logging.error(e)
        else:
            logging.error(f'copy file not exitst {i}')

def write_status(f, s):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(f, 'a') as h:
        print(f'{now}\t{s}', file=h)

def post_pid(taskID):
    try:
        pid = os.getpid()
        post_url(taskID, str(pid), 'http://localhost/task/getTaskPid/')
    except Exception as e:
        logging.error(f'post_url getTaskPid {e}')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--taskID', required=True, help='task ID')
    parser.add_argument('-n', '--name', default='Fastqc', help='task name')
    args = parser.parse_args()

    post_url(args.taskID, args.name)