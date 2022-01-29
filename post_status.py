import os
import shutil
import urllib.parse
import urllib.request
import logging
import json
import argparse

def post_url(taskID, status, url='http://localhost/task/createSubNodeReport/'):
    url_tmp = urllib.parse.urljoin(url, f'{taskID}/{status}')
    logging.info(url_tmp)
    try:
        request = urllib.request.Request(url_tmp)
        respose = json.loads(urllib.request.urlopen(request).read().decode('utf-8'))
        if respose['status']:
            logging.error(f'post url {url_tmp} ')
            logging.error(f'respose {respose}')
        else:
            logging.info(f'post url {url_tmp} ')
            logging.info(f'respose {respose}')
        return respose['status']
    except Exception as e:
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

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--taskID', required=True, help='task ID')
    parser.add_argument('-n', '--name', default='Fastqc', help='task name')
    args = parser.parse_args()

    post_url(args.taskID, args.name)