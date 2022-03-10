import time
import argparse

def logtime(*args):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), end='\t')
    for s in args:
        print(f"{s}", end='\t')
    print()


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument('-s', '--str', help="strings to print", default='', nargs='+')
    args = parse.parse_args()
    logtime(' '.join(args.str))