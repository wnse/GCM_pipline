import os

def mkdir(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except Exception as e:
            print(e)
