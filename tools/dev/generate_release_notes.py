import argparse
import re
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('tag')
args = parser.parse_args()

proc = subprocess.run(["git", "log", "--format=%s", f"{args.tag}.."], capture_output=True, text=True)
data = []
for line in proc.stdout.rstrip().split('\n'):
    m = re.match(r'(.*) \(\#(\d+)\)', line)
    if m is not None:
        data.append(m.groups())
    
for comment, num in sorted(data, key=lambda x: int(x[1])):
    print(f'- {comment} (`#{num} <https://github.com/openmc-dev/openmc/pull/{num}>`_)')
